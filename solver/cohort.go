package solver

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

type Individual struct {
	Genotype []uint8
	Score    *apd.Decimal
}

type Cohort map[string]*Individual

func NewIndividual() *Individual {
	return &Individual{
		Genotype: make([]uint8, 0),
		Score:    apd.New(0, 0),
	}
}

func NewCohort(p *pgs.PGS) Cohort {
	c := make(Cohort)
	c.Populate(p)
	return c
}

func (c Cohort) Populate(p *pgs.PGS) {
	var err error
	filename := fmt.Sprintf("%s/%s.json", params.DataFolder, p.PgsID)
	// If the file doesn't exist, calculate the PRS and save it
	if _, err = os.Stat(filename); os.IsNotExist(err) {
		err = c.RetrieveGenotypes(p)
		if err != nil {
			log.Fatalf("Error retrieving genotypes and calculating PRS: %v", err)
		}
		c.CalculatePRS(p)
		c.SaveToDisk(filename)
		if _, err = os.Stat(fmt.Sprintf("%s/%s.scores", params.DataFolder, p.PgsID)); os.IsNotExist(err) {
			c.SaveScores(fmt.Sprintf("%s/%s.scores", params.DataFolder, p.PgsID))
		}
		// Save scores separately for the ease of reading
		return
	}
	// Otherwise, load the data from disk
	c.LoadFromDisk(filename)
}

func (c Cohort) RetrieveGenotypes(p *pgs.PGS) error {
	var err error
	var output []byte
	var allele []uint8
	var positionSamples, sampleAlleles []string
	var found, locusAdded bool
	for _, locus := range p.Loci {
		found, locusAdded = false, false
		chr, position := tools.SplitLocus(locus)
		query, args := tools.IndividualSnpsQuery(chr, position)
		cmd := exec.Command(query, args...)
		output, err = cmd.Output()
		if err != nil {
			fmt.Println("Error executing bcftools command:", err)
			return err
		}
		lines := strings.Split(string(output), "\n")
		for _, line := range lines[:len(lines)-1] {
			positionSamples = strings.Split(line, "-")
			if positionSamples[0] != locus {
				//fmt.Printf("Locus %s does not match %s\n", locus, positionSamples[0])
				continue
			}
			found = true
			samples := strings.Split(positionSamples[1], "\t")
			samples = samples[:len(samples)-1]
			for _, sample := range samples {
				sampleAlleles = strings.Split(sample, "=")
				if len(sampleAlleles) != 2 {
					fmt.Printf("Locus %s -- error splitting sample: %s\n", locus, sample)
					return err
				}
				idv := sampleAlleles[0]
				snp := sampleAlleles[1]
				snp, err = tools.NormalizeSnp(snp)
				if err != nil {
					fmt.Printf("Error normalizing %s: %v", snp, err)
					continue
				}
				allele, err = tools.SnpToPair(snp)
				if err != nil {
					fmt.Printf("Error converting SNP %s to an allele: %v", snp, err)
					continue
				}
				if _, ok := c[idv]; !ok {
					c[idv] = NewIndividual()
				}
				if locusAdded {
					// If the same locus is repeated several times with different ref/alt alleles
					sum := c[idv].Genotype[len(c[idv].Genotype)-1] + c[idv].Genotype[len(c[idv].Genotype)-2] + allele[0] + allele[1]
					switch {
					case sum == 1:
						c[idv].Genotype[len(c[idv].Genotype)-2] = 1
						c[idv].Genotype[len(c[idv].Genotype)-1] = 0
					case sum >= 2:
						c[idv].Genotype[len(c[idv].Genotype)-2] = 1
						c[idv].Genotype[len(c[idv].Genotype)-1] = 1
					default:
					}
				} else {
					c[idv].Genotype = append(c[idv].Genotype, allele...)
				}
			}
			locusAdded = true
		}
		if !found {
			// If there is no data, treat as zeros for all individuals
			if len(c) == 0 {
				allSamples := All1000GenomeSamples()
				for _, indv := range allSamples {
					c[indv] = NewIndividual()
				}
			}
			for indv := range c {
				c[indv].Genotype = append(c[indv].Genotype, ReferenceAllele, ReferenceAllele)
			}
		}
	}
	return nil
}

func (c Cohort) CalculatePRS(p *pgs.PGS) {
	for idv := range c {
		c[idv].Score = CalculateDecimalScore(p.Context, c[idv].Genotype, p.Weights, p.EffectAlleles)
		//fmt.Printf("Individual %s: %v\n", idv, c[idv].Genotype)
	}
	// TODO: Storing the linear sum for now to avoid precision issues with division
	//// Normalize the score for each individual by the number of loci
	//divisor := new(apd.Decimal).SetInt64(int64(len(p.Loci) * pgs.Ploidy))
	//for idv := range c {
	//	_, err = ctx.Quo(c[idv].Score, c[idv].Score, divisor)
	//	if err != nil {
	//		log.Println("Error dividing score by the number of loci:", err)
	//		return
	//	}
	//}
}

func (c Cohort) SortByScore() []string {
	sortedInd := make([]string, 0, len(c))
	for ind := range c {
		sortedInd = append(sortedInd, ind)
	}
	for i := 0; i < len(sortedInd)-1; i++ {
		minIndex := i
		for j := i + 1; j < len(sortedInd); j++ {
			if c[sortedInd[j]].Score.Cmp(c[sortedInd[minIndex]].Score) == -1 {
				minIndex = j
			}
		}
		if minIndex != i {
			sortedInd[i], sortedInd[minIndex] = sortedInd[minIndex], sortedInd[i]
		}
	}
	return sortedInd
}

func (c Cohort) SaveToDisk(filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Fatalf("Cannot create json file %v", err)
	}
	defer file.Close()
	encoder := json.NewEncoder(file)
	if err := encoder.Encode(c); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func (c Cohort) LoadFromDisk(filename string) {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Cannot open cohort json file %v", err)
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	if err := decoder.Decode(&c); err != nil {
		log.Fatal("Cannot decode json", err)
	}
}

func (c Cohort) SaveScores(filename string) error {
	sortedInd := c.SortByScore()
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	for _, ind := range sortedInd {
		//writer.Write([]string{ind, fmt.Sprintf(fmt.Sprintf("%%.%df", precision), c[ind].Score)})
		writer.Write([]string{ind, c[ind].Score.String()})
	}
	writer.Flush()
	return nil
}

func (c Cohort) LoadScores(filename string, ctx *apd.Context) error {
	var err error
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err != nil {
			return err
		}
		name := record[0]
		score := record[1]
		c[name] = NewIndividual()
		c[name].Score, _, err = ctx.NewFromString(score)
		if err != nil {
			return err
		}
	}
	return nil
}
