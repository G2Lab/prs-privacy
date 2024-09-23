package solver

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strings"
	"sync"

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

func NewCohort(p *pgs.PGS, dataset string) Cohort {
	c := make(Cohort)
	c.Populate(p, dataset)
	return c
}

func (c Cohort) Populate(p *pgs.PGS, dataset string) {
	var err error
	var filename string
	switch dataset {
	case tools.GG:
		filename = fmt.Sprintf("%s/%s.json", params.LocalDataFolder, p.PgsID)
		// If the file doesn't exist, calculate the PRS and save it
		if _, err = os.Stat(filename); os.IsNotExist(err) {
			err = c.RetrieveGenotypes(p, dataset)
			if err != nil {
				log.Fatalf("Error retrieving genotypes and calculating PRS: %v", err)
			}
			c.CalculatePRS(p)
			c.SaveToDisk(filename)
			if _, err = os.Stat(fmt.Sprintf("%s/%s.scores", params.LocalDataFolder, p.PgsID)); os.IsNotExist(err) {
				c.SaveScores(fmt.Sprintf("%s/%s.scores", params.LocalDataFolder, p.PgsID))
			}
			// Save scores separately for the ease of reading
			return
		}
		// Otherwise, load the data from disk
		c.LoadFromDisk(filename)
	case tools.UKB:
		filename = fmt.Sprintf("%s/%s.scores", params.UKBiobankInputFolder, p.PgsID)
		if info, err := os.Stat(filename); os.IsNotExist(err) || info.Size() == 0 {
			err = c.RetrieveGenotypes(p, dataset)
			c.CalculatePRS(p)
			c.SaveScores(filename)
		}
		c.LoadScores(filename, p.Context)
	default:
		log.Fatalf("Unknown dataset: %s", dataset)
	}
}

func (c Cohort) RetrieveGenotypes(p *pgs.PGS, source string) error {
	var err error
	var ok bool
	var output []byte
	var allele []uint8
	var allSamples []string
	var positionSamples, sampleAlleles []string
	var found, locusAdded bool
	var datasets []string
	switch source {
	case tools.GG:
		datasets = []string{tools.GG, tools.RL}
	case tools.UKB:
		datasets = []string{tools.UKB}
	default:
		log.Fatalf("Unknown source: %s", source)
	}
	for _, locus := range p.Loci {
		for _, dataset := range datasets {
			found, locusAdded = false, false
			chr, position := tools.SplitLocus(locus)
			query, args := tools.IndividualSnpsQuery(chr, position, dataset)
			cmd := exec.Command(query, args...)
			output, err = cmd.Output()
			if err != nil {
				fmt.Println("Error executing bcftools command:", err)
				return err
			}
			lines := strings.Split(string(output), "\n")
			for _, line := range lines[:len(lines)-1] {
				positionSamples = strings.Split(line, "*")
				if positionSamples[0] != locus {
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
				switch dataset {
				case tools.GG:
					allSamples = All1000GenomesSamples()
				case tools.RL:
					allSamples = AllRelativeSamples()
				case tools.UKB:
					allSamples = AllUKBiobankSamples()
				default:
					log.Fatalln("Unknown dataset:", dataset)
				}
				for _, indv := range allSamples {
					if _, ok = c[indv]; !ok {
						c[indv] = NewIndividual()
					}
					c[indv].Genotype = append(c[indv].Genotype, ReferenceSNP, ReferenceSNP)
				}
			}
		}
	}
	return nil
}

func (c Cohort) CalculatePRS(p *pgs.PGS) {
	numWorkers := 15
	tasks := make(chan string, 1)
	var wg sync.WaitGroup
	var mu sync.Mutex
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for idv := range tasks {
				score := CalculateDecimalScore(p.Context, c[idv].Genotype, p.Weights, p.EffectAlleles)
				mu.Lock()
				c[idv].Score.Set(score)
				mu.Unlock()
			}
		}()
	}
	for idv := range c {
		tasks <- idv
		//c[idv].Score = CalculateDecimalScore(p.Context, c[idv].Genotype, p.Weights, p.EffectAlleles)
	}
	close(tasks)
	wg.Wait()
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
