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
	filename := fmt.Sprintf("%s/%s.json", params.DataFolder, p.PgsID)
	// If the file doesn't exist, calculate the PRS and save it
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		c.CalculatePRS(p)
		c.SaveToDisk(filename)
		if _, err := os.Stat(fmt.Sprintf("%s/%s.scores", params.DataFolder, p.PgsID)); os.IsNotExist(err) {
			c.SaveScores(fmt.Sprintf("%s/%s.scores", params.DataFolder, p.PgsID))
		}
		// Save scores separately for the ease of reading
		return
	}
	// Otherwise, load the data from disk
	c.LoadFromDisk(filename)
}

func (c Cohort) CalculatePRS(p *pgs.PGS) {
	ctx := p.Context
	var k uint8
	for i, locus := range p.Loci {
		chr, position := tools.SplitLocus(locus)
		query, args := tools.IndividualSnpsQuery(chr, position)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			fmt.Println("Error executing bcftools command:", err)
			continue
		}

		samples := strings.Split(string(output), "\t")
		for _, sample := range samples {
			fields := strings.Split(sample, "=")
			if len(fields) != 2 {
				//fmt.Println("Error splitting sampleFromMap:", sampleFromMap)
				continue
			}
			individ := fields[0]
			snp := fields[1]
			snp, err = tools.NormalizeSnp(snp)
			if err != nil {
				fmt.Printf("Error normalizing %s: %v", snp, err)
				continue
			}
			allele, err := tools.SnpToPair(snp)
			if err != nil {
				fmt.Printf("Error converting SNP %s to an allele: %v", snp, err)
				continue
			}
			if _, ok := c[individ]; !ok {
				c[individ] = NewIndividual()
			}
			c[individ].Genotype = append(c[individ].Genotype, allele...)
			for k = 0; k < allele[0]+allele[1]; k++ {
				_, err = ctx.Add(c[individ].Score, c[individ].Score, p.Weights[i])
				if err != nil {
					log.Println("Error adding to score:", err)
					return
				}
			}
		}
	}
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
