package main

import (
	"encoding/csv"
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/pgs"
	"log"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/nikirill/prs/tools"
)

type Individual struct {
	Genotype []int
	Score    float64
}

type Cohort map[string]*Individual

func NewIndividual() *Individual {
	return &Individual{
		Genotype: make([]int, 0),
		Score:    0.0,
	}
}

func NewCohort(pgs *pgs.PGS) Cohort {
	c := make(Cohort)
	c.Populate(pgs)
	return c
}

func (c Cohort) Populate(pgs *pgs.PGS) {
	filename := fmt.Sprintf("%s.json", pgs.PgsID)
	// If the file doesn't exist, calculate the PRS and save it
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		c.CalculatePRS(pgs)
		c.SaveToDisk(filename)
		if _, err := os.Stat(fmt.Sprintf("%s.scores", pgs.PgsID)); os.IsNotExist(err) {
			c.SaveScores(fmt.Sprintf("%s.scores", pgs.PgsID))
		}
		// Save scores separately for the ease of reading
		return
	}
	// Otherwise, load the data from disk
	c.LoadFromDisk(filename)
}

func (c Cohort) CalculatePRS(pgs *pgs.PGS) {
	for i, locus := range pgs.Loci {
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
			c[individ].Score += float64(allele[0]+allele[1]) * pgs.Weights[i]
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
			if c[sortedInd[j]].Score < c[sortedInd[minIndex]].Score {
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
		writer.Write([]string{ind, fmt.Sprintf("%0.17f", c[ind].Score)})
	}
	writer.Flush()
	return nil
}

func (c Cohort) LoadScores(filename string) error {
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
		if v, err := strconv.ParseFloat(score, 64); err == nil {
			c[name].Score = v
		} else {
			return err
		}
	}
	return nil
}
