package main

import (
	"encoding/csv"
	"fmt"
	"github.com/nikirill/prs/pgs"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/nikirill/prs/tools"
)

type Cohort map[string]map[string]float64

func NewCohort() Cohort {
	c := make(Cohort)
	return c
}

func (c Cohort) CalculatePRS(pgs *pgs.PGS) {
	for _, variant := range pgs.Variants {
		chr := variant.GetHmChr()
		position := variant.GetHmPos()
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
			individual := fields[0]
			snp := fields[1]
			value, err := tools.SnpToValue(snp)
			if err != nil {
				fmt.Println(err)
				continue
			}
			if _, ok := c[individual]; !ok {
				c[individual] = make(map[string]float64)
			}
			key := fmt.Sprintf("%s:%s", chr, position)
			c[individual][key] = value
			c[individual][SCORE] += value * variant.GetWeight()
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
			if c[sortedInd[j]][SCORE] < c[sortedInd[minIndex]][SCORE] {
				minIndex = j
			}
		}
		if minIndex != i {
			sortedInd[i], sortedInd[minIndex] = sortedInd[minIndex], sortedInd[i]
		}
	}
	return sortedInd
}

func (c Cohort) SaveScores(filename string) error {
	sortedInd := c.SortByScore()
	file, err := os.Create(filename + ".scores")
	if err != nil {
		return err
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	for _, ind := range sortedInd {
		writer.Write([]string{ind, fmt.Sprintf("%0.17f", c[ind][SCORE])})
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
		individual := record[0]
		score := record[1]
		c[individual] = make(map[string]float64)
		if v, err := strconv.ParseFloat(score, 64); err == nil {
			c[individual][SCORE] = v
		} else {
			return err
		}
	}
	return nil
}
