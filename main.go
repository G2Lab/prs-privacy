package main

import (
	"fmt"
	"log"
	"os/exec"
	"strings"

	"github.com/nikirill/prs/utils"
)

func main() {
	pgs := NewPGS()
	err := pgs.Load("PGS000073_hmPOS_GRCh38.txt")
	if err != nil {
		log.Println("Error:", err)
		return
	}

	individuals := make(map[string]map[string]float64)
	for _, variant := range pgs.Variants {
		chr := variant.GetHmChr()
		position := variant.GetHmPos()
		query, args := utils.IndividualSnpsQuery(chr, position)
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
				fmt.Println("Error splitting sample:", sample)
				continue
			}
			individual := fields[0]
			snp := fields[1]
			value, err := utils.SnpToValue(snp)
			if err != nil {
				fmt.Println(err)
				continue
			}
			if _, ok := individuals[individual]; !ok {
				individuals[individual] = make(map[string]float64)
			}
			key := fmt.Sprintf("%s:%s", chr, position)
			individuals[individual][key] = value
			individuals[individual]["score"] += value * variant.GetWeight()
		}
	}

	sortedInd := make([]string, 0, len(individuals))
	for ind := range individuals {
		sortedInd = append(sortedInd, ind)
	}
	sortByScore(sortedInd, individuals)

	for _, ind := range sortedInd {
		fmt.Println(ind, individuals[ind]["score"])
	}
}

func sortByScore(sortedInd []string, individuals map[string]map[string]float64) {
	for i := 0; i < len(sortedInd)-1; i++ {
		minIndex := i
		for j := i + 1; j < len(sortedInd); j++ {
			if individuals[sortedInd[j]]["score"] < individuals[sortedInd[minIndex]]["score"] {
				minIndex = j
			}
		}
		if minIndex != i {
			sortedInd[i], sortedInd[minIndex] = sortedInd[minIndex], sortedInd[i]
		}
	}
}
