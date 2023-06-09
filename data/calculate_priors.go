package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"sync"

	"github.com/nikirill/prs/tools"
)

const (
	numCPUs       = 8
	readBatchSize = 500
)

func main() {
	// first we send a command without -r and obtain all the positions in the file
	// then we query row by row and compute frequencies for each allele
	alleles := []string{"0|0", "0|1", "1|0", "1|1"}
	nRoutines := 0
	wg := &sync.WaitGroup{}
	//for chromosome := 1; chromosome <= 22; chromosome++ {
	for chromosome := 22; chromosome <= 22; chromosome++ {
		wg.Add(1)
		go func(chr int) {
			// Prepare a file to write results
			file, err := os.Create(fmt.Sprintf("stats/data/chr%d.csv", chr))
			if err != nil {
				log.Fatalf("Error creating file chr%d.csv: %s\n", chr, err)
			}
			writer := csv.NewWriter(file)
			// Write header
			err = writer.Write(append([]string{"position"}, alleles...))
			if err != nil {
				log.Fatalf("Error writing header to file chr%d.csv: %s\n", chr, err)
			}
			// Get all the SNP positions for this chromosome
			query, args := tools.AllChrPositionsQuery(strconv.Itoa(chr))
			cmd := exec.Command(query, args...)
			output, err := cmd.Output()
			if err != nil {
				log.Printf("Error querying %d chromosome: %s\n", chr, err)
				return
			}
			positions := strings.Split(string(output), "\t")
			end := 0
			for i := 0; i < len(positions); i += readBatchSize {
				end = i + readBatchSize
				if end >= len(positions) {
					end = len(positions) - 1
				}
				query, args := tools.RangeSnpValuesQuery(strconv.Itoa(chr), positions[i], positions[end])
				cmd := exec.Command(query, args...)
				output, err := cmd.Output()
				if err != nil {
					log.Printf("Error querying %d:%s-%s: %s\n", chr, positions[i], positions[end], err)
					return
				}
				lines := strings.Split(string(output), "\n")
				for k, line := range lines[:len(lines)-1] {
					samples := strings.Split(line, "\t")
					counts := make(map[string]int)
					for _, sample := range samples[:len(samples)-1] {
						if _, exists := counts[sample]; exists {
							counts[sample] += 1
						} else {
							counts[sample] = 1
						}
					}
					total := 0
					for _, allele := range alleles {
						if count, exists := counts[allele]; exists {
							total += count
						}
					}
					freq := make([]string, len(alleles))
					for j, allele := range alleles {
						freq[j] = fmt.Sprintf("%.5f", float64(counts[allele])/float64(total))
					}
					err = writer.Write(append([]string{fmt.Sprintf("%s", positions[i+k])}, freq...))
					if err != nil {
						log.Printf("Error writing to file chr%d.csv: %s\n", chr, err)
						return
					}
				}
			}
			writer.Flush()
			err = file.Close()
			if err != nil {
				log.Printf("Error closing file chr%d.csv: %s\n", chr, err)
			}
			wg.Done()
		}(chromosome)
		if nRoutines++; nRoutines == numCPUs {
			wg.Wait()
			nRoutines = 0
		}
	}
	wg.Wait()
}
