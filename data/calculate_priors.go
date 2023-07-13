package main

import (
	"encoding/csv"
	"flag"
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
	numCPUs       = 22
	readBatchSize = 10000
)

var alleleVariants = []string{"0", "1"}

const MAJOR_ALLELE = "0"

func getAllChromosomePositions(chr int) ([]string, error) {
	// Get all the SNP positions for this chromosome
	query, args := tools.AllChrPositionsQuery(strconv.Itoa(chr))
	cmd := exec.Command(query, args...)
	output, err := cmd.Output()
	if err != nil {
		log.Printf("Error querying %d chromosome for positions: %s\n", chr, err)
		return nil, err
	}
	positions := strings.Split(string(output), "\t")
	// Every returned element is position+Tab so the last element after split by tab is an empty space
	return positions[:len(positions)-1], nil
}

func calculateMAF() {
	// first we send a command without -r and obtain all the positions in the file
	// then we query row by row and compute frequencies for each allele
	nRoutines := 0
	wg := &sync.WaitGroup{}
	for chromosome := 1; chromosome <= 22; chromosome++ {
		wg.Add(1)
		go func(chr int) {
			// Prepare a file to write results
			file, err := os.Create(fmt.Sprintf("data/prior/chr%d.csv", chr))
			if err != nil {
				log.Fatalf("Error creating file chr%d.csv: %v\n", chr, err)
			}
			writer := csv.NewWriter(file)
			// Write header
			err = writer.Write([]string{"position", "major allele likelihood"})
			if err != nil {
				log.Fatalf("Error writing header to file chr%d.csv: %v\n", chr, err)
			}
			// Get all the SNP positions for this chromosome
			positions, err := getAllChromosomePositions(chr)
			if err != nil {
				return
			}
			end := 0
			for i := 0; i < len(positions); i += readBatchSize {
				end = i + readBatchSize - 1
				if end >= len(positions) {
					end = len(positions) - 1
				}
				query, args := tools.RangeSnpValuesQuery(strconv.Itoa(chr), positions[i], positions[end])
				cmd := exec.Command(query, args...)
				batch, err := cmd.Output()
				if err != nil {
					log.Printf("Error querying %d:%s-%s: %v\n", chr, positions[i], positions[end], err)
					return
				}
				lines := strings.Split(string(batch), "\n")
				for k, line := range lines[:len(lines)-1] {
					samples := strings.Split(line, "\t")
					counts := make(map[string]int)
					for _, sample := range samples[:len(samples)-1] {
						alleles := strings.Split(sample, "|")
						for _, allele := range alleles {
							normalized, err := tools.NormalizeAllele(allele)
							if err != nil {
								if normalized != "." {
									log.Printf("%v: %s, snp %s\n", err, sample, positions[i+k])
								}
								continue
							}
							if _, exists := counts[normalized]; exists {
								counts[normalized] += 1
							} else {
								counts[normalized] = 1
							}
						}
					}
					total := 0
					for _, allele := range alleleVariants {
						if count, exists := counts[allele]; exists {
							total += count
						}
					}
					err = writer.Write(append([]string{fmt.Sprintf("%s", positions[i+k])}, fmt.Sprintf("%.5f",
						float64(counts[MAJOR_ALLELE])/float64(total))))
					if err != nil {
						log.Printf("Error writing to file chr%d.csv: %v\n", chr, err)
						return
					}
				}
			}
			writer.Flush()
			err = file.Close()
			if err != nil {
				log.Printf("Error closing file chr%d.csv: %v\n", chr, err)
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

func main() {
	maf := flag.Bool("maf", false, "calculate MAF")
	flag.Parse()

	if *maf {
		calculateMAF()
	}
}
