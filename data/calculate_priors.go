package main

import (
	"encoding/csv"
	"encoding/json"
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
	readBatchSize = 50000
)

var alleleVariants = []string{"0", "1"}

const (
	MajorAllele = "0"
	MinorAllele = "1"
)

func getAllChromosomePositions(chr int) ([]string, error) {
	// Get all the SNP positions for this chromosome
	query, args := tools.AllChrPositionsQuery(strconv.Itoa(chr), tools.GG)
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

// Calculate Minor Allele Frequency (effect allele) for each SNP
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
			err = writer.Write([]string{"position", "minor allele likelihood"})
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
				query, args := tools.RangeSnpValuesQuery(strconv.Itoa(chr), positions[i], positions[end], tools.GG)
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
						float64(counts[MinorAllele])/float64(total))))
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
			fmt.Printf("Chromosome %d complete\n", chr)
			wg.Done()
		}(chromosome)
		if nRoutines++; nRoutines == numCPUs {
			wg.Wait()
			nRoutines = 0
		}
	}
	wg.Wait()
}

func readSamplePopulations() {
	fileName := "data/igsr_samples.tsv"

	// Open the TSV inputFile
	inputFile, err := os.Open(fileName)
	if err != nil {
		fmt.Println("Error:", err)
		return
	}
	defer inputFile.Close()

	// Create a CSV reader with tab as the delimiter
	reader := csv.NewReader(inputFile)
	reader.Comma = '\t'

	// Read the header to get column names
	header, err := reader.Read()
	if err != nil {
		fmt.Println("Error reading header:", err)
		return
	}

	sampleIdx := indexOf(header, "Sample name")
	superPopIdx := indexOf(header, "Superpopulation code")

	if sampleIdx == -1 || superPopIdx == -1 {
		fmt.Println("Required fields not found in the header.")
		return
	}

	fullData := make(map[string]string)
	var record []string
	// Read the remaining records and extract fields
	for {
		record, err = reader.Read()
		if err != nil {
			break // End of inputFile
		}
		fullData[record[sampleIdx]] = record[superPopIdx]
	}

	//// Check that we have all the sample info for the 1000 Genomes dataset
	//sampleFile, err := os.Open("data/1000genome-samples.csv")
	//if err != nil {
	//	fmt.Println("Error:", err)
	//	return
	//}
	//defer sampleFile.Close()
	//
	//// Create a CSV reader with tab as the delimiter
	//reader = csv.NewReader(sampleFile)
	//
	//// Read the header to get column names
	//sampleList, err := reader.Read()
	//if err != nil {
	//	fmt.Println("Error reading sample list:", err)
	//	return
	//}
	//reducedData := make(map[string]string)
	//for _, sample := range sampleList {
	//	if _, exists := fullData[sample]; !exists {
	//		fmt.Printf("Sample %s not found in the superpopulation fullData\n", sample)
	//		continue
	//	}
	//	reducedData[sample] = fullData[sample]
	//}

	// Save the results to a JSON file
	outputFile, err := os.Create("data/superpopulations.json")
	if err != nil {
		fmt.Println("Error creating JSON outputFile:", err)
		return
	}
	defer outputFile.Close()

	// Encode the fullData as JSON and write to the outputFile
	encoder := json.NewEncoder(outputFile)
	err = encoder.Encode(fullData)
	if err != nil {
		fmt.Println("Error encoding fullData to JSON:", err)
	}
}

// indexOf finds the index of a value in a slice, or returns -1 if not found
func indexOf(slice []string, value string) int {
	for i, v := range slice {
		if v == value {
			return i
		}
	}
	return -1
}

func main() {
	maf := flag.Bool("maf", false, "calculate MAF")
	readSP := flag.Bool("populations", false, "read sample populations")
	flag.Parse()

	if *maf {
		calculateMAF()
	}
	if *readSP {
		readSamplePopulations()
	}
}
