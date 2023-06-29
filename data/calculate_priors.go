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

	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

const (
	numCPUs       = 22
	readBatchSize = 1000
)

var alleles = []string{"0|0", "0|1", "1|0", "1|1"}

func getAllChromosomePositions(chr int) ([]string, error) {
	// Get all the SNP positions for this chromosome
	query, args := tools.AllChrPositionsQuery(strconv.Itoa(chr))
	cmd := exec.Command(query, args...)
	output, err := cmd.Output()
	if err != nil {
		log.Printf("Error querying %d chromosome for positions: %s\n", chr, err)
		return nil, err
	}
	return strings.Split(string(output), "\t"), nil
}

func calculateSNPpriors() {
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
			err = writer.Write(append([]string{"position"}, alleles...))
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
				end = i + readBatchSize
				if end >= len(positions) {
					end = len(positions) - 1
				}
				query, args := tools.RangeSnpValuesQuery(strconv.Itoa(chr), positions[i], positions[end])
				cmd := exec.Command(query, args...)
				output, err := cmd.Output()
				if err != nil {
					log.Printf("Error querying %d:%s-%s: %v\n", chr, positions[i], positions[end], err)
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

func calculatePairwisePriors(filename string) {
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(filename)
	if err != nil {
		log.Println("Error:", err)
		return
	}

	// retrieve the variants for all the individuals at the SNPs of interest
	individuals := make(map[string]map[string]int)
	for _, locus := range p.Loci {
		chr, position := strings.Split(locus, ":")[0], strings.Split(locus, ":")[1]
		query, args := tools.IndividualSnpsQuery(chr, position)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			fmt.Printf("Error querying locus %s: %v\n", locus, err)
			continue
		}
		samples := strings.Split(string(output), "\t")
		for _, sample := range samples {
			fields := strings.Split(sample, "=")
			if len(fields) != 2 {
				continue
			}
			individual := fields[0]
			snp := fields[1]
			if _, exists := individuals[individual]; !exists {
				individuals[individual] = make(map[string]int)
			}
			value, err := tools.SnpToValue(snp)
			if err != nil {
				fmt.Println(err)
				continue
			}
			individuals[individual][locus] = int(value)
		}
	}

	// calculate the pairwise frequencies
	frequency := make([][]float64, len(p.Loci)*len(pgs.GENOTYPES))
	for i := range frequency {
		frequency[i] = make([]float64, len(p.Loci)*len(pgs.GENOTYPES))
	}
	// the result consists of len(pgs.GENOTYPES) x len(p.GENOTYPES) blocks which all add up to len(individuals)
	for i, locus1 := range p.Loci {
		for j, locus2 := range p.Loci {
			if i == j {
				continue
			}
			for _, individual := range individuals {
				frequency[i*len(pgs.GENOTYPES)+individual[locus1]][j*len(pgs.GENOTYPES)+individual[locus2]] += 1
			}
		}
	}

	// normalize and save to a file
	title := strings.Split(filename, "_")[0]
	file, err := os.Create(fmt.Sprintf("data/prior/%s.pairwise", title))
	if err != nil {
		log.Fatalf("Error creating file %s.pairwise: %v\n", title, err)
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	var likelihood float64
	for i := range frequency {
		tmp := make([]string, len(frequency[i]))
		for j := 0; j < len(frequency[i]); j += len(pgs.GENOTYPES) {
			sum := 0.0
			for k := 0; k < len(pgs.GENOTYPES); k++ {
				sum += frequency[i][j+k]
			}
			for k := 0; k < len(pgs.GENOTYPES); k++ {
				if sum != 0 {
					likelihood = frequency[i][j+k] / sum
				} else {
					likelihood = 0
				}
				tmp[j+k] = strconv.FormatFloat(likelihood, 'f', 5, 64)
			}
		}
		err = writer.Write(tmp)
		if err != nil {
			log.Printf("Error writing to file %s.pairwise: %v\n", title, err)
			return
		}
		writer.Flush()
	}
}

func main() {
	prior := flag.Bool("prior", false, "calculate individual SNP priors")
	pairwise := flag.Bool("pairwise", false, "calculate pairwise priors")
	flag.Parse()

	if *prior {
		calculateSNPpriors()
	}
	if *pairwise {
		calculatePairwisePriors("PGS000040_hmPOS_GRCh38.txt")
	}
}
