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
	"gonum.org/v1/gonum/stat"
)

const (
	numCPUs       = 22
	readBatchSize = 10000
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
	positions := strings.Split(string(output), "\t")
	// Every returned element is position+Tab so the last element after split by tab is an empty space
	return positions[:len(positions)-1], nil
}

func calculatePriors() {
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
						if strings.Contains(sample, ".") {
							continue
						}
						sample, err = tools.NormalizeSnp(sample)
						if err != nil {
							log.Printf("%v: %s, snp %s\n", err, sample, positions[i+k])
							continue
						}
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

func calculateCorrelation(filename string) {
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(filename)
	if err != nil {
		log.Println("Error:", err)
		return
	}

	// retrieve the variants for all the individuals at the SNPs of interest
	var snp string
	genotypes := make([][]float64, len(p.Loci))
	for i, locus := range p.Loci {
		genotypes[i] = make([]float64, 0)
		chr, position := strings.Split(locus, ":")[0], strings.Split(locus, ":")[1]
		samples, err := tools.GetSnpsAtPosition(chr, position)
		if err != nil {
			log.Printf("Error querying for SNPs at %s: %v", locus, err)
			continue
		}
		for _, sample := range samples {
			snp, err = tools.NormalizeSnp(sample)
			if err != nil {
				log.Printf("%v: %s, snp %s\n", err, sample, locus)
				genotypes[i] = append(genotypes[i], []float64{0, 0}...)
				continue
			}
			diploids, err := tools.SnpToPair(snp)
			if err != nil {
				log.Printf("%v: %s, snp %s\n", err, sample, locus)
				genotypes[i] = append(genotypes[i], []float64{0, 0}...)
				continue
			}
			genotypes[i] = append(genotypes[i], diploids...)
		}
	}

	correlations := make([][]float64, len(p.Loci))
	for i := range p.Loci {
		correlations[i] = make([]float64, len(p.Loci))
		for j := range p.Loci {
			correlations[i][j] = stat.Correlation(genotypes[i], genotypes[j], nil)
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
	for i := range correlations {
		tmp := make([]string, len(correlations[i]))
		for j := 0; j < len(tmp); j++ {
			tmp[j] = strconv.FormatFloat(correlations[i][j], 'f', 5, 64)
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
		calculatePriors()
	}
	if *pairwise {
		calculateCorrelation("PGS000040_hmPOS_GRCh38.txt")
	}
}
