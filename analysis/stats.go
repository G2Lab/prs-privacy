package main

import (
	"bufio"
	"flag"
	"fmt"
	"github.com/montanaflynn/stats"
	"github.com/nikirill/prs/pgs"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

const catalog = "catalog"

func listFiles(folderPath string) ([]string, error) {
	dir, err := os.Open(folderPath)
	if err != nil {
		return nil, err
	}
	defer dir.Close()

	fileInfos, err := dir.Readdir(0)
	if err != nil {
		return nil, err
	}

	var fileNames []string
	for _, fileInfo := range fileInfos {
		if fileInfo.IsDir() || fileInfo.Name()[:3] != "PGS" {
			continue
		}
		fileNames = append(fileNames, fileInfo.Name())
	}

	// Sort the file names
	sort.Strings(fileNames)

	return fileNames, nil
}

func numVariantsStats(fileNames []string) {
	var err error
	numVariants := make([]float64, 0)
	for _, fileName := range fileNames {
		fmt.Println(fileName)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalog, fileName))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		numVariants = append(numVariants, float64(p.VariantCount))
	}
	var mean, median float64
	median, err = stats.Median(numVariants)
	if err != nil {
		log.Println("Median error:", err)
	}
	mean, err = stats.Mean(numVariants)
	if err != nil {
		log.Println("Mean error:", err)
	}
	fmt.Printf("Median and mean number of SNPs per PGS: %d, %d\n", int(median), int(mean))
}

func lociOverlapStats(fileNames []string) {
	sizeLimit := 200
	var pgsID, chr, pos, locus string
	var headerInProgress bool
	var variantCount int
	var fieldnames []string
	overlaps := make(map[string][]string)
	variantNum := make(map[string]int)
	pgsLoci := make(map[string][]string)
filesLoop:
	for _, fileName := range fileNames {
		//fmt.Printf("%s ", fileName)
		file, err := os.Open(filepath.Join(catalog, fileName))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		defer file.Close()

		headerInProgress = true
		scanner := bufio.NewScanner(file)
	scannerLoop:
		for scanner.Scan() {
			line := scanner.Text()

			// If it is the header
			if strings.HasPrefix(line, "#") {
				fields := strings.SplitN(line[1:], "=", 2)
				switch strings.ToLower(fields[0]) {
				case "pgs_id":
					pgsID = fields[1]
				case "variants_number":
					if variantCount, err = strconv.Atoi(fields[1]); err == nil {
						// We skip PGS with too many variants
						if variantCount > sizeLimit {
							continue filesLoop
						}
						variantNum[pgsID] = variantCount
					} else {
						log.Printf("Error parsing variants number %s: %s", fields[1], err)
					}
				}
				continue scannerLoop
			}
			if headerInProgress {
				fieldnames = strings.Split(line, "\t")
				headerInProgress = false
				continue scannerLoop
			}

			values := strings.Split(line, "\t")
			for i, value := range values {
				switch fieldnames[i] {
				case "hm_chr":
					if strings.Contains(value, "_") {
						value = strings.Split(value, "_")[0]
						if value == "Un" {
							for j := 0; j < len(fieldnames); j++ {
								if fieldnames[j] == "chr_name" {
									value = values[j]
									break
								}
							}
						}
					}
					// If there is no mapping, we skip the variant
					if len(value) == 0 {
						//fmt.Printf("No mapping for variant %s\n", values[0])
						continue scannerLoop
					}
					chr = value
				case "hm_pos":
					pos = value
				default:
					continue
				}
			}
			locus = fmt.Sprintf("%s:%s", chr, pos)
			if _, ok := overlaps[locus]; !ok {
				overlaps[locus] = make([]string, 0)
			}
			overlaps[locus] = append(overlaps[locus], pgsID)
			if _, ok := pgsLoci[pgsID]; !ok {
				pgsLoci[pgsID] = make([]string, 0)
			}
			pgsLoci[pgsID] = append(pgsLoci[pgsID], locus)
		}
		if err = scanner.Err(); err != nil {
			log.Println("Scanner error:", err)
			return
		}
	}

	// We count the number of overlaps with more than one PGS
	var small bool
	var loci []string
	candidates := make(map[string][]string)
	for l, pgsIDs := range overlaps {
		if len(pgsIDs) < 2 || l[:1] == "X" || l[:1] == "Y" {
			continue
		}
		small = false
		for _, pgsID = range pgsIDs {
			if variantNum[pgsID] < 40 {
				small = true
				break
			}
		}
		if !small {
			continue
		}
		matches := 1
		for _, pgsID = range pgsIDs {
			if _, ok := candidates[pgsID]; ok || variantNum[pgsID] < 70 {
				continue
			}
			matchingLoci := make([]string, 0)
			loci = pgsLoci[pgsID]
			for _, locus = range loci {
				if locus == l || locus[:1] == "X" || locus[:1] == "Y" {
					continue
				}
				others := overlaps[locus]
				if len(others) < 2 {
					continue
				}
				for _, other := range others {
					if other == pgsID || variantNum[other] > 40 {
						continue
					}
					matches++
					matchingLoci = append(matchingLoci, locus+"-"+other+"("+strconv.Itoa(variantNum[other])+")")
				}
			}
			if matches > 5 {
				candidates[pgsID] = matchingLoci
			}
		}
		//fmt.Printf("\n%s: ", l)
		//for _, pgsID = range pgsIDs {
		//	fmt.Printf("%s(%d) ", pgsID, variantNum[pgsID])
		//}
	}
	for pgsID = range candidates {
		fmt.Printf("%s (%d): %s\n", pgsID, variantNum[pgsID], strings.Join(candidates[pgsID], ","))
	}
}

func main() {
	expr := flag.String("expr", "", "Experiment type")
	flag.Parse()
	filenames, err := listFiles(catalog)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	switch *expr {
	case "num":
		numVariantsStats(filenames)
	case "overlap":
		lociOverlapStats(filenames)
	}
}
