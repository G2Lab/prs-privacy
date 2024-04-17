package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/montanaflynn/stats"
	"github.com/nikirill/prs/pgs"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"regexp"
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
		numVariants = append(numVariants, float64(p.NumVariants))
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

func getPGSWithFewerVariants(limit int) {
	file, err := os.Open("catalog/pgs_all_metadata_scores.csv")
	if err != nil {
		log.Println("Error opening catalog metadata file:", err)
		return
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		fmt.Println("Error reading header:", err)
		return
	}
	numVariantColumn, pgsIdColumn := -1, -1
	for i, field := range header {
		if field == "Number of Variants" {
			numVariantColumn = i
		}
		if field == "Polygenic Score (PGS) ID" {
			pgsIdColumn = i
		}
	}
	ids := make([]string, 0)
	total := 0
	for {
		record, err := reader.Read()
		if err != nil {
			break
		}
		if numVariants, err := strconv.Atoi(record[numVariantColumn]); err == nil && numVariants < limit {
			ids = append(ids, record[pgsIdColumn])
		}
		total++
	}
	catalogFolder := "catalog"
	hasXYchromosomes := make([]string, 0)
	invalidLoci := make([]string, 0)
	tooFewVariants := make([]string, 0)
	impreciseWeights := make([]string, 0)
	filteredIds := make([]string, 0)
	idsToNumVariants := make(map[string]int)
	traits := make(map[string]struct{})
	publications := make(map[string]struct{})
	uniqueLoci := make(map[string]struct{})
	lociTotal := 0
	lociToPgp := make(map[string][]string)
	lociToPgs := make(map[string]map[string]bool)
	pgsToPgp := make(map[string]string)
idLoop:
	for _, id := range ids {
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalogFolder, id+"_hmPOS_GRCh38.txt"))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		if p.NumVariants < 2 {
			tooFewVariants = append(tooFewVariants, id)
			continue idLoop
		}
		for _, locus := range p.Loci {
			if strings.HasPrefix(locus, "X:") || strings.HasPrefix(locus, "Y:") {
				hasXYchromosomes = append(hasXYchromosomes, id)
				continue idLoop
			}
			if !isNumber(strings.Split(locus, ":")[1]) || !isNumberInRange(strings.Split(locus, ":")[0], 1, 22) {
				invalidLoci = append(invalidLoci, id)
				continue idLoop
			}
		}
		if p.WeightPrecision < uint32(math.Floor(math.Log2(float64(p.NumVariants)))) {
			impreciseWeights = append(impreciseWeights, id)
			continue idLoop
		}
		filteredIds = append(filteredIds, id)
		idsToNumVariants[id] = p.NumVariants
		pgsToPgp[id] = p.PgpID
		traits[p.TraitEFO] = struct{}{}
		publications[p.PgpID] = struct{}{}
		lociTotal += len(p.Loci)
		for _, locus := range p.Loci {
			uniqueLoci[locus] = struct{}{}
			if _, ok := lociToPgp[locus]; !ok {
				lociToPgp[locus] = make([]string, 0)
			}
			lociToPgp[locus] = append(lociToPgp[locus], p.PgpID)
			if _, ok := lociToPgs[locus]; !ok {
				lociToPgs[locus] = make(map[string]bool)
			}
			lociToPgs[locus][id] = true
		}
	}
	//fmt.Println(len(ids), "PGS with less than", limit, "variants")
	fmt.Println("Total PGS:", total)
	fmt.Println("PGS with less than", limit, "variants:", len(ids))
	fmt.Println("PGS with X or Y chromosomes:", len(hasXYchromosomes), hasXYchromosomes)
	fmt.Println("PGS with less than 2 variants:", len(tooFewVariants), tooFewVariants)
	fmt.Println("PGS with imprecise weights:", len(impreciseWeights), impreciseWeights)
	fmt.Println("PGS with invalid loci:", len(invalidLoci), invalidLoci)
	fmt.Println("Filtered PGS:", len(filteredIds))
	fmt.Println("Unique traits:", len(traits))
	fmt.Println("Unique publications:", len(publications))
	fmt.Printf("Total loci: %d, unique %d\n", lociTotal, len(uniqueLoci))

	lociPgpFile, err := os.OpenFile("results/loci_pgp_coverage.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening loci file:", err)
		return
	}
	defer lociPgpFile.Close()
	encoder := json.NewEncoder(lociPgpFile)
	err = encoder.Encode(lociToPgp)
	if err != nil {
		log.Println("Error encoding pgp loci:", err)
	}

	lociPgsFile, err := os.OpenFile("results/loci_pgs_coverage.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening loci file:", err)
		return
	}
	defer lociPgsFile.Close()
	encoder = json.NewEncoder(lociPgsFile)
	err = encoder.Encode(lociToPgs)
	if err != nil {
		log.Println("Error encoding pgs loci:", err)
	}

	idsFile, err := os.OpenFile("results/filtered_pgs.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening ids file:", err)
		return
	}
	defer idsFile.Close()
	encoder = json.NewEncoder(idsFile)
	err = encoder.Encode(idsToNumVariants)
	if err != nil {
		log.Println("Error encoding filtered ids:", err)
	}

	pgsToPgpFile, err := os.OpenFile("results/pgs_pgp.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening ids file:", err)
		return
	}
	defer idsFile.Close()
	encoder = json.NewEncoder(pgsToPgpFile)
	err = encoder.Encode(pgsToPgp)
	if err != nil {
		log.Println("Error encoding pgs to pgp:", err)
	}
}

func copyFilteredPGSFiles() {
	file, err := os.Open("results/filtered_pgs.json")
	if err != nil {
		log.Println("Error opening filtered ids file:", err)
		return
	}
	defer file.Close()

	decoder := json.NewDecoder(file)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding filtered ids:", err)
		return
	}
	for id := range idsToNumVariants {
		src := filepath.Join("catalog", id+"_hmPOS_GRCh38.txt")
		dst := filepath.Join("filtered", id+"_hmPOS_GRCh38.txt")
		err = copyFile(src, dst)
		if err != nil {
			log.Println("Error copying file:", err)
		}
	}

}

func validateBases() {
	file, err := os.Open("results/filtered_pgs.json")
	if err != nil {
		log.Println("Error opening filtered ids file:", err)
		return
	}
	defer file.Close()

	decoder := json.NewDecoder(file)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding filtered ids:", err)
		return
	}
	discarded := make([]string, 0)
	validated := make(map[string]int)
	for id, num := range idsToNumVariants {
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join("filtered", id+"_hmPOS_GRCh38.txt"))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		err = p.LoadStats()
		if err != nil {
			discarded = append(discarded, p.PgsID)
			continue
		}
		err = copyFile(filepath.Join("filtered", id+"_hmPOS_GRCh38.txt"),
			filepath.Join("inputs", id+"_hmPOS_GRCh38.txt"))
		if err != nil {
			log.Println("Error copying file:", err)
		}
		validated[p.PgsID] = num
	}
	fmt.Printf("Discarded %d: %v\n", len(discarded), discarded)
	fmt.Printf("Validated %d\n", len(validated))

	validatedFile, err := os.OpenFile("results/validated_pgs.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening validated file:", err)
		return
	}
	defer validatedFile.Close()
	encoder := json.NewEncoder(validatedFile)
	err = encoder.Encode(validated)
	if err != nil {
		log.Println("Error encoding validated ids:", err)
	}
}

func makeLociIndex() {
	file, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()

	decoder := json.NewDecoder(file)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}
	loci := make(map[string][]string)
	for id := range idsToNumVariants {
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join("inputs", id+"_hmPOS_GRCh38.txt"))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		for _, locus := range p.Loci {
			if _, ok := loci[locus]; !ok {
				loci[locus] = make([]string, 0)
			}
			loci[locus] = append(loci[locus], p.PgsID)
		}
	}

	// Sort each entry in the map by the number of variants
	for locus := range loci {
		sort.Slice(loci[locus], func(i, j int) bool {
			return idsToNumVariants[loci[locus][i]] < idsToNumVariants[loci[locus][j]]
		})
	}

	validatedLoci, err := os.OpenFile("results/validated_loci.json", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0777)
	if err != nil {
		log.Println("Error opening validated file:", err)
		return
	}
	defer validatedLoci.Close()
	encoder := json.NewEncoder(validatedLoci)
	err = encoder.Encode(loci)
	if err != nil {
		log.Println("Error encoding validated loci:", err)
	}

}

func copyFile(sourceFile, destFile string) error {
	// Open the source file for reading
	source, err := os.Open(sourceFile)
	if err != nil {
		return err
	}
	defer source.Close()

	// Create the destination file for writing
	destination, err := os.Create(destFile)
	if err != nil {
		return err
	}
	defer destination.Close()

	// Copy the contents from the source file to the destination file
	_, err = io.Copy(destination, source)
	if err != nil {
		return err
	}

	return nil
}

// isNumber checks if a string represents a valid number using a regular expression
func isNumber(str string) bool {
	numberRegex := regexp.MustCompile(`^[\-+]?(\d+(\.\d+)?|\.\d+)$`)
	return numberRegex.MatchString(str)
}

// isNumberInRange checks if a string represents a valid number and falls within the specified range
func isNumberInRange(str string, min, max int) bool {
	num, err := strconv.Atoi(str)
	if err != nil {
		return false
	}
	return num >= min && num <= max
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

func removeFiles() {
	file, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()

	decoder := json.NewDecoder(file)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}
	for id := range idsToNumVariants {
		err = copyFile(filepath.Join("inputs", id+"_hmPOS_GRCh38_stats.txt"), filepath.Join("inputs", id+"_hmPOS_GRCh38.txt"))
		if err != nil {
			log.Println("Error copying file:", err)
		}
		err = os.Remove(filepath.Join("inputs", id+"_hmPOS_GRCh38_stats.txt"))
		if err != nil {
			log.Println("Error removing file:", err)
		}
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
	case "limit":
		getPGSWithFewerVariants(50)
	case "copy":
		copyFilteredPGSFiles()
	case "validate":
		validateBases()
	case "loci":
		makeLociIndex()
	}
}
