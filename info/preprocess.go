package main

import (
	"bytes"
	"compress/gzip"
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"net/http"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/nikirill/prs/data"
	"github.com/nikirill/prs/pgs"
)

func getPGSWithFewerVariants(limit int) {
	ids, err := fewerVariantsPGS(0, limit)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	catalogFolder := "catalog"
	hasXYchromosomes := make([]string, 0)
	invalidLoci := make([]string, 0)
	tooFewVariants := make([]string, 0)
	highDensity := make([]string, 0)
	bigDiff := make([]string, 0)
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
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalogFolder, id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error:", err)
			invalidLoci = append(invalidLoci, id)
			continue idLoop
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
		maxw := p.FindMaxAbsoluteWeight() * math.Pow(10, float64(p.WeightPrecision))
		if float64(p.NumVariants)/Log3(maxw) > 2.5 {
			highDensity = append(highDensity, id)
			fmt.Printf("N=%d, W=%f, N/log3(W)=%.2f\n", p.NumVariants, maxw, float64(p.NumVariants)/Log3(maxw))
			continue idLoop
		}
		if p.WeightPrecision-p.MinPrecision > 5 && p.MinPrecision > 2 && p.MinPrecision < 10 {
			fmt.Printf("------- Big difference in precision: %d, %d\n", p.WeightPrecision, p.MinPrecision)
			bigDiff = append(bigDiff, id)
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
	fmt.Println("PGS with less than", limit, "variants:", len(ids))
	fmt.Println("PGS with X or Y chromosomes:", len(hasXYchromosomes), hasXYchromosomes)
	fmt.Println("PGS with less than 2 variants:", len(tooFewVariants), tooFewVariants)
	fmt.Println("PGS with too high density:", len(highDensity), highDensity)
	fmt.Println("PGS with big precision difference:", len(bigDiff), bigDiff)
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
		src := filepath.Join("catalog", id+"_hmPOS_GRCh37.txt")
		dst := filepath.Join("filtered", id+"_hmPOS_GRCh37.txt")
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
	lociTotal := 0
	uniqueLoci := make(map[string]struct{})
	discarded := make([]string, 0)
	validated := make(map[string]int)
	for id, num := range idsToNumVariants {
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join("catalog", id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		err = p.LoadStats(data.GG)
		if err != nil {
			fmt.Printf("Discarding %s: %v\n", id, err)
			discarded = append(discarded, p.PgsID)
			continue
		}
		err = copyFile(filepath.Join("catalog", id+"_hmPOS_GRCh37.txt"),
			filepath.Join("inputs", id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error copying file:", err)
		}
		validated[p.PgsID] = num
		lociTotal += num
		for _, locus := range p.Loci {
			uniqueLoci[locus] = struct{}{}
		}
	}
	fmt.Printf("Discarded %d: %v\n", len(discarded), discarded)
	fmt.Printf("Validated %d\n", len(validated))
	fmt.Printf("Total loci: %d, unique %d\n", lociTotal, len(uniqueLoci))

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
		err = p.LoadCatalogFile(filepath.Join("inputs", id+"_hmPOS_GRCh37.txt"))
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

func downloadScoreFiles(limit int) {
	ids, err := fewerVariantsPGS(0, limit)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	localFolder := "catalog"
	ftpServer := "https://ftp.ebi.ac.uk"
	ftpFolder := "/pub/databases/spot/pgs/scores/"
	var out *os.File
	var url string
	for _, id := range ids {
		fmt.Printf("==== %s ====\n", id)
		// Create the output file
		filename := filepath.Join(localFolder, id+"_hmPOS_GRCh37.txt")
		out, err = os.Create(filename)
		if err != nil {
			log.Println("Error creating output file:", err)
			return
		}
		defer out.Close()

		// Send GET request to download the file
		url = fmt.Sprintf("%s%s%s/ScoringFiles/Harmonized/%s_hmPOS_GRCh37.txt.gz", ftpServer, ftpFolder, id, id)
		resp, err := http.Get(url)
		if err != nil {
			log.Println("Error sending GET request:", err)
		}
		defer resp.Body.Close()

		// Check if response status code is OK
		if resp.StatusCode != http.StatusOK {
			log.Printf("%s: Unexpected HTTP status code %d\n", id, resp.StatusCode)
			continue
		}

		// Decompress the HTTP response body
		decompressedData, err := decompressHTTPResponse(resp.Body)
		if err != nil {
			fmt.Printf("%s: Error decompressing file %v\n", id, err)
			continue
		}
		// Write the decompressed data to a file
		_, err = out.Write(decompressedData)
		if err != nil {
			fmt.Printf("%s: Error writing decompressed data to file %v\n", id, err)
			continue
		}
	}
}

func decompressHTTPResponse(body io.Reader) ([]byte, error) {
	// Create a gzip reader for the HTTP response body
	gzReader, err := gzip.NewReader(body)
	if err != nil {
		return nil, err
	}
	defer gzReader.Close()

	// Read and decompress the data from the gzip reader
	var buf bytes.Buffer
	_, err = io.Copy(&buf, gzReader)
	if err != nil {
		return nil, err
	}

	return buf.Bytes(), nil
}

func fewerVariantsPGS(lowerLimit, upperLimit int) ([]string, error) {
	file, err := os.Open("catalog/pgs_all_metadata_scores.csv")
	if err != nil {
		log.Println("Error opening catalog metadata file:", err)
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	header, err := reader.Read()
	if err != nil {
		fmt.Println("Error reading header:", err)
		return nil, err
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
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		if numVariants, err := strconv.Atoi(record[numVariantColumn]); err == nil && lowerLimit < numVariants &&
			numVariants < upperLimit {
			ids = append(ids, record[pgsIdColumn])
		}
	}
	return ids, nil
}

func ancestryTrainingDistribution() {
	pgsFile, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids pgsFile:", err)
		return
	}
	defer pgsFile.Close()
	decoder := json.NewDecoder(pgsFile)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}

	metafile, err := os.Open("catalog/pgs_all_metadata_scores.csv")
	if err != nil {
		log.Println("Error opening catalog metadata metafile:", err)
		return
	}
	defer metafile.Close()

	reader := csv.NewReader(metafile)
	header, err := reader.Read()
	if err != nil {
		fmt.Println("Error reading header:", err)
		return
	}
	ancestryColumn, pgsIdColumn := -1, -1
	for i, field := range header {
		if field == "Ancestry Distribution (%) - Source of Variant Associations (GWAS)" {
			ancestryColumn = i
		}
		if field == "Polygenic Score (PGS) ID" {
			pgsIdColumn = i
		}
	}
	ancestryData := make(map[string]float64)
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return
		}
		pgsId := record[pgsIdColumn]
		if _, ok := idsToNumVariants[pgsId]; !ok {
			continue
		}
		portions := strings.Split(record[ancestryColumn], "|")
		for _, portion := range portions {
			ancestryAndPercentage := strings.Split(portion, ":")
			if len(ancestryAndPercentage) != 2 {
				continue
			}
			percentage, err := strconv.ParseFloat(ancestryAndPercentage[1], 64)
			if err != nil {
				log.Println("Error converting percentage to int:", err)
				continue
			}
			ancestryData[ancestryAndPercentage[0]] += percentage
		}
	}
	var totalPercentage float64 = 0
	for _, percentage := range ancestryData {
		totalPercentage += percentage
	}
	for ancestry, percentage := range ancestryData {
		fmt.Printf("%s: %.1f\n", ancestry, percentage*100/totalPercentage)
	}
}

func weightPrecisionDistribution() {
	pgsFile, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids pgsFile:", err)
		return
	}
	defer pgsFile.Close()
	decoder := json.NewDecoder(pgsFile)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}

	precisions := make(map[uint32]int)
	for id := range idsToNumVariants {
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join("inputs", id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error:", err)
			return
		}
		precisions[p.WeightPrecision]++
		err = p.LoadStats(data.GG)
		if err != nil {
			log.Printf("%s: Error loading stats: %v\n", p.PgsID, err)
		}
	}
	// Sort and print the distribution of weight precisions
	var keys []uint32
	total := 0
	for k := range precisions {
		keys = append(keys, k)
		total += precisions[k]
	}
	sort.Slice(keys, func(i, j int) bool {
		return precisions[keys[i]] > precisions[keys[j]]
	})
	for _, k := range keys {
		fmt.Printf("%d: %.0f\n", k, float64(precisions[k])*100/float64(total))
	}
}

func Log3(x float64) float64 {
	return math.Log2(x) / math.Log2(3)
}

func readSamplePopulations() {
	fileName := "info/igsr_samples.tsv"

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

	// Save the results to a JSON file
	outputFile, err := os.Create("info/superpopulations.json")
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
	expr := flag.String("expr", "", "Experiment type")
	flag.Parse()
	switch *expr {
	case "download":
		downloadScoreFiles(500)
	case "limit":
		getPGSWithFewerVariants(50)
	case "copy":
		copyFilteredPGSFiles()
	case "validate":
		validateBases()
	case "index":
		makeLociIndex()
	case "ancestry":
		ancestryTrainingDistribution()
	case "precision":
		weightPrecisionDistribution()
	case "populations":
		readSamplePopulations()
	}
}
