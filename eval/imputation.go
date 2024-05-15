package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/tools"
	"log"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func prepareVCF() {
	dummyInfo := "\t.\t.\t.\t100\tPASS\t.\tGT"
	var err error
	var file *os.File
	var decoder *json.Decoder
	folder := "results/sequential"
	dir, err := os.ReadDir(folder)
	if err != nil {
		log.Fatalf("Cannot read directory %s: %v", folder, err)
	}

	type Guess struct {
		Individual string
		SNPs       map[string]uint8
	}
	input := make(map[string]map[string]uint8)
	for _, object := range dir {
		if !object.IsDir() && strings.HasPrefix(object.Name(), "guesses") && strings.HasSuffix(object.Name(), ".json") {
			guesses := make([]*Guess, 0)
			file, err = os.Open(filepath.Join(folder, object.Name()))
			if err != nil {
				log.Fatalf("Cannot open file %s: %v", object.Name(), err)
			}
			defer file.Close()
			decoder = json.NewDecoder(file)
			if err = decoder.Decode(&guesses); err != nil {
				log.Fatalf("Cannot decode json file %s: %v", object.Name(), err)
			}
			for _, g := range guesses {
				input[g.Individual] = g.SNPs
			}
		}
	}

	vcfFile, err := os.Open("chr22.vcf")
	if err != nil {
		log.Fatalf("Cannot open VCF file: %v", err)
	}
	defer vcfFile.Close()
	imputeFile, err := os.OpenFile("chr22-impute.vcf", os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Cannot create impute file: %v", err)
	}
	defer imputeFile.Close()

	var locus, chr, pos string
	var posInt int
	allPositions := make([]int, 0)
	example := "HG00124"
	for locus := range input[example] {
		chr, pos = tools.SplitLocus(locus)
		if chr != "22" {
			continue
		}
		posInt, err = strconv.Atoi(pos)
		if err != nil {
			log.Fatalf("Cannot convert position %s to integer: %v", pos, err)
		}
		allPositions = append(allPositions, posInt)
	}
	sort.Ints(allPositions)
	fmt.Printf("Number of guessed positions in chr 22: %d\n", len(allPositions))
	fmt.Printf("All positions: %v\n", allPositions)

	var fields []string
	var newLine string
	var snp uint8
	var ok bool
	var ptr, totalSNPs = 0, 0
	//accuracy := make(map[string]float32)
	individuals := make([]string, 0)
	scanner := bufio.NewScanner(vcfFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			_, err = imputeFile.WriteString(line + "\n")
			if err != nil {
				log.Fatalf("Cannot write to impute file: %v", err)
			}
			if strings.HasPrefix(line, "#CHROM") {
				fields := strings.Split(line, "\t")
				for _, field := range fields[9:] {
					individuals = append(individuals, field)
				}
			}
			continue
		}
		fields = strings.Split(line, "\t")
		chr = fields[0]
		pos = fields[1]
		posInt, err = strconv.Atoi(pos)
		if err != nil {
			log.Fatalf("Cannot convert position %s to integer: %v", pos, err)
		}
		for ptr < len(allPositions) && posInt > allPositions[ptr] {
			// Locus is not in the original VCF file, but we have it guessed
			fmt.Printf("Locus %s:%d is not in the original VCF file, current %s\n", chr, allPositions[ptr], pos)
			newLine = "22\t" + strconv.Itoa(allPositions[ptr]) + dummyInfo
			for _, idv := range individuals {
				if _, ok = input[idv]; !ok {
					newLine += "\t./."
					continue
				}
				if snp, ok = input[idv]["22:"+strconv.Itoa(allPositions[ptr])]; ok {
					switch snp {
					case 0:
						newLine += "\t0/0"
					case 1:
						newLine += "\t1/0"
					case 2:
						newLine += "\t1/1"
					default:
						log.Fatalf("Invalid SNP value: %d", snp)
					}
				} else {
					log.Printf("No SNP for individual %s at locus %s\n", idv, "22:"+strconv.Itoa(allPositions[ptr]))
					newLine += "\t./."
				}
			}
			_, err = imputeFile.WriteString(newLine + "\n")
			if err != nil {
				log.Fatalf("Cannot write to impute file %s: %v", newLine, err)
			}
			totalSNPs++
			ptr++
		}
		if ptr < len(allPositions) && posInt == allPositions[ptr] {
			ptr++
		}
		locus = chr + ":" + pos
		newLine = strings.Join(fields[:9], "\t")
		for _, idv := range individuals {
			if _, ok = input[idv]; !ok {
				newLine += "\t./."
				continue
			}
			if snp, ok = input[idv][locus]; ok {
				switch snp {
				case 0:
					newLine += "\t0/0"
				case 1:
					newLine += "\t1/0"
				case 2:
					newLine += "\t1/1"
				default:
					log.Fatalf("Invalid SNP value: %d", snp)
				}
			} else {
				newLine += "\t./."
			}
		}
		_, err = imputeFile.WriteString(newLine + "\n")
		if err != nil {
			log.Fatalf("Cannot write to impute file %s: %v", newLine, err)
		}
		totalSNPs++
	}
	if err := scanner.Err(); err != nil {
		log.Printf("Scanning error: %v\n", err)
	}
	fmt.Printf("Total number of SNPs: %d\n", totalSNPs)
}

func imputedAccuracy() {
	originalFile, err := os.Open("chr22.vcf")
	if err != nil {
		log.Fatalf("Cannot open VCF file: %v", err)
	}
	defer originalFile.Close()
	var fields []string
	var pos, normalizedSnp string
	var snp uint8
	var individuals []string
	trueSNPs := make(map[string]map[string]uint8)
	scanner := bufio.NewScanner(originalFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			if strings.HasPrefix(line, "#CHROM") {
				fields := strings.Split(line, "\t")
				for _, field := range fields[9:] {
					individuals = append(individuals, field)
				}
			}
			continue
		}
		fields = strings.Split(line, "\t")
		pos = fields[1]
		for i, idv := range individuals {
			if _, ok := trueSNPs[idv]; !ok {
				trueSNPs[idv] = make(map[string]uint8)
			}
			normalizedSnp, err = tools.NormalizeSnp(fields[9+i])
			if err != nil {
				log.Fatalf("Cannot normalize SNP %s: %v", fields[9+i], err)
			}
			snp, err = tools.SnpToSum(normalizedSnp)
			if err != nil {
				log.Fatalf("Cannot convert SNP %s to sum: %v", normalizedSnp, err)
			}
			trueSNPs[idv][pos] = snp
		}
	}

	imputedFile, err := os.Open("chr22-result.vcf")
	if err != nil {
		log.Fatalf("Cannot open imputed VCF file: %v", err)
	}
	defer imputedFile.Close()
	var alleles []string
	var parsed float64
	imputedSNPs := make(map[string]map[string]uint8)
	scanner = bufio.NewScanner(imputedFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			continue
		}
		fields = strings.Split(line, "\t")
		pos = fields[1]
		for i, idv := range individuals {
			if _, ok := imputedSNPs[idv]; !ok {
				imputedSNPs[idv] = make(map[string]uint8)
			}
			alleles = strings.Split(fields[9+i], ",")
			snp = 0
			for _, allele := range alleles {
				parsed, err = strconv.ParseFloat(allele, 64)
				if err != nil {
					log.Fatalf("Cannot parse allele %s: %v", allele, err)
				}
				if parsed >= 0.5 {
					snp += 1
				}
			}
			imputedSNPs[idv][pos] = snp
		}
	}

	//originalLoci := []string{"17649774", "19872645", "24159307", "24171305", "26164079", "28552698", "28934313",
	//"29121087", "29300306", "30254994", "30416527", "30531091", "30592069", "35700467", "37319009", "37469591",
	//"37534034", "38477930", "38544298", "38545619", "38545942", "38563471", "39332623", "39542292", "39546145",
	//"40932041", "41009707", "41023304", "41621714", "43500212", 44324727 44324730 44324855 44340904 44342116 45529171 46615880 50356693 50971266 51156933}
	for idv := range trueSNPs {
		numMatches, totalSNPs := 0, 0
		for pos := range trueSNPs[idv] {
			if trueSNPs[idv][pos] == imputedSNPs[idv][pos] {
				numMatches++
			}
			totalSNPs++
		}
		fmt.Printf("Individual %s: %d matches out of %d, accuracy %.3f\n", idv, numMatches, totalSNPs, float64(numMatches)/float64(totalSNPs))
	}
}
