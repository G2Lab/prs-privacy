package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/tools"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

var (
	preImputeUnzippedFilename  = "imputation/chr%d-%s-preimpute.vcf"
	postImputeUnzippedFilename = "imputation/chr%d-%s-postimpute.vcf"
	preImputeZippedFilename    = preImputeUnzippedFilename + ".gz"
	postImputeZippedFilename   = postImputeUnzippedFilename + ".gz"
	groundTruthFilename        = "imputation/chr%d-truth.vcf"
)

func imputeWorkflow() {
	for _, chrId := range []int{6} {
		//for _, ancestry := range pgs.POPULATIONS {
		for _, ancestry := range []string{"EUR"} {
			fmt.Printf("===== %s =====\n", ancestry)
			fillPreImputeVCF(chrId, ancestry)
			compressVCF(chrId, ancestry)
			indexCompressedVCF(chrId, ancestry, false)
			impute(ancestry, chrId)
			indexCompressedVCF(chrId, ancestry, true)
		}
	}
}

func impute(ancestry string, chrId int) {
	prg := "minimac4"
	args := []string{
		fmt.Sprintf("data/references/%d.msav", chrId),
		"-t", "4",
		"--sample-ids-file", fmt.Sprintf("data/1000g_%s_no_relatives.txt", ancestry),
		fmt.Sprintf(preImputeZippedFilename, chrId, ancestry),
		"-o", fmt.Sprintf(postImputeZippedFilename, chrId, ancestry),
	}
	cmd := exec.Command(prg, args...)
	output, err := cmd.Output()
	fmt.Println(string(output))
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func compressVCF(chrId int, ancestry string) {
	prg := "bgzip"
	args := []string{"-c", fmt.Sprintf(preImputeUnzippedFilename, chrId, ancestry), ">",
		fmt.Sprintf(preImputeZippedFilename, chrId, ancestry)}
	cmd := exec.Command(prg, args...)
	output, err := cmd.Output()
	fmt.Println(string(output))
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func indexCompressedVCF(chrId int, ancestry string, postImpute bool) {
	var filename = preImputeZippedFilename
	if postImpute {
		filename = fmt.Sprintf(postImputeZippedFilename, chrId, ancestry)
	}
	prg := "tabix"
	args := []string{"-p", "vcf", "-f", fmt.Sprintf(filename, chrId, ancestry)}
	cmd := exec.Command(prg, args...)
	output, err := cmd.Output()
	fmt.Println(string(output))
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func fillPreImputeVCF(chrId int, ancestry string) {
	var err error
	var file *os.File
	var decoder *json.Decoder
	folder := "results/sequential"
	dir, err := os.ReadDir(folder)
	if err != nil {
		log.Fatalf("Cannot read directory %s: %v", folder, err)
	}

	ancestries := tools.LoadAncestry()
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
			decoder = json.NewDecoder(file)
			if err = decoder.Decode(&guesses); err != nil {
				log.Fatalf("Cannot decode json file %s: %v", object.Name(), err)
			}
			for _, g := range guesses {
				if ancestry != getIndividualAncestry(g.Individual, ancestries) {
					continue
				}
				input[g.Individual] = g.SNPs
			}
			file.Close()
		}
	}

	truthFile, err := os.Open(fmt.Sprintf(groundTruthFilename, chrId))
	if err != nil {
		log.Fatalf("Cannot open VCF file: %v", err)
	}
	defer truthFile.Close()
	preImputeFile, err := os.OpenFile(fmt.Sprintf(preImputeUnzippedFilename, chrId, ancestry), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Cannot create impute file: %v", err)
	}
	defer preImputeFile.Close()

	var locus, chr, pos string
	var posInt int
	var firstIdv string
	allPositions := make([]int, 0)
	for idv := range input {
		firstIdv = idv
		break
	}
	for locus := range input[firstIdv] {
		chr, pos = tools.SplitLocus(locus)
		if chr != strconv.Itoa(chrId) {
			continue
		}
		posInt, err = strconv.Atoi(pos)
		if err != nil {
			log.Fatalf("Cannot convert position %s to integer: %v", pos, err)
		}
		allPositions = append(allPositions, posInt)
	}
	sort.Ints(allPositions)
	fmt.Printf("Number of guessed positions in chr %d: %d\n", chrId, len(allPositions))
	fmt.Printf("All positions: %v\n", allPositions)

	var fields []string
	var newLine string
	var snp uint8
	var ok bool
	var samplesStartPos int
	var ptr, totalSNPs = 0, 0
	individuals := make([]string, 0)
	scanner := bufio.NewScanner(truthFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			if strings.HasPrefix(line, "#CHROM") {
				fields := strings.Split(line, "\t")
				for i, field := range fields {
					if samplePredicate(field) && getIndividualAncestry(field, ancestries) == ancestry {
						if samplesStartPos == 0 {
							samplesStartPos = i
							line = strings.Join(fields[:samplesStartPos], "\t")
						}
						individuals = append(individuals, field)
						line += "\t" + field
					}
				}
			}
			_, err = preImputeFile.WriteString(line + "\n")
			if err != nil {
				log.Fatalf("Cannot write to impute file: %v", err)
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
			newLine = fmt.Sprintf("%d\t%d\t.\t.\t.\t100\tPASS\t.\tGT", chrId, allPositions[ptr])
			for _, idv := range individuals {
				if _, ok = input[idv]; !ok {
					newLine += "\t./."
					continue
				}
				if snp, ok = input[idv][fmt.Sprintf("%d:%d", chrId, allPositions[ptr])]; ok {
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
					log.Printf("No SNP for individual %s at locus %s\n", idv, fmt.Sprintf("%d:%d", chrId, allPositions[ptr]))
					newLine += "\t./."
				}
			}
			_, err = preImputeFile.WriteString(newLine + "\n")
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
		newLine = strings.Join(fields[:samplesStartPos], "\t")
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
		_, err = preImputeFile.WriteString(newLine + "\n")
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

func imputedAccuracy(chrId int, ancestry string) {
	originalFile, err := os.Open(fmt.Sprintf(groundTruthFilename, chrId))
	if err != nil {
		log.Fatalf("Cannot open VCF file: %v", err)
	}
	defer originalFile.Close()
	var fields []string
	var chrPos, normalizedSnp string
	var snp uint8
	idvOriginalPos := make(map[string]int)
	ancestries := tools.LoadAncestry()
	trueSNPs := make(map[string]map[string]uint8)
	scanner := bufio.NewScanner(originalFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			if strings.HasPrefix(line, "#CHROM") {
				fields := strings.Split(line, "\t")
				for i, field := range fields {
					if samplePredicate(field) && getIndividualAncestry(field, ancestries) == ancestry {
						idvOriginalPos[field] = i
					}
				}
			}
			continue
		}
		fields = strings.Split(line, "\t")
		chrPos = fields[1]
		for idv, pos := range idvOriginalPos {
			if _, ok := trueSNPs[idv]; !ok {
				trueSNPs[idv] = make(map[string]uint8)
			}
			normalizedSnp, err = tools.NormalizeSnp(fields[pos])
			if err != nil {
				log.Fatalf("Cannot normalize SNP %s: %v", fields[pos], err)
			}
			snp, err = tools.SnpToSum(normalizedSnp)
			if err != nil {
				log.Fatalf("Cannot convert SNP %s to sum: %v", normalizedSnp, err)
			}
			trueSNPs[idv][chrPos] = snp
		}
	}
	fmt.Printf("Original positions: %v\n", idvOriginalPos)

	imputedFile, err := os.Open(fmt.Sprintf(postImputeUnzippedFilename, chrId, ancestry))
	if err != nil {
		log.Fatalf("Cannot open imputed VCF file: %v", err)
	}
	defer imputedFile.Close()
	var alleles []string
	var parsed float64
	idvImputedPos := make(map[string]int)
	imputedSNPs := make(map[string]map[string]uint8)
	scanner = bufio.NewScanner(imputedFile)
	for scanner.Scan() {
		line := scanner.Text()
		// If it is the header
		if strings.HasPrefix(line, "#") {
			if strings.HasPrefix(line, "#CHROM") {
				fields := strings.Split(line, "\t")
				for i, field := range fields {
					if _, ok := idvOriginalPos[field]; ok {
						idvImputedPos[field] = i
					}
				}
				fmt.Printf("Imputed positions: %v\n", idvImputedPos)
			}
			continue
		}
		fields = strings.Split(line, "\t")
		chrPos = fields[1]
		for idv, pos := range idvImputedPos {
			if _, ok := imputedSNPs[idv]; !ok {
				imputedSNPs[idv] = make(map[string]uint8)
			}
			if len(fields) < pos {
				log.Printf("Not enough fields in line:\n%s\n", line)
				continue
			}
			alleles = strings.Split(fields[pos], ",")
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
			imputedSNPs[idv][chrPos] = snp
		}
	}

	var exists bool
	for idv := range trueSNPs {
		numMatches, totalSNPs := 0, 0
		for pos := range trueSNPs[idv] {
			if _, exists = imputedSNPs[idv][pos]; exists && trueSNPs[idv][pos] == imputedSNPs[idv][pos] {
				numMatches++
			}
			totalSNPs++
		}
		fmt.Printf("Individual %s: %d matches out of %d, accuracy %.3f\n", idv, numMatches, totalSNPs, float64(numMatches)/float64(totalSNPs))
	}
}

func samplePredicate(input string) bool {
	return strings.HasPrefix(input, "HG") || strings.HasPrefix(input, "NA")
}

func getIndividualAncestry(idv string, ancestry map[string]string) string {
	if _, ok := ancestry[idv]; !ok {
		log.Fatalf("Individual %s not found in ancestry file\n", idv)
	}
	idvAnc := ancestry[idv]
	if strings.Contains(idvAnc, ",") {
		idvAnc = strings.Split(idvAnc, ",")[0]
	}
	return idvAnc
}

func loadIndividualAncestry(idv string) string {
	return getIndividualAncestry(idv, tools.LoadAncestry())
}
