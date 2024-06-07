package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"
)

var (
	preImputeUnzippedFilename  = "imputation/chr%d-%s-preimpute.vcf"
	postImputeUnzippedFilename = "imputation/chr%d-%s-postimpute.vcf"
	preImputeZippedFilename    = preImputeUnzippedFilename + ".gz"
	postImputeZippedFilename   = postImputeUnzippedFilename + ".gz"
	groundTruthFilename        = "imputation/chr%d-truth.vcf"
)

const (
	imputationChunkSize    = 100000
	imputationChunkOverlap = 10000
)

func imputeWorkflow() {
	for chrId := 1; chrId <= 22; chrId++ {
		fmt.Printf("---- %d ----\n", chrId)
		for _, ancestry := range pgs.POPULATIONS {
			fmt.Printf("===== %s =====\n", ancestry)
			fillPreImputeVCF(chrId, ancestry)
			fmt.Printf("Pre-impute VCF filled\n")
			compressVCF(chrId, ancestry)
			fmt.Printf("VCF compressed\n")
			indexCompressedVCF(chrId, ancestry, false)
			fmt.Printf("VCF indexed\n")
			impute(chrId, ancestry)
			fmt.Printf("Imputation complete\n")
			indexCompressedVCF(chrId, ancestry, true)
			fmt.Printf("Imputed VCF indexed\n")
		}
	}
	fmt.Println("Completed")
}

func impute(chrId int, ancestry string) {
	prg := "../minimac4/bin/minimac4"
	args := []string{
		fmt.Sprintf("data/references/%d.msav", chrId),
		"-t", "8",
		"-a",
		"--sample-ids-file", fmt.Sprintf("data/1000g_%s_no_relatives.txt", ancestry),
		"--format", "GT",
		"-c", strconv.Itoa(imputationChunkSize),
		"--overlap", strconv.Itoa(imputationChunkOverlap),
		"--min-ratio", "1e-6",
		"--min-ratio-behavior", "skip",
		fmt.Sprintf(preImputeZippedFilename, chrId, ancestry),
		"-o", fmt.Sprintf(postImputeZippedFilename, chrId, ancestry),
		"-O", "vcf.gz",
	}
	outputFile, err := os.OpenFile(fmt.Sprintf(postImputeZippedFilename, chrId, ancestry), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error creating output file: %v", err)
	}
	defer outputFile.Close()
	err = exec.Command(prg, args...).Run()
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func compressVCF(chrId int, ancestry string) {
	prg := "bgzip"
	cmd := exec.Command(prg, "-c", fmt.Sprintf(preImputeUnzippedFilename, chrId, ancestry))
	outputFile, err := os.OpenFile(fmt.Sprintf(preImputeZippedFilename, chrId, ancestry), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error creating output file: %v", err)
	}
	defer outputFile.Close()
	cmd.Stdout = outputFile
	err = cmd.Run()
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func indexCompressedVCF(chrId int, ancestry string, postImpute bool) {
	var filename = preImputeZippedFilename
	if postImpute {
		filename = postImputeZippedFilename
	}
	prg := "tabix"
	args := []string{"-f", "-p", "vcf", fmt.Sprintf(filename, chrId, ancestry)}
	err := exec.Command(prg, args...).Run()
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
}

func fillPreImputeVCF(chrId int, ancestry string) {
	guessed := loadGuessedGenotypes()
	ancestries := tools.LoadAncestry()
	chr := strconv.Itoa(chrId)
	var pos, ref, alt string
	var posInt int
	individuals := make([]string, 0)
	positionMap := make(map[string]struct{})
	for idv := range guessed {
		if ancestry == getIndividualAncestry(idv, ancestries) {
			individuals = append(individuals, idv)
			for pos = range guessed[idv][chr] {
				positionMap[pos] = struct{}{}
			}
		}
	}
	var err error
	chrPositions := make([]int, 0)
	alleles := make(map[int][]string)
	for pos = range positionMap {
		posInt, err = strconv.Atoi(pos)
		if err != nil {
			log.Fatalf("Cannot convert position %s to integer: %v", pos, err)
		}
		chrPositions = append(chrPositions, posInt)
		ref, alt = retrievePositionAlleles(chr, pos)
		alleles[posInt] = []string{ref, alt}
	}
	sort.Ints(chrPositions)
	fmt.Printf("Number of guessed positions in chr %d: %d\n", chrId, len(chrPositions))
	fmt.Printf("All positions: %v\n", chrPositions)

	formatHeader := "##fileformat=VCFv4.1\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	samplesHeader := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
	for _, idv := range individuals {
		samplesHeader += fmt.Sprintf("\t%s", idv)
	}
	samplesHeader += "\n"
	preImputeFile, err := os.OpenFile(fmt.Sprintf(preImputeUnzippedFilename, chrId, ancestry), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Cannot create impute file: %v", err)
	}
	defer preImputeFile.Close()
	_, err = preImputeFile.WriteString(formatHeader)
	if err != nil {
		log.Fatalf("Cannot write format header to impute file: %v", err)
	}
	_, err = preImputeFile.WriteString(samplesHeader)
	if err != nil {
		log.Fatalf("Cannot write samples header to impute file: %v", err)
	}
	lineTemplate := "%d\t%d\t.\t%s\t%s\t100\tPASS\t.\tGT"

	var snp uint8
	var ok bool
	var line string
	for i, pos := range chrPositions {
		line = fmt.Sprintf(lineTemplate, chrId, chrPositions[i], alleles[pos][0], alleles[pos][1])
		for _, idv := range individuals {
			if _, ok = guessed[idv]; !ok {
				line += "\t./."
				continue
			}
			if snp, ok = guessed[idv][chr][fmt.Sprintf("%d", chrPositions[i])]; ok {
				switch snp {
				case 0:
					line += "\t0/0"
				case 1:
					line += "\t1/0"
				case 2:
					line += "\t1/1"
				default:
					log.Fatalf("Invalid SNP value: %d", snp)
				}
			} else {
				log.Printf("No SNP for individual %s at position %s\n", idv, fmt.Sprintf("%d", chrPositions[i]))
				line += "\t./."
			}
		}
		_, err = preImputeFile.WriteString(line + "\n")
		if err != nil {
			log.Fatalf("Cannot write to impute file %s: %v", line, err)
		}
	}
}

func retrievePositionAlleles(chr, pos string) (string, string) {
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS\t%REF\t%ALT\n",
		"-r", fmt.Sprintf("%s:%s-%s", chr, pos, pos),
		tools.GetChromosomeFilepath(chr, tools.GG),
	}
	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing command: %v", err)
	}
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]
	var fields []string
	for _, line := range lines {
		fields = strings.Split(line, "\t")
		if fields[0] != pos {
			continue
		}
		return fields[1], fields[2]
	}
	return ".", "."
}

func imputationAccuracy() {
	var chrId int
	var err error
	if len(os.Args) > 1 {
		chrId, err = strconv.Atoi(os.Args[1])
		if err != nil {
			log.Fatalf("Cannot convert chromosome ID to integer: %v", err)
		}
	} else {
		chrId = 22
	}
	fmt.Printf("Imputation accuracy for chromosome %d\n", chrId)
	numWorkers := 15
	type Relation struct {
		Target string
		Count  []int
		King   []int
	}
	imputedSNPs := make(map[string]map[string]uint8)
	for _, ancestry := range pgs.POPULATIONS {
		fmt.Printf("Reading imputed SNPs: %d, %s\n", chrId, ancestry)
		chrPath := fmt.Sprintf(postImputeZippedFilename, chrId, ancestry)
		imputedChunk := getAllIndividualAlleles(chrPath, numWorkers)
		for idv := range imputedChunk {
			imputedSNPs[idv] = imputedChunk[idv]
		}
	}
	individuals := make([]string, 0)
	for idv := range imputedSNPs {
		individuals = append(individuals, idv)
	}
	db := solver.All1000GenomesSamples()
	relations := make(map[string][]Relation)
	for idv := range imputedSNPs {
		relations[idv] = make([]Relation, 0)
	}

	// Read all the imputed positions, sort them and divide them into regions
	imputedPositions := getAllImputedPositions(chrId)
	imputedRegions := dividePositionsIntoRegions(imputedPositions, imputationChunkSize)
	fmt.Printf("Number of imputed SNPs: %d\n", len(imputedPositions))
	fmt.Printf("Imputation regions: %v\n", imputedRegions)

	var mu sync.Mutex
	var wg sync.WaitGroup
	var taskWg sync.WaitGroup
	// Channel to send work to goroutines
	type task struct {
		idv   string
		other string
		snps  map[string]map[string]uint8
	}
	taskChan := make(chan task)

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for t := range taskChan {
				rel := Relation{t.other, matchCount(imputedSNPs[t.idv], t.snps[t.other]),
					kingRobustPair(imputedSNPs[t.idv], t.snps[t.other])}
				mu.Lock()
				relations[t.idv] = append(relations[t.idv], rel)
				mu.Unlock()
				taskWg.Done()
			}
		}()
	}
	go func() {
		chunkSize := 200
		for i := 0; i < len(db); i += chunkSize {
			fmt.Printf("==== Chunk %d ====\n", i/chunkSize)
			chunk := db[i:min(i+chunkSize, len(db))]
			chrPath := tools.GetChromosomeFilepath(strconv.Itoa(chrId), tools.GG)
			chunkSnps := getRegionIndividualAlleles(strconv.Itoa(chrId), chrPath, chunk, imputedPositions,
				imputedRegions, numWorkers)
			for idv := range imputedSNPs {
				for other := range chunkSnps {
					taskWg.Add(1)
					taskChan <- task{idv, other, chunkSnps}
				}
			}
			taskWg.Wait()
		}
		close(taskChan)
	}()
	wg.Wait()

	references := make(map[string]map[string]uint8)
	for _, ppl := range pgs.POPULATIONS {
		references[ppl] = allMajorAlleleSample(chrId, ppl)
	}

	ancestries := tools.LoadAncestry()
	relatives := solver.AllRelativeSamples()
	chrStr := strconv.Itoa(chrId)
	relativeSnps := getRegionIndividualAlleles(chrStr, tools.GetChromosomeFilepath(chrStr, tools.RL), relatives,
		imputedPositions, imputedRegions, numWorkers)
	for idv := range imputedSNPs {
		for other := range relativeSnps {
			relations[idv] = append(relations[idv], Relation{other,
				matchCount(imputedSNPs[idv], relativeSnps[other]), kingRobustPair(imputedSNPs[idv], relativeSnps[other])})
		}
		relations[idv] = append(relations[idv], Relation{"reference",
			matchCount(imputedSNPs[idv], references[getIndividualAncestry(idv, ancestries)]),
			kingRobustPair(imputedSNPs[idv], references[getIndividualAncestry(idv, ancestries)])})
	}
	resultFolder := "results/impute"
	filepath := path.Join(resultFolder, fmt.Sprintf("imputed%d.json", chrId))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(relations); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Completed")
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func calculateGenotypeFrequencies() {
	numWorkers := 5
	folder := "data/frequencies"
	var filepath string
	var outFile *os.File
	var err error
	var ok bool
	var positions []int
	var regions [][]string
	var alleles map[string][]uint8
	var frequencies map[string][]float32
	var chrStr string
	for chr := 1; chr <= 22; chr++ {
		fmt.Printf("---- %d ----\n", chr)
		frequencies = make(map[string][]float32)
		chrStr = strconv.Itoa(chr)
		positions = getAllImputedPositions(chr)
		regions = dividePositionsIntoRegions(positions, imputationChunkSize)
		//fmt.Printf("Positions: %v\n", positions)
		//fmt.Printf("Regions: %v\n", regions)
		alleles = getRegionAllelesPerLocus(chrStr, tools.GetChromosomeFilepath(chrStr, tools.GG), regions, numWorkers)
		for locus := range alleles {
			frequencies[locus] = make([]float32, 3)
			for _, allele := range alleles[locus] {
				frequencies[locus][allele]++
			}
			for i := range frequencies[locus] {
				frequencies[locus][i] /= float32(len(alleles[locus]))
			}
		}
		for _, pos := range positions {
			if _, ok = frequencies[strconv.Itoa(pos)]; !ok {
				fmt.Printf("Position %d was missing\n", pos)
				frequencies[strconv.Itoa(pos)] = []float32{(1 - pgs.MissingEAF) * (1 - pgs.MissingEAF),
					pgs.MissingEAF * (1 - pgs.MissingEAF), pgs.MissingEAF * pgs.MissingEAF}
			}
		}
		filepath = path.Join(folder, fmt.Sprintf("%d.json", chr))
		outFile, err = os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatalf("Error opening output file: %v", err)
		}
		encoder := json.NewEncoder(outFile)
		if err = encoder.Encode(frequencies); err != nil {
			log.Fatal("Cannot encode json", err)
		}
		outFile.Close()
		fmt.Printf("Chr %d completed\n", chr)
	}
}

func getAllImputedPositions(chr int) []int {
	imputedMap := make(map[string]struct{})
	var imputedPositionsPerAncestry []string
	for _, ancestry := range pgs.POPULATIONS {
		imputedPositionsPerAncestry = getAllPositions(fmt.Sprintf(postImputeZippedFilename, chr, ancestry))
		for _, pos := range imputedPositionsPerAncestry {
			imputedMap[pos] = struct{}{}
		}
	}
	imputedPositions := make([]int, 0, len(imputedMap))
	var posInt int
	var err error
	for pos := range imputedMap {
		posInt, err = strconv.Atoi(pos)
		if err != nil {
			log.Fatalf("Cannot convert position %s to integer: %v", pos, err)
		}
		imputedPositions = append(imputedPositions, posInt)
	}
	sort.Ints(imputedPositions)
	return imputedPositions
}

func dividePositionsIntoRegions(positions []int, regionLen int) [][]string {
	regions := make([][]string, 0)
	regionStart := 0
	for i := 1; i < len(positions); i++ {
		if positions[i]-positions[i-1] > regionLen {
			regions = append(regions, []string{strconv.Itoa(positions[regionStart]), strconv.Itoa(positions[i-1])})
			regionStart = i
		}
	}
	regions = append(regions, []string{strconv.Itoa(positions[regionStart]), strconv.Itoa(positions[len(positions)-1])})
	return regions
}

func imputationAccuracyAll(chrId int) {
	originalFile, err := os.Open(fmt.Sprintf(groundTruthFilename, chrId))
	if err != nil {
		log.Fatalf("Cannot open VCF file: %v", err)
	}
	defer originalFile.Close()
	var fields []string
	var snp uint8
	var chrPos string
	var alleles []string
	var parsed float64
	imputedSNPs := make(map[string]map[string]uint8)
	for _, ancestry := range pgs.POPULATIONS {
		imputedFile, err := os.Open(fmt.Sprintf(postImputeUnzippedFilename, chrId, ancestry))
		if err != nil {
			log.Fatalf("Cannot open imputed VCF file: %v", err)
		}
		idvImputedPos := make(map[string]int)
		scanner := bufio.NewScanner(imputedFile)
		for scanner.Scan() {
			line := scanner.Text()
			// If it is the header
			if strings.HasPrefix(line, "#") {
				if strings.HasPrefix(line, "#CHROM") {
					fields := strings.Split(line, "\t")
					for i, field := range fields {
						//if _, ok := idvOriginalPos[field]; ok {
						if samplePredicate(field) {
							idvImputedPos[field] = i
						}
					}
				}
				continue
			}
			fields = strings.Split(line, "\t")
			chrPos = fields[1]
			for idv, pos := range idvImputedPos {
				if len(fields) < pos {
					log.Printf("Not enough fields in line:\n%s\n", line)
					continue
				}
				if _, ok := imputedSNPs[idv]; !ok {
					imputedSNPs[idv] = make(map[string]uint8)
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
				if snp > 2 {
					snp = 2
				}
				imputedSNPs[idv][chrPos] = snp
			}
		}
		imputedFile.Close()
	}
	linkingAccuracy(chrId, imputedSNPs)
}

func linkingAccuracy(chrId int, imputedSNPs map[string]map[string]uint8) {
	fmt.Printf("Linking\n")
	type relation struct {
		target string
		score  float32
	}
	db := solver.All1000GenomesSamples()
	related := solver.ReadRelatedIndividuals()
	//ancestry := tools.LoadAncestry()
	//chrAF := getChromosomeAFs(strconv.Itoa(chrId))
	//fmt.Printf("AF loaded\n")

	relations := make(map[string][]relation)
	for idv := range imputedSNPs {
		relations[idv] = make([]relation, 0)
	}
	for _, other := range db {
		fmt.Printf("==== Main: %s ====\n", other)
		otherSnps := getIndividualAlleles(strconv.Itoa(chrId), other, tools.GG)
		for idv := range imputedSNPs {
			relations[idv] = append(relations[idv], relation{other, countSimilarity(imputedSNPs[idv], otherSnps)})
		}
	}

	relatives := solver.AllRelativeSamples()
	for _, relative := range relatives {
		fmt.Printf("==== Extra: %s ====\n", relative)
		otherSnps := getIndividualAlleles(strconv.Itoa(chrId), relative, tools.RL)
		for idv := range imputedSNPs {
			//relations[idv] = append(relations[idv], Relation{relative, similarity(imputedSNPs[idv], otherSnps, chrAF[ancestry[idv]])})
			relations[idv] = append(relations[idv], relation{relative, countSimilarity(imputedSNPs[idv], otherSnps)})
		}
	}

	type Result struct {
		SelfPos      int
		RelativePos  int
		SelfAccuracy float32
		RelAccuracy  float32
	}
	var selfFound, relativeFound bool
	results := make(map[string]*Result)
	for idv := range imputedSNPs {
		fmt.Printf("==== %s ====\n", idv)
		sort.Slice(relations[idv], func(i, j int) bool {
			return relations[idv][i].score > relations[idv][j].score
		})
		relatives = related[idv]
		selfFound, relativeFound = false, false
		results[idv] = new(Result)
		for i, rel := range relations[idv] {
			if rel.target == idv {
				fmt.Printf("Self: %d pos, %.0f score / %.0f top\n", i, rel.score, relations[idv][0].score)
				results[idv].SelfPos = i
				results[idv].SelfAccuracy = rel.score
				selfFound = true
			}
			for _, rltv := range relatives {
				if rel.target == rltv {
					fmt.Printf("Relative %s: %d pos, %.0f score\n", rltv, i, rel.score)
					results[idv].RelativePos = i
					results[idv].RelAccuracy = rel.score
					relativeFound = true
					break
				}
			}
			if selfFound && relativeFound {
				break
			}
		}
	}
	resultFolder := "results/impute"
	filepath := path.Join(resultFolder, "imputed.json")
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Completed")
}

func linkingWithGuessed() {
	type Relation struct {
		Target      string
		Count       []int
		King        []int
		Information []float32
	}
	guessed := loadGuessedGenotypes()
	relatives := solver.AllRelativeSamples()
	mainSamples := solver.All1000GenomesSamples()
	ancestries := tools.LoadAncestry()
	resultFolder := "results/impute"
	var encoder *json.Encoder
	var resFile *os.File
	var err error
	var chrStr, idv, pos, filepath string
	var positionMap map[string]struct{}
	var positions []string
	var relations map[string][]Relation
	var gtpFreqs map[string][]float32
	var references map[string]map[string]uint8
	var otherSnps map[string]map[string]uint8
	for chr := 1; chr <= 22; chr++ {
		fmt.Printf("---- %d ----\n", chr)
		chrStr = strconv.Itoa(chr)
		positionMap = make(map[string]struct{})
		for idv = range guessed {
			for pos = range guessed[idv][chrStr] {
				positionMap[pos] = struct{}{}
			}
		}
		positions = make([]string, 0, len(positionMap))
		for pos = range positionMap {
			positions = append(positions, pos)
		}
		relations = make(map[string][]Relation)
		for idv := range guessed {
			relations[idv] = make([]Relation, 0)
		}
		gtpFreqs = loadGenotypeFrequencies(chr)
		fmt.Println("Frequencies loaded")
		otherSnps = getIndividualPositionAlleles(chrStr, positions, mainSamples, tools.GG)
		fmt.Println("Main samples' alleles loaded")
		for _, other := range mainSamples {
			for idv := range guessed {
				relations[idv] = append(relations[idv], Relation{other,
					matchCount(guessed[idv][chrStr], otherSnps[other]),
					kingRobustPair(guessed[idv][chrStr], otherSnps[other]),
					sharedInformation(guessed[idv][chrStr], otherSnps[other], gtpFreqs)})
			}
		}
		fmt.Println("Main samples' relations calculated")

		otherSnps = getIndividualPositionAlleles(chrStr, positions, relatives, tools.RL)
		for _, relative := range relatives {
			for idv := range guessed {
				relations[idv] = append(relations[idv], Relation{relative,
					matchCount(guessed[idv][chrStr], otherSnps[relative]),
					kingRobustPair(guessed[idv][chrStr], otherSnps[relative]),
					sharedInformation(guessed[idv][chrStr], otherSnps[relative], gtpFreqs)})
			}
		}
		fmt.Println("Relatives' relations calculated")

		references = positionsMajorAlleleSamples(chrStr, positions, pgs.POPULATIONS)
		for idv := range guessed {
			relations[idv] = append(relations[idv], Relation{"reference",
				matchCount(guessed[idv][chrStr], references[getIndividualAncestry(idv, ancestries)]),
				kingRobustPair(guessed[idv][chrStr], references[getIndividualAncestry(idv, ancestries)]),
				sharedInformation(guessed[idv][chrStr], references[getIndividualAncestry(idv, ancestries)],
					gtpFreqs)})
		}
		fmt.Println("References' relations calculated")

		filepath = path.Join(resultFolder, fmt.Sprintf("unimputed%d.json", chr))
		resFile, err = os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatalf("Error opening result file: %v", err)
		}
		encoder = json.NewEncoder(resFile)
		if err = encoder.Encode(relations); err != nil {
			log.Fatal("Cannot encode json", err)
		}
		resFile.Close()
	}
	fmt.Println("Completed")
}

func loadGuessedGenotypes() map[string]map[string]map[string]uint8 {
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
	guessed := make(map[string]map[string]map[string]uint8)
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
				guessed[g.Individual] = make(map[string]map[string]uint8)
				for locus, snp := range g.SNPs {
					chr, pos := tools.SplitLocus(locus)
					if _, ok := guessed[g.Individual][chr]; !ok {
						guessed[g.Individual][chr] = make(map[string]uint8)
					}
					guessed[g.Individual][chr][pos] = snp
				}
			}
			file.Close()
		}
	}
	return guessed
}

func similarity(a, b map[string]uint8, af map[string][]float32) float32 {
	var score float32 = 0.0
	for locus := range a {
		if a[locus] == b[locus] {
			if _, ok := b[locus]; !ok {
				b[locus] = 0
			}
			switch a[locus] {
			case 0:
				score += solver.AfToLikelihood(af[locus][0]) * pgs.Ploidy
			case 1:
				score += solver.AfToLikelihood(af[locus][0]) + solver.AfToLikelihood(af[locus][1])
			case 2:
				score += solver.AfToLikelihood(af[locus][1]) * pgs.Ploidy
			default:
				log.Fatalf("Invalid SNP value: %d", a[locus])
			}
		}
	}
	return score
}

func matchCount(a, b map[string]uint8) []int {
	var numMatches int = 0
	for locus := range a {
		if _, ok := b[locus]; !ok {
			b[locus] = 0
		}
		if a[locus] == b[locus] {
			numMatches++
		}
	}
	return []int{numMatches, len(a)}
}

func sharedInformation(a, b map[string]uint8, genFreq map[string][]float32) []float32 {
	var cumul float32 = 0
	for locus := range a {
		if _, ok := b[locus]; !ok {
			b[locus] = 0
		}
		if a[locus] == b[locus] {
			if _, ok := genFreq[locus]; !ok {
				//log.Printf("Missing frequency for locus %s", locus)
				continue
			}
			cumul += -float32(math.Log2(float64(genFreq[locus][a[locus]])))
		}
	}
	return []float32{cumul, float32(len(a))}
}

func loadGenotypeFrequencies(chr int) map[string][]float32 {
	folder := "data/frequencies"
	filepath := path.Join(folder, fmt.Sprintf("%d.json", chr))
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatalf("Cannot open file %s: %v", filepath, err)
	}
	defer file.Close()
	var frequencies map[string][]float32
	decoder := json.NewDecoder(file)
	if err = decoder.Decode(&frequencies); err != nil {
		log.Fatalf("Cannot decode json file %s: %v", filepath, err)
	}
	return frequencies
}

func countSimilarity(a, b map[string]uint8) float32 {
	var score float32 = 0.0
	for locus := range a {
		if _, ok := b[locus]; !ok {
			b[locus] = 0
		}
		if a[locus] == b[locus] {
			score++
		}
	}
	return score / float32(len(a))
}

func calculateKinship(chrId string, imputedSNPs map[string]map[string]uint8) {
	related := solver.ReadRelatedIndividuals()
	var relativeSNPs map[string]uint8
	for idv := range imputedSNPs {
		for _, relative := range related[idv] {
			relativeSNPs = getIndividualAlleles(chrId, relative, tools.GG)
			fmt.Printf("Individual %s and relative %s: %.3f\n", idv, relative, kingRobust(imputedSNPs[idv], relativeSNPs))
		}
	}
}

func calculateAccuracy(chrId int, ancestry string, trueSNPs, imputedSNPs map[string]map[string]uint8) {
	fmt.Printf("Retrieving AF for chromosome %d\n", chrId)
	chrAF := getChromosomeAF(strconv.Itoa(chrId), ancestry)
	fmt.Println("Retrieval complete")
	var exists bool
	for idv := range trueSNPs {
		numMatches, totalSNPs := 0, 0
		var adjMatch, adjTotal, addition float64 = 0.0, 0.0, 0.0
		for pos := range trueSNPs[idv] {
			switch trueSNPs[idv][pos] {
			case 0:
				addition = 2 / chrAF[pos][0]
			case 1:
				addition = 1/chrAF[pos][0] + 1/chrAF[pos][1]
			case 2:
				addition = 2 / chrAF[pos][1]
			default:
				log.Printf("Invalid SNP value: %d", trueSNPs[idv][pos])
			}
			adjTotal += addition
			totalSNPs++
			if _, exists = imputedSNPs[idv][pos]; exists && trueSNPs[idv][pos] == imputedSNPs[idv][pos] {
				numMatches++
				adjMatch += addition
			}
		}
		fmt.Printf("Individual %s: %d matches out of %d, accuracy %.3f, AF-adjusted acc %.4f\n", idv, numMatches,
			totalSNPs, float64(numMatches)/float64(totalSNPs), adjMatch/adjTotal)
	}
}

func getPositionsAFs(chr string, positions []string, ancestries []string) map[string]map[string]float32 {
	dataset := tools.GG
	// Define the bcftools command
	prg := "bcftools"
	format := "%POS-%" + strings.Join(ancestries, "_AF\t%") + "_AF\n"
	args := []string{
		"query",
		"-f", format,
		tools.GetChromosomeFilepath(chr, dataset),
		"-r", "",
	}
	posAlleleFreqs := make(map[string]map[string]float32)
	for _, ancestry := range ancestries {
		posAlleleFreqs[ancestry] = make(map[string]float32)
	}
	var pos string
	var afFloat float64
	for _, position := range positions {
		args[len(args)-1] = fmt.Sprintf("%s:%s-%s", chr, position, position)
		output, err := exec.Command(prg, args...).Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]
		for _, line := range lines {
			fields := strings.Split(line, "-")
			pos = fields[0]
			if pos != position {
				continue
			}
			afs := strings.Split(fields[1], "\t")
			for i, ancestry := range ancestries {
				for _, afStr := range strings.Split(afs[i], ",") {
					afFloat, err = strconv.ParseFloat(afStr, 64)
					if err != nil {
						log.Printf("Error converting to AF: %s -- %v", line, err)
						continue
					}
					posAlleleFreqs[ancestry][pos] += float32(afFloat)
				}
			}
		}
		for _, ancestry := range ancestries {
			if posAlleleFreqs[ancestry][pos] == 0 {
				posAlleleFreqs[ancestry][pos] = pgs.MissingEAF
			}
			if posAlleleFreqs[ancestry][pos] >= 1 {
				fmt.Printf("AF for %s at position %s: %.3f\n", ancestry, pos, posAlleleFreqs[ancestry][pos])
				posAlleleFreqs[ancestry][pos] = 1 - pgs.MissingEAF
			}
		}
	}
	return posAlleleFreqs
}

func getChromosomeAF(chr, ancestry string) map[string][]float64 {
	dataset := tools.GG

	// Define the bcftools command
	prg := "bcftools"
	args := []string{
		"query",
		"-f", fmt.Sprintf("%%POS\t%%%s_AF\n", ancestry),
		"-r", chr,
		tools.GetChromosomeFilepath(chr, dataset),
	}
	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command: %v", err)
	}
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]
	af := make(map[string][]float64)
	var afTotal float64
	for _, line := range lines {
		fields := strings.Split(line, "\t")
		pos := fields[0]
		afs := strings.Split(fields[1], ",")
		afTotal = 0
		for _, afStr := range afs {
			afFloat, err := strconv.ParseFloat(afStr, 64)
			if err != nil {
				log.Printf("Error converting to AF: %s -- %v", line, err)
				continue
			}
			afTotal += afFloat
		}
		if afTotal == 0 {
			afTotal = pgs.MissingEAF
		}
		if afTotal > 1 {
			log.Printf("AF total is greater than 1: %s\n", line)
		}
		af[pos] = []float64{1 - afTotal, afTotal}
	}
	return af
}

func getChromosomeAFs(chr string) map[string]map[string][]float32 {
	dataset := tools.GG
	prg := "bcftools"
	args := []string{
		"query",
		"-f", fmt.Sprintf("%%POS-%%" + strings.Join(pgs.POPULATIONS, "_AF\\t%%") + "_AF\n"),
		//"-f", fmt.Sprintf("%%POS\t%%%s_AF\n", ancestry),
		"-r", chr,
		tools.GetChromosomeFilepath(chr, dataset),
	}
	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command %v: %v", args, err)
	}
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]
	af := make(map[string]map[string][]float32)
	for _, ppl := range pgs.POPULATIONS {
		af[ppl] = make(map[string][]float32)
	}
	var afTotal float32
	var pos string
	var fields, pplAFs []string
	for _, line := range lines {
		fields = strings.Split(line, "-")
		pos = fields[0]
		pplAFs = strings.Split(fields[1], "\t")
		for i, ppl := range pgs.POPULATIONS {
			afs := strings.Split(pplAFs[i], ",")
			afTotal = 0
			for _, afStr := range afs {
				afFloat, err := strconv.ParseFloat(afStr, 64)
				if err != nil {
					log.Printf("Error converting to AF: %s -- %v", line, err)
					continue
				}
				afTotal += float32(afFloat)
			}
			if afTotal == 0 {
				afTotal = pgs.MissingEAF
			}
			if afTotal > 1 {
				log.Printf("AF total is greater than 1: %s\n", line)
				afTotal = 1
			}
			af[ppl][pos] = []float32{1 - afTotal, afTotal}
		}
	}
	return af
}

func getIndividualPositionAlleles(chr string, positions []string, individuals []string, dataset string) map[string]map[string]uint8 {
	prg := "bcftools"
	args := []string{
		"query",
		"-s", strings.Join(individuals, ","),
		"-f", "%POS-[%SAMPLE=%GT\t]\n",
		tools.GetChromosomeFilepath(chr, dataset),
		"-r", "",
	}
	var pos, gt, idv string
	var ok bool
	var snp uint8
	alleles := make(map[string]map[string]uint8)
	for _, position := range positions {
		args[len(args)-1] = fmt.Sprintf("%s:%s-%s", chr, position, position)
		output, err := exec.Command(prg, args...).Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]
		for _, line := range lines {
			fields := strings.Split(line, "-")
			if len(fields) < 2 {
				log.Printf("Not enough fields in line: %s\n", line)
				continue
			}
			pos = fields[0]
			if pos != position {
				//log.Printf("Position mismatch: %s != %s\n", pos, position)
				continue
			}
			samples := strings.Split(fields[1], "\t")
			samples = samples[:len(samples)-1]
			for _, sample := range samples {
				fields = strings.Split(sample, "=")
				idv = fields[0]
				if _, ok = alleles[idv]; !ok {
					alleles[idv] = make(map[string]uint8)
				}
				gt, err = tools.NormalizeSnp(fields[1])
				if err != nil {
					log.Printf("Error normalizing SNP: %s -- %v\n", line, err)
					continue
				}
				snp, err = tools.SnpToSum(gt)
				if err != nil {
					log.Printf("Error converting SNP to sum: %s -- %v\n", line, err)
				}
				alleles[idv][pos] = snp
			}
		}
	}
	return alleles
}

func getIndividualAlleles(chr, idv, dataset string) map[string]uint8 {
	prg := "bcftools"
	args := []string{
		"query",
		"-s", idv,
		"-f", "%POS\t[%GT]\n",
		"-r", chr,
		tools.GetChromosomeFilepath(chr, dataset),
	}
	cmd := exec.Command(prg, args...)
	output, err := cmd.Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command: %v", err)
	}
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]
	alleles := make(map[string]uint8)
	var pos, gt string
	var snp uint8
	for _, line := range lines {
		fields := strings.Split(line, "\t")
		if len(fields) < 2 {
			log.Printf("Not enough fields in line: %s\n", line)
			continue
		}
		pos = fields[0]
		gt, err = tools.NormalizeSnp(fields[1])
		if err != nil {
			log.Printf("Error normalizing SNP: %s -- %v\n", line, err)
			continue
		}
		snp, err = tools.SnpToSum(gt)
		if err != nil {
			log.Printf("Error converting SNP to sum: %s -- %v\n", line, err)
		}
		alleles[pos] = snp
	}
	return alleles
}

func getRegionIndividualAlleles(chr, filepath string, individuals []string, includedPositions []int,
	imputedRegions [][]string, numThreads int) map[string]map[string]uint8 {
	alleles := make(map[string]map[string]uint8)
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS-[%SAMPLE=%GT\t]\n",
		filepath,
	}
	if len(individuals) > 0 {
		args = append(args, "-s")
		args = append(args, strings.Join(individuals, ","))
	}
	args = append(args, "-r")
	var wg sync.WaitGroup
	var mu sync.Mutex
	var regArg string
	fmt.Printf("Receiving bcftools output: \n")
	for r, region := range imputedRegions {
		regArg = fmt.Sprintf("%s:%s-%s", chr, region[0], region[1])
		args[len(args)-1] = regArg
		output, err := exec.Command(prg, args...).Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		fmt.Printf("%d/%d\t", r+1, len(imputedRegions))
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]

		lineChan := make(chan string, numThreads)
		for i := 0; i < numThreads; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				var pos, gt, idv string
				var snp uint8
				var ok bool
				var samples []string
				for line := range lineChan {
					fields := strings.Split(line, "-")
					if len(fields) < 2 {
						log.Printf("not enough fields in line: %s", line)
						return
					}
					pos = fields[0]
					samples = strings.Split(fields[1], "\t")
					samples = samples[:len(samples)-1]
					for _, sample := range samples {
						fields = strings.Split(sample, "=")
						idv = fields[0]

						mu.Lock()
						if _, ok = alleles[idv]; !ok {
							alleles[idv] = make(map[string]uint8)
						}
						mu.Unlock()

						gt = fields[1]
						var err error
						gt, err = tools.NormalizeSnp(gt)
						if err != nil {
							log.Printf("error normalizing SNP: %s -- %v", line, err)
							continue
						}
						snp, err = tools.SnpToSum(gt)
						if err != nil {
							log.Printf("error converting SNP to sum: %s -- %v", line, err)
							continue
						}

						mu.Lock()
						alleles[idv][pos] = snp
						mu.Unlock()
					}
				}
			}()
		}

		go func() {
			for _, line := range lines {
				lineChan <- line
			}
			close(lineChan)
		}()
		wg.Wait()
	}
	fmt.Printf("\nBcftools output processed\n")
	if len(includedPositions) > 0 {
		var ok bool
		for pos := range includedPositions {
			for idv := range alleles {
				if _, ok = alleles[idv][strconv.Itoa(pos)]; ok {
					break
				}
				alleles[idv][strconv.Itoa(pos)] = 0
			}

		}
	}
	fmt.Printf("Missing loci filled\n")
	return alleles
}

func getRegionAllelesPerLocus(chr, filepath string, imputedRegions [][]string, numThreads int) map[string][]uint8 {
	alleles := make(map[string][]uint8)
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS-[%GT\t]\n",
		filepath,
		"-r", "",
	}
	fmt.Printf("Receiving bcftools output: \n")
	var wg sync.WaitGroup
	var mu sync.Mutex
	var regArg string
	for r, region := range imputedRegions {
		regArg = fmt.Sprintf("%s:%s-%s", chr, region[0], region[1])
		args[len(args)-1] = regArg
		output, err := exec.Command(prg, args...).Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		fmt.Printf("%d/%d\t", r+1, len(imputedRegions))
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]
		lineChan := make(chan string, numThreads)
		for i := 0; i < numThreads; i++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				var pos, gtp string
				var snp uint8
				var ok bool
				var samples []string
				for line := range lineChan {
					fields := strings.Split(line, "-")
					if len(fields) < 2 {
						log.Printf("not enough fields in line: %s", line)
						return
					}
					pos = fields[0]
					mu.Lock()
					if _, ok = alleles[pos]; ok {
						alleles[pos] = make([]uint8, 0)
					}
					mu.Unlock()
					samples = strings.Split(fields[1], "\t")
					samples = samples[:len(samples)-1]
					snps := make([]uint8, len(samples))
					for j, sample := range samples {
						var err error
						gtp, err = tools.NormalizeSnp(sample)
						if err != nil {
							log.Printf("error normalizing SNP: %s -- %v", line, err)
							continue
						}
						snp, err = tools.SnpToSum(gtp)
						if err != nil {
							log.Printf("error converting SNP to sum: %s -- %v", line, err)
							continue
						}
						snps[j] = snp
					}
					mu.Lock()
					alleles[pos] = append(alleles[pos], snps...)
					mu.Unlock()
				}
			}()
		}

		go func() {
			for _, line := range lines {
				lineChan <- line
			}
			close(lineChan)
		}()
		wg.Wait()
	}
	fmt.Printf("\nBcftools output processed\n")
	return alleles
}

func getAllIndividualAlleles(filepath string, numThreads int) map[string]map[string]uint8 {
	alleles := make(map[string]map[string]uint8)
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS-[%SAMPLE=%GT\t]\n",
		filepath,
	}
	var wg sync.WaitGroup
	var mu sync.Mutex
	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command: %v", err)
	}
	fmt.Printf("Bcftools output received\n")
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]

	lineChan := make(chan string, numThreads)
	for i := 0; i < numThreads; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			var pos, gt, idv string
			var snp uint8
			var ok bool
			var samples []string
			for line := range lineChan {
				fields := strings.Split(line, "-")
				if len(fields) < 2 {
					log.Printf("not enough fields in line: %s", line)
					return
				}
				pos = fields[0]
				samples = strings.Split(fields[1], "\t")
				samples = samples[:len(samples)-1]
				for _, sample := range samples {
					fields = strings.Split(sample, "=")
					idv = fields[0]

					mu.Lock()
					if _, ok = alleles[idv]; !ok {
						alleles[idv] = make(map[string]uint8)
					}
					mu.Unlock()

					gt = fields[1]
					var err error
					gt, err = tools.NormalizeSnp(gt)
					if err != nil {
						log.Printf("error normalizing SNP: %s -- %v", line, err)
						continue
					}
					snp, err = tools.SnpToSum(gt)
					if err != nil {
						log.Printf("error converting SNP to sum: %s -- %v", line, err)
						continue
					}

					mu.Lock()
					alleles[idv][pos] = snp
					mu.Unlock()
				}
			}
		}()
	}

	go func() {
		for _, line := range lines {
			lineChan <- line
		}
		close(lineChan)
	}()
	wg.Wait()
	fmt.Printf("Bcftools output processed\n")
	return alleles
}

func allMajorAlleleSample(chrId int, ancestry string) map[string]uint8 {
	af := getChromosomeAF(strconv.Itoa(chrId), ancestry)
	alleles := make(map[string]uint8)
	for pos := range af {
		if af[pos][0] > af[pos][1] {
			alleles[pos] = 0
		} else {
			alleles[pos] = 2
		}
	}
	return alleles
}

func positionsMajorAlleleSamples(chr string, positions []string, ancestries []string) map[string]map[string]uint8 {
	af := getPositionsAFs(chr, positions, ancestries)
	alleles := make(map[string]map[string]uint8)
	for _, ancestry := range ancestries {
		alleles[ancestry] = make(map[string]uint8)
		for pos := range af[ancestry] {
			if af[ancestry][pos] >= 0.5 {
				alleles[ancestry][pos] = 2
			} else {
				alleles[ancestry][pos] = 0
			}
		}
	}
	return alleles
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

func getAllPositions(filepath string) []string {
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS\n",
		filepath,
	}
	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command: %v", err)
	}
	lines := strings.Split(string(output), "\n")
	lines = lines[:len(lines)-1]
	return lines
}
