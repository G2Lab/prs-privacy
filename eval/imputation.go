package main

import (
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
	postImputeTestingFilename  = "imputation/testing/chr%d-%s-postimpute.vcf.gz"
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
			//fillTruthVCF(chrId, ancestry)
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

func fillTruthVCF(chrId int, ancestry string) {
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
	stringedPositions := make([]string, len(chrPositions))
	for i, pos := range chrPositions {
		stringedPositions[i] = strconv.Itoa(pos)
	}

	trueGtp := getIndividualPositionSNPs(chr, stringedPositions, individuals, tools.RL)

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
			if _, ok = trueGtp[idv]; !ok {
				line += "\t./."
				continue
			}
			if snp, ok = trueGtp[idv][stringedPositions[i]]; ok {
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
				log.Printf("No SNP for individual %s at position %s\n", idv, stringedPositions[i])
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

type Relation struct {
	Target      string
	Count       []int
	King        []int
	Information []float32
}

func newRelation(target string, imputed, compared map[string]uint8, freqs map[string][]float32) *Relation {
	return &Relation{
		Target:      target,
		Count:       matchCount(imputed, compared),
		King:        kingRobustPair(imputed, compared),
		Information: mutualInformation(imputed, compared, freqs),
	}
}

func linkingWithImputation() {
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
	relations := make(map[string][]*Relation)
	for idv := range imputedSNPs {
		relations[idv] = make([]*Relation, 0)
	}

	// Read all the imputed positions, sort them and divide them into regions
	imputedPositions := getAllImputedPositions(chrId)
	imputedRegions := dividePositionsIntoRegions(imputedPositions, imputationChunkSize)
	fmt.Printf("Number of imputed SNPs: %d\n", len(imputedPositions))
	fmt.Printf("Imputation regions: %v\n", imputedRegions)
	gtpFrequencies := loadGenotypeFrequencies(chrId)

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
	chunkSize := 313

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for t := range taskChan {
				rel := newRelation(t.other, imputedSNPs[t.idv], t.snps[t.other], gtpFrequencies)
				mu.Lock()
				relations[t.idv] = append(relations[t.idv], rel)
				mu.Unlock()
				taskWg.Done()
			}
		}()
	}
	go func() {
		for i := 0; i < len(db); i += chunkSize {
			fmt.Printf("==== Chunk %d ====\n", i/chunkSize)
			chunk := db[i:min(i+chunkSize, len(db))]
			chrPath := tools.GetChromosomeFilepath(strconv.Itoa(chrId), tools.GG)
			chunkSnps := getRegionIndividualSNPs(strconv.Itoa(chrId), chrPath, chunk, imputedPositions,
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

	stringedImputedPositions := make([]string, len(imputedPositions))
	for i, pos := range imputedPositions {
		stringedImputedPositions[i] = strconv.Itoa(pos)
	}
	chrStr := strconv.Itoa(chrId)
	ancestries := tools.LoadAncestry()
	references := positionsMajorAlleleSamples(chrStr, imputedRegions, pgs.POPULATIONS)

	relatives := solver.AllRelativeSamples()
	relativeSnps := getRegionIndividualSNPs(chrStr, tools.GetChromosomeFilepath(chrStr, tools.RL), relatives,
		imputedPositions, imputedRegions, numWorkers)
	for idv := range imputedSNPs {
		for other := range relativeSnps {
			relations[idv] = append(relations[idv], newRelation(other, imputedSNPs[idv], relativeSnps[other], gtpFrequencies))
		}
		relations[idv] = append(relations[idv], newRelation("reference", imputedSNPs[idv],
			references[getIndividualAncestry(idv, ancestries)], gtpFrequencies))
	}
	resultFolder := "results/linking"
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

func linkingWithGuessed() {
	guessed := loadGuessedGenotypes()
	relatives := solver.AllRelativeSamples()
	mainSamples := solver.All1000GenomesSamples()
	ancestries := tools.LoadAncestry()
	resultFolder := "results/linking"
	var encoder *json.Encoder
	var resFile *os.File
	var err error
	var j int
	var chrStr, idv, pos, filepath string
	var positionMap map[string]struct{}
	var positions []string
	var regions [][]string
	var relations map[string][]*Relation
	var gtpFrequencies map[string][]float32
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
		relations = make(map[string][]*Relation)
		for idv := range guessed {
			relations[idv] = make([]*Relation, 0)
		}
		gtpFrequencies = loadGenotypeFrequencies(chr)
		fmt.Println("Frequencies loaded")
		otherSnps = getIndividualPositionSNPs(chrStr, positions, mainSamples, tools.GG)
		fmt.Println("Main samples' alleles loaded")
		for _, other := range mainSamples {
			for idv := range guessed {
				relations[idv] = append(relations[idv], newRelation(other, guessed[idv][chrStr], otherSnps[other], gtpFrequencies))
			}
		}
		fmt.Println("Main samples' relations calculated")

		otherSnps = getIndividualPositionSNPs(chrStr, positions, relatives, tools.RL)
		for _, relative := range relatives {
			for idv := range guessed {
				relations[idv] = append(relations[idv], newRelation(relative, guessed[idv][chrStr], otherSnps[relative], gtpFrequencies))
			}
		}
		fmt.Println("Relatives' relations calculated")

		regions = make([][]string, len(positions))
		for j, pos = range positions {
			regions[j] = []string{pos, pos}
		}
		references = positionsMajorAlleleSamples(chrStr, regions, pgs.POPULATIONS)
		for idv := range guessed {
			relations[idv] = append(relations[idv], newRelation("reference", guessed[idv][chrStr],
				references[getIndividualAncestry(idv, ancestries)], gtpFrequencies))
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

func evaluateImputation() {
	var chrId int
	var err error
	//if len(os.Args) > 1 {
	//	chrId, err = strconv.Atoi(os.Args[1])
	//	if err != nil {
	//		log.Fatalf("Cannot convert chromosome ID to integer: %v", err)
	//	}
	//} else {
	//	chrId = 22
	//}
	for chrId = 1; chrId <= 22; chrId++ {
		chrStr := strconv.Itoa(chrId)
		fmt.Printf("Imputation evaluation for chromosome %d\n", chrId)
		numWorkers := 4
		imputedSNPs := make(map[string]map[string]uint8)
		for _, ancestry := range pgs.POPULATIONS {
			fmt.Printf("Reading imputed SNPs: %d, %s\n", chrId, ancestry)
			chrPath := fmt.Sprintf(postImputeTestingFilename, chrId, ancestry)
			imputedChunk := getAllIndividualAlleles(chrPath, numWorkers)
			for idv := range imputedChunk {
				imputedSNPs[idv] = imputedChunk[idv]
			}
		}
		individuals := make([]string, 0)
		for idv := range imputedSNPs {
			individuals = append(individuals, idv)
		}

		guessed := loadGuessedGenotypes()
		guessedPositionsMap := make(map[int]struct{})
		var posInt int
		for idv := range guessed {
			for chr := range guessed[idv] {
				if chr != strconv.Itoa(chrId) {
					continue
				}
				for pos := range guessed[idv][chr] {
					posInt, err = strconv.Atoi(pos)
					if err != nil {
						fmt.Printf("Cannot convert position %s to integer: %v\n", pos, err)
					}
					guessedPositionsMap[posInt] = struct{}{}
				}
			}
		}

		// Read all the imputed positions, sort them and divide them into regions
		imputedPositions := getAllImputedPositions(chrId)
		imputedRegions := dividePositionsIntoRegions(imputedPositions, imputationChunkSize)
		fmt.Printf("Number of imputed SNPs: %d\n", len(imputedPositions))
		fmt.Printf("Imputation regions: %v\n", imputedRegions)

		trueSNPs := getRegionIndividualSNPs(chrStr, tools.GetChromosomeFilepath(chrStr, tools.RL), individuals,
			imputedPositions, imputedRegions, numWorkers)

		gtpFrequencies := loadGenotypeFrequencies(chrId)
		observations := make(map[int]float32)
		var wg sync.WaitGroup
		type result struct {
			pos   int
			value float32
		}
		tasks := make(chan int, numWorkers)
		results := make(chan result, len(imputedPositions))
		// Worker function
		worker := func(tasks <-chan int, results chan<- result) {
			var strPos, idv string
			var ok bool
			var majorGtp uint8
			for t := range tasks {
				strPos = strconv.Itoa(t)
				confusionMatrix := make([][]int, 3)
				for i := 0; i < 3; i++ {
					confusionMatrix[i] = make([]int, 3)
				}

				for idv = range imputedSNPs {
					if _, ok = trueSNPs[idv][strPos]; ok {
						confusionMatrix[trueSNPs[idv][strPos]][imputedSNPs[idv][strPos]]++
					}
				}

				// If all the true genotypes are major, skip
				majorGtp = majorGenotype(gtpFrequencies[strPos])
				if confusionMatrix[majorGtp][majorGtp] == len(trueSNPs) {
					continue
				}
				results <- result{pos: t, value: macroAveragedF1(confusionMatrix)}
			}
		}

		for w := 0; w < numWorkers; w++ {
			wg.Add(1)
			go func() {
				defer wg.Done()
				worker(tasks, results)
			}()
		}
		go func() {
			//chunkSize := len(imputedPositions) / 100
			//for i, pos := range imputedPositions {
			for _, pos := range imputedPositions {
				tasks <- pos
				//if i%chunkSize == 0 {
				//	fmt.Printf("Progress: %d%%\n", i/chunkSize)
				//}
			}
			close(tasks)
			wg.Wait()
			close(results)
		}()

		for res := range results {
			observations[res.pos] = res.value
		}
		fmt.Println("Observations calculated")

		accuracyByDistance := make(map[int][]float32)
		var regionGuessedPositions []int
		var regionStart, regionEnd int
		var closestDistance int
		for _, region := range imputedRegions {
			regionGuessedPositions = make([]int, 0)
			regionStart, err = strconv.Atoi(region[0])
			if err != nil {
				log.Fatalf("Cannot convert region start to integer: %v", err)
			}
			regionEnd, err = strconv.Atoi(region[1])
			if err != nil {
				log.Fatalf("Cannot convert region end to integer: %v", err)
			}
			for pos := range guessedPositionsMap {
				if regionStart <= pos && pos <= regionEnd {
					regionGuessedPositions = append(regionGuessedPositions, pos)
				}
			}
			sort.Ints(regionGuessedPositions)

			for _, pos := range imputedPositions {
				if regionStart > pos {
					continue
				}
				if pos > regionEnd {
					break
				}
				closestDistance = math.MaxInt32
				for _, guessedPos := range regionGuessedPositions {
					if dist := abs(pos - guessedPos); dist < closestDistance {
						closestDistance = dist
					}
				}
				if _, ok := accuracyByDistance[closestDistance]; !ok {
					accuracyByDistance[closestDistance] = make([]float32, 0)
				}
				accuracyByDistance[closestDistance] = append(accuracyByDistance[closestDistance], observations[pos])
			}
		}

		// Save results
		filepath := path.Join("results/impute", fmt.Sprintf("%d.json", chrId))
		resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
		if err != nil {
			log.Fatalf("Error opening result file: %v", err)
		}
		defer resFile.Close()
		encoder := json.NewEncoder(resFile)
		if err = encoder.Encode(accuracyByDistance); err != nil {
			log.Fatal("Cannot encode json", err)
		}
	}
	fmt.Println("Completed")
}

func majorGenotype(freq []float32) uint8 {
	var major uint8 = 0
	for i := 1; i < len(freq); i++ {
		if freq[i] > freq[major] {
			major = uint8(i)
		}
	}
	return major
}

func abs(a int) int {
	if a >= 0 {
		return a
	}
	return -a
}

func macroAveragedF1(confusionMatrix [][]int) float32 {
	precision := make([]float32, len(confusionMatrix))
	recall := make([]float32, len(confusionMatrix))
	f1Scores := make([]float32, len(confusionMatrix))

	var tp, fp, fn float32
	// Calculate precision and recall for each class
	for i := 0; i < len(confusionMatrix); i++ {
		tp = float32(confusionMatrix[i][i])
		fp = 0
		fn = 0
		for j := 0; j < len(confusionMatrix[i]); j++ {
			if j != i {
				fp += float32(confusionMatrix[j][i])
				fn += float32(confusionMatrix[i][j])
			}
		}
		if tp+fp != 0 {
			precision[i] = tp / (tp + fp)
		} else {
			// If there are no instances of the class and none is recognized as such, the precision is 1
			precision[i] = 1.0
		}
		if tp+fn != 0 {
			recall[i] = tp / (tp + fn)
		} else {
			recall[i] = 1.0
		}
		if precision[i]+recall[i] != 0 {
			f1Scores[i] = 2 * (precision[i] * recall[i]) / (precision[i] + recall[i])
		} else {
			f1Scores[i] = 0.0
		}
	}

	// Calculate the macro-averaged F1 score
	var macroF1 float32
	for _, f1 := range f1Scores {
		macroF1 += f1
	}
	//for _, pr := range precision {
	//	macroF1 += pr
	//}
	macroF1 /= 3

	return macroF1
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

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func calculateGenotypeFrequencies() {
	numWorkers := 5
	missingGenotypeFrequencies := []float32{(1 - pgs.MissingEAF) * (1 - pgs.MissingEAF),
		2 * pgs.MissingEAF * (1 - pgs.MissingEAF), pgs.MissingEAF * pgs.MissingEAF}
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
				switch frequencies[locus][i] {
				case 0:
					frequencies[locus][i] = 1 / float32(len(alleles[locus]))
				case 1:
					frequencies[locus][i] = 1 - 1/float32(len(alleles[locus]))
				default:

				}
			}
		}
		for _, pos := range positions {
			if _, ok = frequencies[strconv.Itoa(pos)]; !ok {
				fmt.Printf("Position %d was missing\n", pos)
				frequencies[strconv.Itoa(pos)] = missingGenotypeFrequencies
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
		//imputedPositionsPerAncestry = getAllPositions(fmt.Sprintf(postImputeZippedFilename, chr, ancestry))
		imputedPositionsPerAncestry = getAllPositions(fmt.Sprintf(postImputeTestingFilename, chr, ancestry))
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

func mutualInformation(a, b map[string]uint8, genFreq map[string][]float32) []float32 {
	var cumul float64 = 0
	for locus := range a {
		if _, ok := b[locus]; !ok {
			b[locus] = 0
		}
		if a[locus] == b[locus] {
			if _, ok := genFreq[locus]; !ok {
				//log.Printf("Missing frequency for locus %s", locus)
				continue
			}
			cumul += -math.Log2(float64(genFreq[locus][a[locus]]))
		}
	}
	if math.IsInf(cumul, 1) || math.IsInf(cumul, -1) {
		cumul = math.MaxFloat32
	}
	return []float32{float32(cumul), float32(len(a))}
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
			relativeSNPs = getAllIndividualSNPs(chrId, relative, tools.GG)
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

func getPositionsAFs(chr string, regions [][]string, ancestries []string) map[string]map[string]float32 {
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
	fmt.Printf("Receiving bcftools output: \n")
	for r, region := range regions {
		args[len(args)-1] = fmt.Sprintf("%s:%s-%s", chr, region[0], region[1])
		output, err := exec.Command(prg, args...).Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		fmt.Printf("%d/%d\t", r+1, len(regions))
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]
		for _, line := range lines {
			fields := strings.Split(line, "-")
			pos = fields[0]
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
				//fmt.Printf("AF for %s at position %s: %.3f\n", ancestry, pos, posAlleleFreqs[ancestry][pos])
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

func getIndividualPositionSNPs(chr string, positions []string, individuals []string, dataset string) map[string]map[string]uint8 {
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

func getAllIndividualSNPs(chr, idv, dataset string) map[string]uint8 {
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

func getRegionIndividualSNPs(chr, filepath string, individuals []string, includedPositions []int,
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
	args = append(args, "")
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

func positionsMajorAlleleSamples(chr string, regions [][]string, ancestries []string) map[string]map[string]uint8 {
	af := getPositionsAFs(chr, regions, ancestries)
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
