package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
	"log"
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

func imputeWorkflow() {
	//for _, chrId := range []int{22} {
	for chrId := 1; chrId <= 21; chrId++ {
		for _, ancestry := range pgs.POPULATIONS {
			//for _, ancestry := range []string{"EUR"} {
			fmt.Printf("===== %s =====\n", ancestry)
			fillPreImputeVCF(chrId, ancestry)
			fmt.Printf("Pre-impute VCF filled\n")
			compressVCF(chrId, ancestry)
			fmt.Printf("VCF compressed\n")
			indexCompressedVCF(chrId, ancestry, false)
			fmt.Printf("VCF indexed\n")
			impute(ancestry, chrId)
			fmt.Printf("Imputation complete\n")
			indexCompressedVCF(chrId, ancestry, true)
			fmt.Printf("Imputed VCF indexed\n")
		}
	}
}

func impute(ancestry string, chrId int) {
	prg := "../minimac4/bin/minimac4"
	args := []string{
		fmt.Sprintf("data/references/%d.msav", chrId),
		"-t", "4",
		"-a",
		"--sample-ids-file", fmt.Sprintf("data/1000g_%s_no_relatives.txt", ancestry),
		"--min-ratio", "0",
		"--format", "GT",
		fmt.Sprintf(preImputeZippedFilename, chrId, ancestry),
		//"-o", fmt.Sprintf(postImputeUnzippedFilename, chrId, ancestry),
		//"-O", "vcf.gz",
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
	guessed := make(map[string]map[string]uint8)
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
				guessed[g.Individual] = g.SNPs
			}
			file.Close()
		}
	}

	var chr, pos, ref, alt string
	var posInt int
	individuals := make([]string, 0)
	chrPositions := make([]int, 0)
	alleles := make(map[int][]string)
	for idv := range guessed {
		if ancestry == getIndividualAncestry(idv, ancestries) {
			individuals = append(individuals, idv)
		}
	}
	for locus := range guessed[individuals[0]] {
		chr, pos = tools.SplitLocus(locus)
		if chr != strconv.Itoa(chrId) {
			continue
		}
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

	formatHeader := "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
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
			if snp, ok = guessed[idv][fmt.Sprintf("%d:%d", chrId, chrPositions[i])]; ok {
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
				log.Printf("No SNP for individual %s at locus %s\n", idv, fmt.Sprintf("%d:%d", chrId, chrPositions[i]))
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

func imputationAccuracy(chrId int) {
	type relation struct {
		target string
		count  float32
		king   float32
	}
	type Result struct {
		Self          map[string]int // count, king
		Relative      map[string]int // count, king
		CountAccuracy map[string]float32
		KingScore     map[string]float32
	}
	newResult := func() *Result {
		return &Result{
			Self:          make(map[string]int),
			Relative:      make(map[string]int),
			CountAccuracy: make(map[string]float32),
			KingScore:     make(map[string]float32),
		}
	}

	imputedSNPs := make(map[string]map[string]uint8)
	for _, ancestry := range pgs.POPULATIONS {
		fmt.Printf("Reading imputed SNPs: %d, %s\n", chrId, ancestry)
		imputedChunk := getAllIndividualAlleles(fmt.Sprintf(postImputeZippedFilename, chrId, ancestry), nil)
		for idv := range imputedChunk {
			imputedSNPs[idv] = imputedChunk[idv]
		}
	}
	individuals := make([]string, 0)
	for idv := range imputedSNPs {
		individuals = append(individuals, idv)
	}
	db := solver.All1000GenomesSamples()
	relations := make(map[string][]relation)
	for idv := range imputedSNPs {
		relations[idv] = make([]relation, 0)
	}

	var mu sync.Mutex
	var wg sync.WaitGroup
	// Channel to send work to goroutines
	type task struct {
		idv   string
		other string
		snps  map[string]map[string]uint8
	}
	taskChan := make(chan task)

	// Start 11 worker goroutines
	numWorkers := 15
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for t := range taskChan {
				similarity := countSimilarity(imputedSNPs[t.idv], t.snps[t.other])
				kingValue := kingRobust(imputedSNPs[t.idv], t.snps[t.other])
				rel := relation{t.other, similarity, kingValue}

				mu.Lock()
				relations[t.idv] = append(relations[t.idv], rel)
				mu.Unlock()
			}
		}()
	}
	go func() {
		chunkSize := 250
		for i := 0; i < len(db); i += chunkSize {
			fmt.Printf("==== Chunk %d ====\n", i/chunkSize)
			chunk := db[i:min(i+chunkSize, len(db))]
			chunkSnps := getAllIndividualAlleles(tools.GetChromosomeFilepath(strconv.Itoa(chrId), tools.GG), chunk)
			fmt.Printf("Chunk loaded\n")
			for idv := range imputedSNPs {
				for other := range chunkSnps {
					taskChan <- task{idv, other, chunkSnps}
				}
			}
		}
		close(taskChan)
	}()
	wg.Wait()
	//for i := 0; i < len(db); i += chunkSize {
	//	fmt.Printf("==== Chunk %d ====\n", i/chunkSize)
	//	chunk := db[i:min(i+chunkSize, len(db))]
	//	chunkSnps := getAllIndividualAlleles(tools.GetChromosomeFilepath(strconv.Itoa(chrId), tools.GG), chunk)
	//	fmt.Printf("Chunk loaded\n")
	//	for idv := range imputedSNPs {
	//		for other := range chunkSnps {
	//			relations[idv] = append(relations[idv], relation{other,
	//				countSimilarity(imputedSNPs[idv], chunkSnps[other]), kingRobust(imputedSNPs[idv], chunkSnps[other])})
	//		}
	//	}
	//}

	relatives := solver.AllRelativeSamples()
	relativeSnps := getAllIndividualAlleles(tools.GetChromosomeFilepath(strconv.Itoa(chrId), tools.RL), relatives)
	for idv := range imputedSNPs {
		for other := range relativeSnps {
			relations[idv] = append(relations[idv], relation{other,
				countSimilarity(imputedSNPs[idv], relativeSnps[other]), kingRobust(imputedSNPs[idv], relativeSnps[other])})
		}
	}

	related := solver.ReadRelatedIndividuals()
	var selfFound, relativeFound bool
	results := make(map[string]*Result)
	for idv := range imputedSNPs {
		fmt.Printf("==== %s ====\n", idv)
		results[idv] = newResult()
		for _, metric := range []string{"count", "king"} {
			fmt.Printf("---- Metric: %s\n", metric)
			sort.Slice(relations[idv], func(i, j int) bool {
				if metric == "count" {
					return relations[idv][i].count > relations[idv][j].count
				}
				return relations[idv][i].king > relations[idv][j].king
			})
			relatives = related[idv]
			selfFound, relativeFound = false, false
			for i, rel := range relations[idv] {
				if rel.target == idv {
					fmt.Printf("Self: %d pos, %.0f count score, %.0f king score\n", i, rel.count, rel.king)
					results[idv].Self[metric] = i
					results[idv].CountAccuracy["self"] = rel.count
					selfFound = true
				}
				for _, rltv := range relatives {
					if rel.target == rltv {
						fmt.Printf("Relative %s: %d pos, %.0f count score, %.0f king score\n", rltv, i, rel.count, rel.king)
						results[idv].Relative[metric] = i
						results[idv].KingScore["relative"] = rel.king
						relativeFound = true
						break
					}
				}
				if selfFound && relativeFound {
					break
				}
			}
		}
	}
	resultFolder := "results/impute"
	filepath := path.Join(resultFolder, "imputed22.json")
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

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
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
			//relations[idv] = append(relations[idv], relation{relative, similarity(imputedSNPs[idv], otherSnps, chrAF[ancestry[idv]])})
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

	type relation struct {
		target string
		score  float32
	}
	chrId := "22"
	related := solver.ReadRelatedIndividuals()
	relations := make(map[string][]relation)
	for idv := range guessed {
		relations[idv] = make([]relation, 0)
	}

	positions := make([]string, 0)
	for idv := range guessed {
		for pos := range guessed[idv][chrId] {
			positions = append(positions, pos)
		}
		break
	}
	for _, other := range solver.All1000GenomesSamples() {
		fmt.Printf("Main: %s\n", other)
		otherSnps := getIndividualPositionAlleles(chrId, other, tools.GG, positions)
		for idv := range guessed {
			relations[idv] = append(relations[idv], relation{other, countSimilarity(guessed[idv][chrId], otherSnps)})
		}
	}

	relatives := solver.AllRelativeSamples()
	for _, relative := range relatives {
		fmt.Printf("Extra: %s\n", relative)
		otherSnps := getIndividualPositionAlleles(chrId, relative, tools.RL, positions)
		for idv := range guessed {
			relations[idv] = append(relations[idv], relation{relative, countSimilarity(guessed[idv][chrId], otherSnps)})
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
	for idv := range guessed {
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
	filepath := path.Join(resultFolder, "unimputed.json")
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

func getIndividualPositionAlleles(chr, idv, dataset string, positions []string) map[string]uint8 {
	prg := "bcftools"
	var args []string
	var pos, gt string
	var snp uint8
	alleles := make(map[string]uint8)
	for _, position := range positions {
		args = []string{
			"query",
			"-s", idv,
			"-f", "%POS\t[%GT]\n",
			"-r", fmt.Sprintf("%s:%s-%s", chr, position, position),
			tools.GetChromosomeFilepath(chr, dataset),
		}
		cmd := exec.Command(prg, args...)
		output, err := cmd.Output()
		if err != nil {
			log.Fatalf("Error executing bcftools command: %v", err)
		}
		lines := strings.Split(string(output), "\n")
		lines = lines[:len(lines)-1]
		for _, line := range lines {
			fields := strings.Split(line, "\t")
			if len(fields) < 2 {
				log.Printf("Not enough fields in line: %s\n", line)
				continue
			}
			pos = fields[0]
			if pos != position {
				//log.Printf("Position mismatch: %s != %s\n", pos, position)
				continue
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
			alleles[pos] = snp
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

func getAllIndividualAlleles(filepath string, individuals []string) map[string]map[string]uint8 {
	alleles := make(map[string]map[string]uint8)
	var pos, gt, idv string
	var err error
	var snp uint8
	var samples []string
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
		samples = strings.Split(fields[1], "\t")
		samples = samples[:len(samples)-1]
		for _, sample := range samples {
			fields = strings.Split(sample, "=")
			idv = fields[0]
			if _, ok := alleles[fields[0]]; !ok {
				alleles[fields[0]] = make(map[string]uint8)
			}
			gt = fields[1]
			gt, err = tools.NormalizeSnp(gt)
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
