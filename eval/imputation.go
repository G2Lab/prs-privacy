package main

import (
	"bytes"
	"encoding/json"
	"fmt"
	"log"
	"math"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
)

const (
	perAncestryAf = "ancestry"
	globalAf      = "global"
)

func prepareIBD() {
	guessedSnps := loadGuessedGenotypes(perAncestryAf)
	guessedIndividuals := make([]string, 0)
	for idv := range guessedSnps {
		guessedIndividuals = append(guessedIndividuals, idv)
	}
	sort.Strings(guessedIndividuals)
	preparePlinkInput(guessedIndividuals, guessedSnps)
	compressVcfFile("ibd/plink.vcf", "ibd/plink.vcf.gz")
	indexVcfFile("ibd/plink.vcf.gz")
	plinkGCTA("ibd/plink.vcf.gz")
	plinkKing("ibd/plink.vcf.gz")
}

func preparePlinkInput(guessedIndividuals []string, guessed map[string]map[string]map[string]uint8) {
	alleles := make(map[string]map[string][]string)
	positions := make(map[string][]int)
	var ref, alt, chr string
	var ok bool
	// Get all the guessed positions and corresponding alleles
	for chrId := 1; chrId <= 22; chrId++ {
		positionsMap := make(map[string]struct{})
		chr = strconv.Itoa(chrId)
		for idv := range guessed {
			for pos := range guessed[idv][chr] {
				positionsMap[pos] = struct{}{}
			}
		}
		positions[chr] = make([]int, 0)
		alleles[chr] = make(map[string][]string)
		for pos := range positionsMap {
			posInt, _ := strconv.Atoi(pos)
			positions[chr] = append(positions[chr], posInt)
			ref, alt = retrievePositionAlleles(chr, pos)
			alleles[chr][pos] = []string{ref, alt}
		}
		sort.Ints(positions[chr])
	}
	// Retrieve the true genotypes for all the guessed guessedIndividuals
	mainIndividuals, relativeIndividuals := make([]string, 0), make([]string, 0)
	for _, idv := range solver.All1000GenomesSamples() {
		if _, ok = guessed[idv]; ok {
			mainIndividuals = append(mainIndividuals, idv)
		}
	}
	for _, idv := range solver.AllRelativeSamples() {
		if _, ok = guessed[idv]; ok {
			relativeIndividuals = append(relativeIndividuals, idv)
		}
	}
	truthIndividuals := append(mainIndividuals, relativeIndividuals...)
	trueSNPs := make(map[string]map[string]map[string]string)
	var retrieved map[string]map[string]string
	for chrId := 1; chrId <= 22; chrId++ {
		chr = strconv.Itoa(chrId)
		strPos := make([]string, len(positions[chr]))
		for i, pos := range positions[chr] {
			strPos[i] = strconv.Itoa(pos)
		}
		retrieved = getIndividualPositionStrings(chr, strPos, mainIndividuals, tools.GG)
		transferSnps(trueSNPs, retrieved, chr)
		retrieved = getIndividualPositionStrings(chr, strPos, relativeIndividuals, tools.RL)
		transferSnps(trueSNPs, retrieved, chr)
		fmt.Printf("Chr %s: loaded true trueSNPs\n", chr)
	}
	references := getMajorGenotypesForGuessedLoci()

	formatHeader := "##fileformat=VCFv4.1\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
	lineTemplate := "%s\t%s\t.\t%s\t%s\t100\tPASS\t.\tGT"
	// Construct VCF
	samplesHeader := "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
	for _, idv := range guessedIndividuals {
		samplesHeader += fmt.Sprintf("\t%s", "$"+idv)
	}
	for anc := range references {
		samplesHeader += fmt.Sprintf("\t%s", "$"+anc)
	}
	for _, idv := range truthIndividuals {
		samplesHeader += fmt.Sprintf("\t%s", idv)
	}
	samplesHeader += "\n"
	vcfFile, err := os.OpenFile("ibd/plink.vcf", os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Cannot create plink ibd file: %v", err)
	}
	defer vcfFile.Close()
	_, err = vcfFile.WriteString(formatHeader)
	if err != nil {
		log.Fatalf("Cannot write format header to plink ibd file: %v", err)
	}
	_, err = vcfFile.WriteString(samplesHeader)
	if err != nil {
		log.Fatalf("Cannot write samples header to plink ibd file: %v", err)
	}
	var snp uint8
	for chrId := 1; chrId <= 22; chrId++ {
		chr = strconv.Itoa(chrId)
	positionLoop:
		for _, posId := range positions[chr] {
			pos := strconv.Itoa(posId)
			if _, ok = alleles[chr]; !ok || alleles[chr][pos][0] == "." || alleles[chr][pos][1] == "." {
				//if _, ok = alleles[chr]; !ok {
				continue
			}
			if _, ok = alleles[chr][pos]; !ok {
				alleles[chr][pos] = []string{".", "."}
			}
			line := fmt.Sprintf(lineTemplate, chr, pos, alleles[chr][pos][0], alleles[chr][pos][1])
			for _, idv := range guessedIndividuals {
				if snp, ok = guessed[idv][chr][pos]; ok {
					switch snp {
					case 0:
						line += "\t0|0"
					case 1:
						line += "\t1/0"
					case 2:
						line += "\t1|1"
					default:
						log.Fatalf("Invalid SNP value: %d", snp)
					}
				} else {
					continue positionLoop
					//line += "\t.|."
				}
			}
			for anc := range references {
				if snp, ok = references[anc][chr][pos]; ok {
					switch snp {
					case 0:
						line += "\t0|0"
					case 1:
						line += "\t1/0"
					case 2:
						line += "\t1|1"
					default:
						log.Fatalf("Invalid SNP value: %d", snp)
					}
				} else {
					continue positionLoop
					//line += "\t.|."
				}
			}
			for _, idv := range truthIndividuals {
				if _, ok = trueSNPs[idv][chr][pos]; ok {
					line += fmt.Sprintf("\t%s", trueSNPs[idv][chr][pos])
				} else {
					line += "\t0|0"
				}
			}
			_, err = vcfFile.WriteString(line + "\n")
			if err != nil {
				log.Fatalf("Cannot write to plink ibd file %s: %v", line, err)
			}
		}
	}
}

func compressVcfFile(pathIn, pathOut string) {
	prg := "bgzip"
	cmd := exec.Command(prg, "-c", pathIn)
	outputFile, err := os.OpenFile(pathOut, os.O_CREATE|os.O_WRONLY, 0644)
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

func indexVcfFile(path string) {
	runCommand("bcftools", []string{"index", path})
}

func transferSnps(to map[string]map[string]map[string]string, from map[string]map[string]string, chr string) {
	for idv := range from {
		if _, ok := to[idv]; !ok {
			to[idv] = make(map[string]map[string]string)
		}
		to[idv][chr] = from[idv]
	}
}

func runCommand(prg string, args []string) {
	cmd := exec.Command(prg, args...)
	var out bytes.Buffer
	var stderr bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &stderr
	err := cmd.Run()
	if err != nil {
		fmt.Printf("Output: %s\n Error: %s\n", out.String(), stderr.String())
		log.Fatalf("Error executing command: %v", err)
	}
}

func plinkGCTA(path string) {
	runCommand("plink", []string{"--vcf", path, "--make-rel", "square", "--out", "ibd/plink"})
}

func plinkKing(path string) {
	runCommand("plink", []string{"--vcf", path, "--make-king", "square", "--out", "ibd/plink"})
}

func retrievePositionAlleles(chr, pos string) (string, string) {
	ref, alt := retrievePositionAllelesDataset(chr, pos, tools.GG)
	if ref == "." || alt == "." {
		ref, alt = retrievePositionAllelesDataset(chr, pos, tools.RL)
	}
	return ref, alt
}

func retrievePositionAllelesDataset(chr, pos, dataset string) (string, string) {
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS\t%REF\t%ALT\n",
		"-r", fmt.Sprintf("%s:%s-%s", chr, pos, pos),
		tools.GetChromosomeFilepath(chr, dataset),
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

type IdvResult struct {
	Individual        string
	Ancestry          string
	SNPs              int
	GuessAccuracy     float32
	ReferenceAccuracy float32
}

type LocusResult struct {
	CorrectGuesses map[string]int
	TotalGuesses   map[string]int
	EAF            map[string]float32
	GWASEAF        float32
	SmallestPGS    int
	Density        float32
	EffectWeight   float32
}

func newLocusResult() *LocusResult {
	return &LocusResult{
		CorrectGuesses: make(map[string]int),
		TotalGuesses:   make(map[string]int),
		EAF:            make(map[string]float32),
	}
}

func guessAccuracy(afType string) {
	var guessed map[string]map[string]map[string]uint8
	var idvOutputFilename, lociOutputFilename string
	switch afType {
	case perAncestryAf:
		guessed = loadGuessedGenotypes(perAncestryAf)
		idvOutputFilename = "individuals.json"
		lociOutputFilename = "loci.json"
	case globalAf:
		guessed = loadGuessedGenotypes(globalAf)
		idvOutputFilename = "individuals_gaf.json"
		lociOutputFilename = "loci_gaf.json"
	default:
		log.Fatalf("Invalid allele frequency type: %s", afType)
	}
	relatives := solver.AllRelativeSamples()
	mainSamples := solver.All1000GenomesSamples()
	var err error
	var chr, idv, pos, filepath string
	var positionMap map[string]struct{}
	var positions []string
	trueSnps := make(map[string]map[string]map[string]uint8)
	for _, idv = range solver.All1000GenomesAndRelativeSamples() {
		trueSnps[idv] = make(map[string]map[string]uint8)
	}
	for chrId := 1; chrId <= 22; chrId++ {
		chr = strconv.Itoa(chrId)
		positionMap = make(map[string]struct{})
		for idv = range guessed {
			for pos = range guessed[idv][chr] {
				positionMap[pos] = struct{}{}
			}
		}
		positions = make([]string, 0, len(positionMap))
		for pos = range positionMap {
			positions = append(positions, pos)
		}
		retrieved := getIndividualPositionSNPs(chr, positions, mainSamples, tools.GG)
		for idv = range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
		retrieved = getIndividualPositionSNPs(chr, positions, relatives, tools.RL)
		for idv = range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
	}
	//
	references := getMajorGenotypesForGuessedLoci()
	ancestries := tools.LoadAncestry()
	idvResults := make([]*IdvResult, 0)
	locusResults := make(map[string]*LocusResult)
	var locus string
	for idv = range guessed {
		res := IdvResult{Individual: idv, Ancestry: pgs.GetIndividualAncestry(idv, ancestries),
			SNPs: 0, GuessAccuracy: 0, ReferenceAccuracy: 0}
		for chr = range guessed[idv] {
			for pos = range guessed[idv][chr] {
				// Individual accuracy
				if guessed[idv][chr][pos] == trueSnps[idv][chr][pos] {
					res.GuessAccuracy++
				}
				if guessed[idv][chr][pos] == references[res.Ancestry][chr][pos] {
					res.ReferenceAccuracy++
				}
				res.SNPs++
				// Locus accuracy
				locus = tools.MergeLocus(chr, pos)
				if _, ok := locusResults[locus]; !ok {
					locusResults[locus] = newLocusResult()
				}
				locusResults[locus].TotalGuesses[res.Ancestry]++
				if guessed[idv][chr][pos] == trueSnps[idv][chr][pos] {
					locusResults[locus].CorrectGuesses[res.Ancestry]++
				}
			}
		}
		res.GuessAccuracy /= float32(res.SNPs)
		res.ReferenceAccuracy /= float32(res.SNPs)
		idvResults = append(idvResults, &res)
	}
	resultFolder := "results/guessAccuracy"
	filepath = path.Join(resultFolder, idvOutputFilename)
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(idvResults); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Saved accuracy results")
	//
	pgsFile, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer pgsFile.Close()
	decoder := json.NewDecoder(pgsFile)
	var pgsToNumVariants map[string]int
	err = decoder.Decode(&pgsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}
	allPgs := make([]string, 0)
	for pgsID := range pgsToNumVariants {
		allPgs = append(allPgs, pgsID)
	}
	sort.Slice(allPgs, func(i, j int) bool {
		return pgsToNumVariants[allPgs[i]] < pgsToNumVariants[allPgs[j]]
	})
	observedLoci := make(map[string]struct{})
	var ok bool
	for _, pgsID := range allPgs {
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		err = p.LoadStats(tools.GG)
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return
		}
		maxw := p.FindMaxAbsoluteWeight()
		minw := p.FindMinAbsoluteWeight()
		for i, locus := range p.Loci {
			if _, ok = observedLoci[locus]; ok {
				continue
			}
			observedLoci[locus] = struct{}{}
			for ppl := range p.PopulationStats {
				locusResults[locus].EAF[ppl] = p.PopulationStats[ppl].AF[i][p.EffectAlleles[i]]
			}
			locusResults[locus].GWASEAF = float32(p.StudyEAF[i])
			locusResults[locus].SmallestPGS = pgsToNumVariants[pgsID]
			locusResults[locus].Density = float32(p.NumVariants) /
				float32(tools.Log3(p.FindMaxAbsoluteWeight()*math.Pow(10, float64(p.WeightPrecision))))
			weight, err := p.Weights[i].Float64()
			if err != nil {
				log.Printf("Error converting weight %s to float64: %v\n", p.Weights[i].String(), err)
				continue
			}
			locusResults[locus].EffectWeight = float32((math.Abs(weight) - minw) / (maxw - minw))
		}
	}
	eafFile, err := os.OpenFile(path.Join(resultFolder, lociOutputFilename), os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Println("Error opening eaf file:", err)
		return
	}
	defer eafFile.Close()
	encoder = json.NewEncoder(eafFile)
	if err = encoder.Encode(locusResults); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Completed")
}

// Find the most frequent genotype for each position
func getMajorGenotypesForGuessedLoci() map[string]map[string]map[string]uint8 {
	// Find the most frequent genotype for each position
	gtpFrequencies := loadGenotypeFrequenciesForGuessed()
	references := make(map[string]map[string]map[string]uint8)
	var gtp uint8
	for anc := range gtpFrequencies {
		references[anc] = make(map[string]map[string]uint8)
		for chrStr := range gtpFrequencies[anc] {
			references[anc][chrStr] = make(map[string]uint8)
			for pos := range gtpFrequencies[anc][chrStr] {
				gtp = 0
				for j := 1; j <= 2; j++ {
					if gtpFrequencies[anc][chrStr][pos][j] > gtpFrequencies[anc][chrStr][pos][gtp] {
						gtp = uint8(j)
					}
				}
				references[anc][chrStr][pos] = gtp
			}
		}
	}
	return references
}

func loadGuessedGenotypes(afType string) map[string]map[string]map[string]uint8 {
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
	var prefix string
	switch afType {
	case perAncestryAf:
		prefix = "guesses"
	case globalAf:
		prefix = "af-guesses"
	default:
		log.Fatalf("Unknown AF type: %s", afType)
	}
	for _, object := range dir {
		if !object.IsDir() && strings.HasPrefix(object.Name(), prefix) && strings.HasSuffix(object.Name(), ".json") {
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

func calculateGenotypeFrequenciesForGuessed() {
	guessed := loadGuessedGenotypes(perAncestryAf)
	missingGenotypeFrequencies := []float32{(1 - pgs.MissingEAF) * (1 - pgs.MissingEAF),
		2 * pgs.MissingEAF * (1 - pgs.MissingEAF), pgs.MissingEAF * pgs.MissingEAF}
	folder := "data/frequencies"
	samples := solver.All1000GenomesSamples()
	var filepath string
	var outFile *os.File
	var err error
	var ok bool
	var positions []string
	var chrStr, pos, idv, anc string
	var positionMap map[string]struct{}
	frequencies := make(map[string]map[string]map[string][]float32)
	ancestries := tools.LoadAncestry()
	for _, ppl := range pgs.ANCESTRIES {
		frequencies[ppl] = make(map[string]map[string][]float32)
		for chr := 1; chr <= 22; chr++ {
			frequencies[ppl][strconv.Itoa(chr)] = make(map[string][]float32)
		}
	}
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
		retrieved := getIndividualPositionSNPs(chrStr, positions, samples, tools.GG)
		for idv = range retrieved {
			anc = pgs.GetIndividualAncestry(idv, ancestries)
			for pos = range retrieved[idv] {
				if _, ok := frequencies[anc][chrStr][pos]; !ok {
					frequencies[anc][chrStr][pos] = make([]float32, 3)
				}
				frequencies[anc][chrStr][pos][retrieved[idv][pos]]++
			}
		}
		for anc = range frequencies {
			for pos = range frequencies[anc][chrStr] {
				for i := range frequencies[anc][chrStr][pos] {
					frequencies[anc][chrStr][pos][i] /= float32(len(retrieved))
				}
			}
		}

		for _, ppl := range pgs.ANCESTRIES {
			for _, pos := range positions {
				if _, ok = frequencies[ppl][chrStr][pos]; !ok {
					fmt.Printf("%s: position %s:%s is missing\n", ppl, chrStr, pos)
					frequencies[ppl][chrStr][pos] = missingGenotypeFrequencies
				}
			}
		}
	}
	filepath = path.Join(folder, "guessed.json")
	outFile, err = os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening output file: %v", err)
	}
	encoder := json.NewEncoder(outFile)
	if err = encoder.Encode(frequencies); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	outFile.Close()
	fmt.Printf("=== Completed ===\n")
}

func loadGenotypeFrequenciesForGuessed() map[string]map[string]map[string][]float32 {
	folder := "data/frequencies"
	filepath := path.Join(folder, "guessed.json")
	file, err := os.Open(filepath)
	if err != nil {
		log.Fatalf("Cannot open file %s: %v", filepath, err)
	}
	defer file.Close()
	var frequencies map[string]map[string]map[string][]float32
	decoder := json.NewDecoder(file)
	if err = decoder.Decode(&frequencies); err != nil {
		log.Fatalf("Cannot decode json file %s: %v", filepath, err)
	}
	return frequencies
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

func getIndividualPositionStrings(chr string, positions []string, individuals []string, dataset string) map[string]map[string]string {
	prg := "bcftools"
	args := []string{
		"query",
		"-s", strings.Join(individuals, ","),
		"-f", "%POS-[%SAMPLE=%GT\t]\n",
		tools.GetChromosomeFilepath(chr, dataset),
		"-r", "",
	}
	var pos, idv string
	var ok bool
	snps := make(map[string]map[string]string)
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
				if _, ok = snps[idv]; !ok {
					snps[idv] = make(map[string]string)
				}
				snps[idv][pos] = fields[1]
			}
		}
	}
	return snps
}
