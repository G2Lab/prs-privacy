package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
	"gonum.org/v1/gonum/stat/distuv"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"path/filepath"
	"regexp"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
var memprofile = flag.String("memprofile", "", "write cpu profile to file")

const (
	DeterminismLimit       = 34
	ScratchSolvingSnpLimit = 44
	ConfidenceThreshold    = 5
)

func main() {
	expr := flag.String("e", "", "Experiment type")
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	if *memprofile != "" {
		f, err := os.Create(*memprofile)
		if err != nil {
			log.Fatal(err)
		}
		defer f.Close()
		runtime.MemProfileRate = 1
		defer pprof.WriteHeapProfile(f)
	}
	switch *expr {
	case "sequence":
		sequenceSolving()
	case "sequenceaf":
		sequenceSolvingAF()
	case "findall":
		findAllSolutions()
	case "impute":
		imputeWorkflow()
	case "accuracy":
		guessAccuracy(perAncestryAf)
	case "gaccuracy":
		guessAccuracy(globalAf)
	case "genfreq":
		calculateGenotypeFrequenciesOnlyGuessed()
	case "predict":
		predictPRS()
	case "king":
		linkingWithKing()
	case "ibd":
		prepareIBD()
	case "findprs":
		findUnsolvablePRSWithOverlap()
	case "uniqueness_gg":
		uniquenessExperiment(tools.GG)
	case "uniqueness_uk":
		uniquenessExperiment(tools.UKB)
	case "scores":
		calculateScoresAndStats()
	}
	//kinshipExperiment()
	//kingTest()
	//calculateGenotypeFrequencies()
	//imputeWorkflow()
	//linkingWithImputation()
	//linkingWithGuessed()
	//evaluateImputation()
}

func predictPRS() {
	type Result struct {
		Known      int
		Total      int
		Max        string
		Min        string
		Means      map[string]float64
		Stds       map[string]float64
		Ancestries []string
		Predicted  []string
		Real       []string
	}
	results := make(map[string]*Result)
	guessed := loadGuessedGenotypes(perAncestryAf)
	individuals := solver.All1000GenomesAndRelativeSamples()
	populations := tools.LoadAncestry()

	pgsIDs := findUnsolvablePRSWithOverlap()
	for _, pgsID := range pgsIDs {
		p := pgs.NewPGS()
		err := p.LoadCatalogFile(path.Join("catalog", pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		p.LoadStats(tools.GG)

		minScore, maxScore := p.FindMinAndMaxScores()
		means, stds := p.EstimateMeanAndStd()
		results[p.PgsID] = &Result{
			Known:      0,
			Total:      len(p.Loci),
			Max:        maxScore.String(),
			Min:        minScore.String(),
			Means:      means,
			Stds:       stds,
			Ancestries: make([]string, 0),
			Predicted:  make([]string, 0),
			Real:       make([]string, 0),
		}
		cohort := solver.NewCohort(p, tools.GG)
		var majorFreq float32
		var guess, majorGtp uint8
		var anc string
		var ok bool
		divisor := new(apd.Decimal).SetInt64(int64(len(p.Weights) * pgs.Ploidy))
		for k, idv := range individuals {
			predictedScore := apd.New(0, 0)
			anc = pgs.GetIndividualAncestry(idv, populations)
			for i, locus := range p.Loci {
				chr, pos := tools.SplitLocus(locus)
				guess, ok = guessed[idv][chr][pos]
				if !ok {
					// we do not have a guess for the locus, so we assume the most common genotype
					majorGtp = 0
					majorFreq = p.PopulationStats[anc].GF[i][majorGtp]
					for j, freq := range p.PopulationStats[anc].GF[i] {
						if freq > majorFreq {
							majorGtp = uint8(j)
						}
					}
					guess = majorGtp
				} else if k == 0 {
					results[p.PgsID].Known++
				}
				switch {
				case (guess == 0 && p.EffectAlleles[i] == 0) || (guess == 2 && p.EffectAlleles[i] == 1):
					p.Context.Add(predictedScore, predictedScore, p.Weights[i])
					p.Context.Add(predictedScore, predictedScore, p.Weights[i])
				case guess == 1:
					p.Context.Add(predictedScore, predictedScore, p.Weights[i])
				default:
					if guess != 0 && guess != 2 {
						log.Printf("Unknown guess genotype at %s: %d\n", locus, guess)
					}
				}
			}
			_, err := p.Context.Quo(predictedScore, predictedScore, divisor)
			if err != nil {
				log.Println("Error normalizing the score:", err)
				return
			}
			results[p.PgsID].Ancestries = append(results[p.PgsID].Ancestries, anc)
			results[p.PgsID].Predicted = append(results[p.PgsID].Predicted, predictedScore.String())
			results[p.PgsID].Real = append(results[p.PgsID].Real, cohort[idv].Score.String())
		}
	}
	filePath := path.Join("results/predict", "prediction.json")
	resFile, err := os.OpenFile(filePath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func percentile(num, lower, upper float64) float64 {
	if num < lower {
		return 0
	}
	if num > upper {
		return 100
	}
	return (num - lower) / (upper - lower)
}

func findUnsolvablePRSWithOverlap() []string {
	guessedPerIndividual := loadGuessedGenotypes(perAncestryAf)
	guessedLoci := make(map[string]struct{})
	for idv := range guessedPerIndividual {
		for chr := range guessedPerIndividual[idv] {
			for pos := range guessedPerIndividual[idv][chr] {
				guessedLoci[fmt.Sprintf("%s:%s", chr, pos)] = struct{}{}
			}
		}
	}
	ids, err := fewerVariantsPGS(50, 200)
	if err != nil {
		log.Println("Error:", err)
		return nil
	}
	catalogFolder := "catalog"
	unknownPgs := make(map[string][]int)
idLoop:
	for _, id := range ids {
		//fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalogFolder, id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error:", err)
			continue
		}
		for _, locus := range p.Loci {
			if strings.HasPrefix(locus, "X:") || strings.HasPrefix(locus, "Y:") || strings.HasPrefix(locus, "XY:") {
				continue idLoop
			}
			if !isNumber(strings.Split(locus, ":")[1]) || !isNumberInRange(strings.Split(locus, ":")[0], 1, 22) {
				continue idLoop
			}
		}
		unknownPgs[id] = make([]int, 2)
		for _, locus := range p.Loci {
			if _, ok := guessedLoci[locus]; ok {
				unknownPgs[id][0]++
			}
			unknownPgs[id][1]++
		}
	}
	results := make([]string, 0)
	for id, counts := range unknownPgs {
		if counts[0] > 0 && float64(counts[0])/float64(counts[1]) > 0.1 {
			p := pgs.NewPGS()
			fmt.Printf("==== %s ====\n", id)
			p.LoadCatalogFile(filepath.Join(catalogFolder, id+"_hmPOS_GRCh37.txt"))
			err = p.LoadStats(tools.GG)
			if err != nil {
				continue
			}
			results = append(results, id)
		}
	}
	sort.Slice(results, func(i, j int) bool {
		return float64(unknownPgs[results[i]][0])/float64(unknownPgs[results[i]][1]) >
			float64(unknownPgs[results[j]][0])/float64(unknownPgs[results[j]][1])
	})
	for _, id := range results {
		fmt.Printf("%s: %d/%d\n", id, unknownPgs[id][0], unknownPgs[id][1])
	}
	fmt.Printf("%d PRSs found\n", len(results))
	return results
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

func kingTest() {
	fmt.Println("King test")
	related := solver.ReadRelatedIndividuals()
	individuals := solver.All1000GenomesAndRelativeSamples()
	db := make(map[string]map[string]uint8)
	for _, idv := range individuals {
		db[idv] = make(map[string]uint8)
	}

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
	file, err = os.Open("results/validated_loci.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder = json.NewDecoder(file)
	var lociToPgs map[string][]string
	err = decoder.Decode(&lociToPgs)
	if err != nil {
		log.Println("Error decoding validated loci:", err)
		return
	}

	allPgs := make([]string, 0, len(idsToNumVariants))
	for id, _ := range idsToNumVariants {
		allPgs = append(allPgs, id)
	}
	pgsToLoci := make(map[string]map[string]struct{})
	for locus, pgsIDs := range lociToPgs {
		for _, pgsID := range pgsIDs {
			if _, ok := pgsToLoci[pgsID]; !ok {
				pgsToLoci[pgsID] = make(map[string]struct{})
			}
			pgsToLoci[pgsID][locus] = struct{}{}
		}
	}

	thresholds := []int{1000, 2000, 3000}
	thresholdCounter := 0
	type Result struct {
		TruePhi  float32
		HighPhi  float32
		Position int
	}
	results := make(map[int][]*Result)
	type relation struct {
		target string
		phi    float32
	}

	var pgsID string
	firstIdv := individuals[0]
	guessThreshold := 3500
	targetSamples := solver.All1000GenomesSamples()
	for {
		if len(allPgs) == 0 {
			break
		}
		sort.Slice(allPgs, func(i, j int) bool {
			return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
		})
		pgsID = allPgs[0]
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		err = p.LoadStats(tools.GG)
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return
		}
		cohort := solver.NewCohort(p, tools.GG)
		for _, locus := range p.Loci {
			for idv := range db {
				if _, ok := db[idv][locus]; ok {
					continue
				}
				db[idv][locus] = cohort[idv].Genotype[0] + cohort[idv].Genotype[1]
			}
		}
		if len(db[firstIdv]) > thresholds[thresholdCounter] {
			fmt.Printf("=== SNPs covered %d\n", len(db[firstIdv]))
			results[thresholds[thresholdCounter]] = make([]*Result, 0)
			for idv, relatives := range related {
				relations := make([]relation, 0)
				for _, other := range targetSamples {
					relations = append(relations, relation{other, kingRobust(db[idv], db[other])})
				}
				sort.Slice(relations, func(i, j int) bool {
					return relations[i].phi > relations[j].phi
				})

				for _, relative := range relatives {
					pos := 0
					for i, rel := range relations {
						if rel.target == relative {
							pos = i
							break
						}
					}
					results[thresholds[thresholdCounter]] = append(results[thresholds[thresholdCounter]], &Result{
						TruePhi:  kingRobust(db[idv], db[relative]),
						HighPhi:  relations[0].phi,
						Position: pos,
					})
				}
			}
			thresholdCounter++
		}
		if len(db[firstIdv]) > guessThreshold || thresholdCounter > len(thresholds)-1 {
			break
		}
		allPgs = allPgs[1:]
		for _, locus := range p.Loci {
			for _, id := range lociToPgs[locus] {
				if id == pgsID {
					continue
				}
				if _, ok := pgsToLoci[id][locus]; ok {
					delete(pgsToLoci[id], locus)
				}
			}
		}
	}
	resultFolder := "results/kinship"
	filepath := path.Join(resultFolder, "truth.json")
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

	//var highestPhi float32
	//var highestPhiIdv string
	//var kr float32
	//for idv, relatives := range related {
	//	fmt.Printf("### %s\n", idv)
	//	highestPhi = 0
	//	highestPhiIdv = ""
	//	for other := range db {
	//		if other == idv {
	//			continue
	//		}
	//		//kh = kingHomo(db[idv], db[other], alfreq[idv])
	//		kr = kingRobust(db[idv], db[other])
	//		if kr > highestPhi {
	//			highestPhi = kr
	//			highestPhiIdv = other
	//		}
	//	}
	//	fmt.Printf("--- %s->%s: %.3f\n", idv, highestPhiIdv, highestPhi)
	//	for _, relative := range relatives {
	//		fmt.Printf("++++ %s->%s: %.3f, %.3f\n", idv, relative, kingHomo(db[idv], db[relative], alfreq[idv]),
	//			kingRobust(db[idv], db[relative]))
	//	}
	//}

}

func kinshipExperiment() {
	fmt.Printf("Kinship experiment\n")
	related := solver.ReadRelatedIndividuals()
	relatives := make([]string, 0)
	for ind, _ := range related {
		relatives = append(relatives, ind)
	}
	sort.Strings(relatives)
	//fmt.Println(relatives)
	//relative := relatives[0]
	//relative := "HG03948"

	individuals := solver.All1000GenomesSamples()
	db := make(map[string]map[string]uint8)
	for _, sample := range individuals {
		db[sample] = make(map[string]uint8)
	}

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
	file, err = os.Open("results/validated_loci.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder = json.NewDecoder(file)
	var lociToPgs map[string][]string
	err = decoder.Decode(&lociToPgs)
	if err != nil {
		log.Println("Error decoding validated loci:", err)
		return
	}

	type Result struct {
		TruePhi  float32
		HighPhi  float32
		Position int
		Accuracy float32
	}
	results := make([]*Result, 0)
	type relation struct {
		target string
		phi    float32
	}

	var pgsID string
	populations := tools.LoadAncestry()
	guessThreshold := 2500

	chunkNum, chunkSize := getChunkInfo(len(relatives))
	fmt.Println(relatives[chunkNum*chunkSize : (chunkNum+1)*chunkSize])
	var relative string
	for c := chunkNum * chunkSize; c < (chunkNum+1)*chunkSize; c++ {
		if c >= len(relatives) {
			break
		}
		relative = relatives[c]
		fmt.Printf("--------- %s --------\n", relative)

		allPgs := make([]string, 0, len(idsToNumVariants))
		for id, _ := range idsToNumVariants {
			allPgs = append(allPgs, id)
		}
		pgsToLoci := make(map[string]map[string]struct{})
		for locus, pgsIDs := range lociToPgs {
			for _, pgsID := range pgsIDs {
				if _, ok := pgsToLoci[pgsID]; !ok {
					pgsToLoci[pgsID] = make(map[string]struct{})
				}
				pgsToLoci[pgsID][locus] = struct{}{}
			}
		}

		indPop := populations[relative]
		if strings.Contains(indPop, ",") {
			indPop = strings.Split(indPop, ",")[0]
		}
		var solutions [][]uint8
		guessedSnps := make(map[string]uint8)
		trueSnps := make(map[string]uint8)
		for {
			if len(allPgs) == 0 {
				break
			}
			sort.Slice(allPgs, func(i, j int) bool {
				if len(pgsToLoci[allPgs[i]]) == len(pgsToLoci[allPgs[j]]) {
					return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
				}
				return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
			})
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats(tools.GG)
			if err != nil {
				log.Printf("Error loading stats: %v\n", err)
				return
			}
			recoveredSnps := make(map[int]uint8)
			for l, locus := range p.Loci {
				if guess, ok := guessedSnps[locus]; ok {
					recoveredSnps[l] = guess
				}
			}
			fmt.Printf("Total SNPs %d, unknown %d\n", len(p.Loci), len(p.Loci)-len(recoveredSnps))
			fmt.Printf("Recovered snps: %v\n", recoveredSnps)
			cohort := solver.NewCohort(p, tools.GG)
			for _, locus := range p.Loci {
				for idv := range db {
					if _, ok := db[idv][locus]; ok {
						continue
					}
					db[idv][locus] = cohort[idv].Genotype[0] + cohort[idv].Genotype[1]
				}
			}
			solutions = findSolutions(p, cohort, relative, indPop, recoveredSnps)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				fmt.Println("^^^ No solutions with extra loci, calculating from scratch ^^^")
				solutions = findSolutions(p, cohort, relative, indPop, make(map[int]uint8))
			}
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[relative].Genotype))
				for i, locus := range p.Loci {
					if _, ok := guessedSnps[locus]; ok {
						continue
					}
					guessedSnps[locus] = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
				}
			} else {
				fmt.Println("!!!!!!! Still could not find a solution !!!!!!!")
			}
			for i, locus := range p.Loci {
				if _, ok := trueSnps[locus]; ok {
					continue
				}
				trueSnps[locus] = cohort[relative].Genotype[pgs.Ploidy*i] + cohort[relative].Genotype[pgs.Ploidy*i+1]
			}
			fmt.Printf("Guessed %d\n", len(guessedSnps))
			if len(guessedSnps) > guessThreshold {
				var accuracy float32 = 0.0
				for locus, snp := range trueSnps {
					if guessed, ok := guessedSnps[locus]; ok && snp == guessed {
						accuracy += 1
					}
				}
				fmt.Printf("### CountAccuracy: %.3f\n", accuracy/float32(len(trueSnps)))

				relations := make([]relation, 0)
				for _, other := range individuals {
					relations = append(relations, relation{other, kingRobust(guessedSnps, db[other])})
				}
				sort.Slice(relations, func(i, j int) bool {
					return relations[i].phi > relations[j].phi
				})

				pos := 0
				matches := related[relative]
				var matched string
				for i, rel := range relations {
					for _, match := range matches {
						if rel.target == match {
							pos = i
							matched = match
							break
						}
					}
				}
				results = append(results, &Result{
					TruePhi:  kingRobust(guessedSnps, db[matched]),
					HighPhi:  relations[0].phi,
					Position: pos,
					Accuracy: accuracy / float32(len(trueSnps)),
				})

				break
			}
			allPgs = allPgs[1:]
			for _, locus := range p.Loci {
				for _, id := range lociToPgs[locus] {
					if id == pgsID {
						continue
					}
					if _, ok := pgsToLoci[id][locus]; ok {
						delete(pgsToLoci[id], locus)
					}
				}
			}
		}
	}
	resultFolder := "results/kinship"
	filepath := path.Join(resultFolder, fmt.Sprintf("chunk%d.json", chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func calculateScoresAndStats() {
	file, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	var pgsToNumVariants map[string]int
	err = decoder.Decode(&pgsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}
	allPgs := make([]string, 0, len(pgsToNumVariants))
	for id, _ := range pgsToNumVariants {
		allPgs = append(allPgs, id)
	}
	sort.Slice(allPgs, func(i, j int) bool {
		if pgsToNumVariants[allPgs[i]] == pgsToNumVariants[allPgs[j]] {
			return allPgs[i] < allPgs[j]
		}
		return pgsToNumVariants[allPgs[i]] < pgsToNumVariants[allPgs[j]]
	})
	pgsNum, err := strconv.Atoi(os.Args[2])
	if err != nil {
		log.Fatalf("Error parsing pgsNum %s: %v", os.Args[2], err)
	}
	if pgsNum < 0 || pgsNum >= len(allPgs) {
		log.Printf("Too high pgsNum %d\n", pgsNum)
		return
	}
	fmt.Printf("======== %s ========\n", allPgs[pgsNum])
	p := pgs.NewPGS()
	err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, allPgs[pgsNum]+"_hmPOS_GRCh37.txt"))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats(tools.UKB)
	solver.NewCohort(p, tools.UKB)
	fmt.Printf("%s complete\n", allPgs[pgsNum])
}

func uniquenessExperiment(dataset string) {
	fmt.Printf("--- Uniqueness experiment ---\n")
	type Result struct {
		PgsID                       string
		NumVariants                 int
		TotalPresentScores          int
		TotalPossibleScores         int
		RealPercentageUnique        float64
		PredictedPercentageUnique   float64
		RealMedianAnonymitySet      int
		PredictedMedianAnonymitySet int
	}

	file, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	var pgsToNumVariants map[string]int
	err = decoder.Decode(&pgsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids:", err)
		return
	}
	allPgs := make([]string, 0, len(pgsToNumVariants))
	for id, _ := range pgsToNumVariants {
		allPgs = append(allPgs, id)
	}
	sort.Slice(allPgs, func(i, j int) bool {
		return pgsToNumVariants[allPgs[i]] < pgsToNumVariants[allPgs[j]]
	})

	numWorkers := 3
	tasks := make(chan string, 1)
	results := make([]*Result, 0)
	var mutex sync.Mutex
	var wg sync.WaitGroup
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for pgsID := range tasks {
				var sc string
				var sf float64
				var realNumUnique int
				p := pgs.NewPGS()
				err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
				if err != nil {
					log.Printf("Error loading catalog file: %v\n", err)
					return
				}
				p.LoadStats(dataset)
				fmt.Printf("======== %s (%d) ========\n", p.PgsID, len(p.Loci))
				cohort := solver.NewCohort(p, dataset)
				fmt.Printf("%s: cohort size %d\n", p.PgsID, len(cohort))
				scores := make(map[string]int)
				scoresAsStrings := make([]string, 0)
				scoresAsFloats := make([]float64, 0)
				for idv := range cohort {
					sc = cohort[idv].Score.String()
					scores[sc]++
					scoresAsStrings = append(scoresAsStrings, sc)
					sf, err = cohort[idv].Score.Float64()
					if err != nil {
						log.Println("Error converting score to float:", err)
						return
					}
					scoresAsFloats = append(scoresAsFloats, sf)
				}
				realNumUnique = 0
				anonsets := make([]int, 0)
				for _, count := range scores {
					if count == 1 {
						realNumUnique++
					}
					anonsets = append(anonsets, count)
				}
				realMedianAnonSize := median(anonsets)
				minScore, maxScore := p.FindMinAndMaxScores()
				minScoreFloat, _ := minScore.Float64()
				maxScoreFloat, _ := maxScore.Float64()
				numPossibleSubsetSums := calculateNumPossibleSums(p.Weights)
				predictedMeans, predictedStds := p.EstimateMeanAndStd()
				predictedNumUniquePredictedStats, predictedMedianAnonSize, numPossibleScores :=
					estimateNumUniqueScores(len(cohort), numPossibleSubsetSums, predictedMeans["ALL"],
						predictedStds["ALL"], minScoreFloat, maxScoreFloat, p.WeightPrecision, 13)
				if numPossibleScores == -1 || math.IsNaN(predictedNumUniquePredictedStats) {
					fmt.Printf("Too many possible scores for %s\n", pgsID)
					continue
				}
				result := &Result{
					PgsID:                       pgsID,
					NumVariants:                 pgsToNumVariants[pgsID],
					TotalPresentScores:          len(scores),
					TotalPossibleScores:         numPossibleScores,
					RealPercentageUnique:        float64(realNumUnique) * 100 / float64(len(cohort)),
					PredictedPercentageUnique:   predictedNumUniquePredictedStats * 100 / float64(len(cohort)),
					RealMedianAnonymitySet:      realMedianAnonSize,
					PredictedMedianAnonymitySet: predictedMedianAnonSize,
				}
				mutex.Lock()
				results = append(results, result)
				mutex.Unlock()
				fmt.Printf("%s -- Real ratio: %.2f%%, predicted ratio: %.2f%%\nReal median set size: %d, "+
					"predicted set size: %d\n", pgsID, result.RealPercentageUnique, result.PredictedPercentageUnique,
					result.RealMedianAnonymitySet, result.PredictedMedianAnonymitySet)
			}
		}()
	}
	for _, pgsID := range allPgs {
		tasks <- pgsID
	}
	fmt.Println("---- All tasks sent ----")
	close(tasks)
	wg.Wait()
	for _, res := range results {
		fmt.Println(*res)
	}
	resultFolder := "results/uniqueness"
	filepath := path.Join(resultFolder, fmt.Sprintf("scores_%s.json", dataset))
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

func estimateNumUniqueScores(M int, numSubsetSums int, mu, sigma, minVal, maxVal float64, precision uint32,
	powerBound int) (float64, int, int) {
	normalDist := distuv.Normal{Mu: mu, Sigma: sigma}
	step := math.Pow10(-int(precision))
	numPossibleScoresDuePrecision := int((maxVal-minVal)/step) + 1
	numPossibleScores := numPossibleScoresDuePrecision
	if numSubsetSums > 0 && numSubsetSums < numPossibleScoresDuePrecision {
		fmt.Printf("Number of subsums: %d < In-between scores: %d\n", numSubsetSums, numPossibleScoresDuePrecision)
		step = (maxVal - minVal) / float64(numSubsetSums)
		numPossibleScores = numSubsetSums
	} else {
		fmt.Printf("Num subsums: %d > In-between scores: %d\n", numSubsetSums, numPossibleScoresDuePrecision)
	}
	if numPossibleScores > int(math.Pow10(powerBound)) {
		return -1, -1, -1
	}

	type result struct {
		unique float64
		sets   map[uint16]uint16
	}
	numWorkers := 8
	segmentSize := (maxVal - minVal) / float64(numWorkers)
	results := make(chan result, numWorkers)
	var wg sync.WaitGroup
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			var lowerBound, upperBound, prob float64
			var averageCount float64
			var roundedCount uint16
			anonSets := make(map[uint16]uint16)
			segmentTotalUnique := 0.0
			for x := minVal + float64(i)*segmentSize; x < minVal+float64(i+1)*segmentSize; x += step {
				lowerBound = x - step/2
				upperBound = x + step/2
				if x == minVal {
					lowerBound = x
				}
				if x == maxVal {
					upperBound = x
				}
				prob = normalDist.CDF(upperBound) - normalDist.CDF(lowerBound)
				averageCount = prob * float64(M)
				segmentTotalUnique += averageCount * math.Pow(1-prob, float64(M-1))
				roundedCount = uint16(math.Round(averageCount))
				if averageCount < 1 {
					randomValue := rand.Float64()
					if randomValue >= averageCount {
						continue
					}
					roundedCount = 1
				}
				if math.IsNaN(averageCount) {
					continue
				}
				if _, exists := anonSets[roundedCount]; exists {
					anonSets[roundedCount]++
				} else {
					anonSets[roundedCount] = 1
				}
			}
			results <- result{segmentTotalUnique, anonSets}
		}(i)
	}
	wg.Wait()
	close(results)
	totalUnique := 0.0
	allSets := make(map[uint16]uint16)
	for res := range results {
		totalUnique += res.unique
		for sz, count := range res.sets {
			if _, exists := allSets[sz]; exists {
				allSets[sz] += count
			} else {
				allSets[sz] = count
			}
		}
	}
	allSizes := make([]uint16, 0, len(allSets))
	var totalCount uint16 = 0
	for sz, count := range allSets {
		allSizes = append(allSizes, sz)
		totalCount += count
	}
	// Sorting the anonymity-set sizes
	sort.Slice(allSizes, func(i, j int) bool {
		return allSizes[i] < allSizes[j]
	})
	cumulativeCount := uint16(0)
	halfTotal := totalCount / 2
	for _, sz := range allSizes {
		cumulativeCount += allSets[sz]
		if cumulativeCount >= halfTotal {
			return totalUnique, int(sz), numPossibleScores
		}
	}
	if len(allSizes) == 0 {
		return totalUnique, 0, numPossibleScores
	}
	return totalUnique, int(allSizes[len(allSizes)-1]), numPossibleScores
}

//func estimateMedianAnonymitySetSize(M int, numSubsetSums int, mu, sigma, minVal, maxVal float64, precision uint32) float64 {
//	normalDist := distuv.Normal{
//		Mu:    mu,
//		Sigma: sigma,
//	}
//	step := math.Pow10(-int(precision))
//	numPossibleScoresDuePrecision := int((maxVal-minVal)/step) + 1
//	if numSubsetSums > 0 && numSubsetSums < numPossibleScoresDuePrecision {
//		step = (maxVal - minVal) / float64(numSubsetSums)
//	}
//	discreteMean := math.Round(mu/step) * step
//	pMean := normalDist.CDF(discreteMean+step/2) - normalDist.CDF(discreteMean-step/2)
//	return float64(M) * pMean
//}

func median(values []int) int {
	sort.Ints(values)
	n := len(values)
	if n%2 == 0 {
		return (values[n/2-1] + values[n/2]) / 2
	}
	return values[n/2]
}

func calculateNumPossibleSums(weights []*apd.Decimal) int {
	weightRepetitions := make(map[string]int)
	for i := range weights {
		w := weights[i].String()
		weightRepetitions[w]++
	}
	numPossibleSums := 1
	for _, count := range weightRepetitions {
		numPossibleSums *= 2*count + 1
	}
	return numPossibleSums
}

func CalculateMeanAndStd(scores []float64) (float64, float64) {
	if len(scores) == 0 {
		return 0, 0
	}
	sum := 0.0
	for _, score := range scores {
		sum += score
	}
	mean := sum / float64(len(scores))

	varianceSum := 0.0
	for _, score := range scores {
		varianceSum += math.Pow(score-mean, 2)
	}
	std := math.Sqrt(varianceSum / float64(len(scores)))

	return mean, std
}

func sequenceSolving() {
	thresholds := []int{500, 1000, 1500, 2000, 2500}
	fmt.Printf("Sequential solving for %d\n", thresholds)
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
	file, err = os.Open("results/validated_loci.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder = json.NewDecoder(file)
	var lociToPgs map[string][]string
	err = decoder.Decode(&lociToPgs)
	if err != nil {
		log.Println("Error decoding validated loci:", err)
		return
	}

	individuals := solver.All1000GenomesAndRelativeSamples()
	sort.Strings(individuals)

	chunkNum, chunkSize := getChunkInfo(len(individuals))
	fmt.Printf("Chunk %d, size %d\n", chunkNum, chunkSize)
	end := (chunkNum + 1) * chunkSize
	if end > len(individuals) {
		end = len(individuals)
	}
	fmt.Println(individuals[chunkNum*chunkSize : end])

	type Result struct {
		Individual string
		Ancestry   string
		Accuracies map[int]float32
	}
	type Guess struct {
		Individual string
		Ancestry   string
		SNPs       map[string]uint8
	}
	results := make([]*Result, 0)
	guesses := make([]*Guess, 0)
	populations := tools.LoadAncestry()
	var pgsID string
	var individual string
	for c := chunkNum * chunkSize; c < (chunkNum+1)*chunkSize; c++ {
		if c >= len(individuals) {
			break
		}
		individual = individuals[c]
		fmt.Printf("--------- %s --------\n", individual)

		allPgs := make([]string, 0, len(idsToNumVariants))
		for id, _ := range idsToNumVariants {
			allPgs = append(allPgs, id)
		}
		pgsToLoci := make(map[string]map[string]struct{})
		for locus, pgsIDs := range lociToPgs {
			for _, pgsID := range pgsIDs {
				if _, ok := pgsToLoci[pgsID]; !ok {
					pgsToLoci[pgsID] = make(map[string]struct{})
				}
				pgsToLoci[pgsID][locus] = struct{}{}
			}
		}

		sort.Slice(allPgs, func(i, j int) bool {
			if len(pgsToLoci[allPgs[i]]) == len(pgsToLoci[allPgs[j]]) {
				return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
			}
			return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
		})
		fmt.Println(allPgs)
		indPop := pgs.GetIndividualAncestry(individual, populations)
		guessedSnps := make(map[string]uint8)
		guessedRefs := make(map[string]string)
		guessConfidence := make(map[string]int)
		trueSnps := make(map[string]uint8) // Count only the snps that were "guessed", i.e., no solutions does not Count
		accuracies := make(map[int]float32)
		var numMatches float32
		thresholdPos := 0
		for {
			if len(allPgs) == 0 {
				break
			}
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats(tools.GG)
			if err != nil {
				log.Printf("Error loading stats: %v\n", err)
				return
			}
			recoveredSnps := make(map[int]uint8)
			recoveredRefs := make(map[int]string)
			includedRefs := make(map[string]struct{})
			for l, locus := range p.Loci {
				if guess, ok := guessedSnps[locus]; ok {
					recoveredSnps[l] = guess
					recoveredRefs[l] = guessedRefs[locus]
					if _, ok = includedRefs[guessedRefs[locus]]; !ok {
						includedRefs[guessedRefs[locus]] = struct{}{}
					}
				}
			}
			fmt.Printf("Total SNPs %d, unknown %d\n", len(p.Loci), len(p.Loci)-len(recoveredSnps))
			fmt.Printf("Recovered snps: %v\n", recoveredSnps)
			fmt.Printf("Recovered references: %v\n", recoveredRefs)
			cohort := solver.NewCohort(p, tools.GG)
			solutions := findSolutions(p, cohort, individual, indPop, recoveredSnps)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				solutions = selfRepair(p, cohort, individual, indPop, guessedSnps, guessedRefs, guessConfidence,
					recoveredSnps, recoveredRefs)
				recoveredSnps = nil
			}
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[individual].Genotype))
				for i, locus := range p.Loci {
					if _, ok := guessedSnps[locus]; ok {
						continue
					}
					guessedSnps[locus] = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
					guessedRefs[locus] = pgsID
				}
				//
				guessConfidence[pgsID] = 1
				fmt.Printf("Increasing guess confidence for ")
				for ref := range includedRefs {
					guessConfidence[ref]++
					// If all or all but one snps have been guessed prior, and they work,
					// we have higher confidence that they are correct
					if len(p.Loci)-len(recoveredSnps) <= 1 {
						guessConfidence[ref]++
					}
					fmt.Printf("%s:%d ", ref, guessConfidence[ref])
				}
				fmt.Println()
			} else {
				fmt.Println("!!!!!!! Still could not find a solution !!!!!!!")
			}
			for i, locus := range p.Loci {
				trueSnps[locus] = cohort[individual].Genotype[pgs.Ploidy*i] + cohort[individual].Genotype[pgs.Ploidy*i+1]
			}
			if thresholdPos < len(thresholds) && len(guessedSnps) >= thresholds[thresholdPos] {
				numMatches = 0
				for locus, guessed := range guessedSnps {
					if trueSnp, ok := trueSnps[locus]; ok && trueSnp == guessed {
						numMatches += 1
					}
				}
				accuracies[thresholds[thresholdPos]] = numMatches / float32(len(guessedSnps))
				fmt.Printf("### CountAccuracy %d: %.3f\n", thresholds[thresholdPos],
					accuracies[thresholds[thresholdPos]])
				thresholdPos++
			}
			fmt.Printf("Guessed %d\n", len(guessedSnps))
			allPgs = allPgs[1:]
			for _, locus := range p.Loci {
				for _, id := range lociToPgs[locus] {
					if id == pgsID {
						continue
					}
					if _, ok := pgsToLoci[id][locus]; ok {
						delete(pgsToLoci[id], locus)
					}
				}
			}
		}
		results = append(results, &Result{Individual: individual, Ancestry: indPop, Accuracies: accuracies})
		guesses = append(guesses, &Guess{Individual: individual, Ancestry: indPop, SNPs: guessedSnps})
	}

	resultFolder := "results/sequential"
	filepath := path.Join(resultFolder, fmt.Sprintf("accuracies%d.json", chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}

	filepath = path.Join(resultFolder, fmt.Sprintf("guesses%d.json", chunkNum))
	guessFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer guessFile.Close()
	encoder = json.NewEncoder(guessFile)
	if err = encoder.Encode(guesses); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func sequenceSolvingAF() {
	thresholds := []int{500, 1000, 1500, 2000, 2500}
	fmt.Printf("Sequential solving for %d\n", thresholds)
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
	file, err = os.Open("results/validated_loci.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return
	}
	defer file.Close()
	decoder = json.NewDecoder(file)
	var lociToPgs map[string][]string
	err = decoder.Decode(&lociToPgs)
	if err != nil {
		log.Println("Error decoding validated loci:", err)
		return
	}

	individuals := getIndividualsSample()
	ppl := os.Args[2]
	fmt.Printf("%s: %s\n", ppl, individuals[ppl])

	type Result struct {
		Individual string
		Ancestry   string
		Accuracies map[int]float32
	}
	type Guess struct {
		Individual string
		Ancestry   string
		SNPs       map[string]uint8
	}
	results := make([]*Result, 0)
	guesses := make([]*Guess, 0)
	populations := tools.LoadAncestry()
	var pgsID string
	var individual string
	for c := 0; c < len(individuals[ppl]); c++ {
		if c >= len(individuals[ppl]) {
			break
		}
		individual = individuals[ppl][c]
		fmt.Printf("--------- %s --------\n", individual)

		allPgs := make([]string, 0, len(idsToNumVariants))
		for id, _ := range idsToNumVariants {
			allPgs = append(allPgs, id)
		}
		pgsToLoci := make(map[string]map[string]struct{})
		for locus, pgsIDs := range lociToPgs {
			for _, pgsID := range pgsIDs {
				if _, ok := pgsToLoci[pgsID]; !ok {
					pgsToLoci[pgsID] = make(map[string]struct{})
				}
				pgsToLoci[pgsID][locus] = struct{}{}
			}
		}

		sort.Slice(allPgs, func(i, j int) bool {
			if len(pgsToLoci[allPgs[i]]) == len(pgsToLoci[allPgs[j]]) {
				return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
			}
			return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
		})
		fmt.Println(allPgs)
		indPop := "ALL"
		guessedSnps := make(map[string]uint8)
		guessedRefs := make(map[string]string)
		guessConfidence := make(map[string]int)
		trueSnps := make(map[string]uint8) // Count only the snps that were "guessed", i.e., no solutions does not Count
		accuracies := make(map[int]float32)
		var numMatches float32
		thresholdPos := 0
		for {
			if len(allPgs) == 0 {
				break
			}
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats(tools.GG)
			if err != nil {
				log.Printf("Error loading stats: %v\n", err)
				return
			}
			recoveredSnps := make(map[int]uint8)
			recoveredRefs := make(map[int]string)
			includedRefs := make(map[string]struct{})
			for l, locus := range p.Loci {
				if guess, ok := guessedSnps[locus]; ok {
					recoveredSnps[l] = guess
					recoveredRefs[l] = guessedRefs[locus]
					if _, ok = includedRefs[guessedRefs[locus]]; !ok {
						includedRefs[guessedRefs[locus]] = struct{}{}
					}
				}
			}
			fmt.Printf("Total SNPs %d, unknown %d\n", len(p.Loci), len(p.Loci)-len(recoveredSnps))
			fmt.Printf("Recovered snps: %v\n", recoveredSnps)
			fmt.Printf("Recovered references: %v\n", recoveredRefs)
			cohort := solver.NewCohort(p, tools.GG)
			solutions := findSolutions(p, cohort, individual, indPop, recoveredSnps)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				solutions = selfRepair(p, cohort, individual, indPop, guessedSnps, guessedRefs, guessConfidence,
					recoveredSnps, recoveredRefs)
				recoveredSnps = nil
			}
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[individual].Genotype))
				for i, locus := range p.Loci {
					if _, ok := guessedSnps[locus]; ok {
						continue
					}
					guessedSnps[locus] = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
					guessedRefs[locus] = pgsID
				}
				//
				guessConfidence[pgsID] = 1
				fmt.Printf("Increasing guess confidence for ")
				for ref := range includedRefs {
					guessConfidence[ref]++
					// If all or all but one snps have been guessed prior, and they work,
					// we have higher confidence that they are correct
					if len(p.Loci)-len(recoveredSnps) <= 1 {
						guessConfidence[ref]++
					}
					fmt.Printf("%s:%d ", ref, guessConfidence[ref])
				}
				fmt.Println()
			} else {
				fmt.Println("!!!!!!! Still could not find a solution !!!!!!!")
			}
			for i, locus := range p.Loci {
				trueSnps[locus] = cohort[individual].Genotype[pgs.Ploidy*i] + cohort[individual].Genotype[pgs.Ploidy*i+1]
			}
			if thresholdPos < len(thresholds) && len(guessedSnps) >= thresholds[thresholdPos] {
				numMatches = 0
				for locus, guessed := range guessedSnps {
					if trueSnp, ok := trueSnps[locus]; ok && trueSnp == guessed {
						numMatches += 1
					}
				}
				accuracies[thresholds[thresholdPos]] = numMatches / float32(len(guessedSnps))
				fmt.Printf("### CountAccuracy %d: %.3f\n", thresholds[thresholdPos],
					accuracies[thresholds[thresholdPos]])
				thresholdPos++
			}
			fmt.Printf("Guessed %d\n", len(guessedSnps))
			allPgs = allPgs[1:]
			for _, locus := range p.Loci {
				for _, id := range lociToPgs[locus] {
					if id == pgsID {
						continue
					}
					if _, ok := pgsToLoci[id][locus]; ok {
						delete(pgsToLoci[id], locus)
					}
				}
			}
		}
		results = append(results, &Result{Individual: individual,
			Ancestry: pgs.GetIndividualAncestry(individual, populations), Accuracies: accuracies})
		guesses = append(guesses, &Guess{Individual: individual,
			Ancestry: pgs.GetIndividualAncestry(individual, populations), SNPs: guessedSnps})
	}

	resultFolder := "results/sequential"
	filepath := path.Join(resultFolder, fmt.Sprintf("af-accuracies-%s.json", ppl))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}

	filepath = path.Join(resultFolder, fmt.Sprintf("af-guesses-%s.json", ppl))
	guessFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer guessFile.Close()
	encoder = json.NewEncoder(guessFile)
	if err = encoder.Encode(guesses); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func getIndividualsSample() map[string][]string {
	type AccuracyInput struct {
		Individual        string  `json:"Individual"`
		Ancestry          string  `json:"Ancestry"`
		SNPs              int     `json:"SNPs"`
		GuessAccuracy     float64 `json:"GuessAccuracy"`
		ReferenceAccuracy float64 `json:"ReferenceAccuracy"`
	}
	accuracyInputs := make([]*AccuracyInput, 0)
	file, err := os.Open("results/guessAccuracy/individuals.json")
	if err != nil {
		log.Println("Error opening validated ids file:", err)
		return nil
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	err = decoder.Decode(&accuracyInputs)
	if err != nil {
		log.Println("Error decoding accuracy inputs:", err)
		return nil
	}

	ancestryGroups := make(map[string][]*AccuracyInput)
	for _, input := range accuracyInputs {
		if _, ok := ancestryGroups[input.Ancestry]; !ok {
			ancestryGroups[input.Ancestry] = make([]*AccuracyInput, 0)
		}
		ancestryGroups[input.Ancestry] = append(ancestryGroups[input.Ancestry], input)
	}
	individuals := make(map[string][]string)
	for ancestry, inputs := range ancestryGroups {
		fmt.Printf("--- %s ---\n", ancestry)
		sort.Slice(inputs, func(i, j int) bool {
			return inputs[i].GuessAccuracy > inputs[j].GuessAccuracy
		})
		individuals[ancestry] = make([]string, 0)
		step := len(inputs) / 10
		for i := 0; i < 10; i++ {
			index := i * step
			if index >= len(inputs) {
				index = len(inputs) - 1
			}
			individuals[ancestry] = append(individuals[ancestry], inputs[index].Individual)
			fmt.Printf("%s: %s\n", inputs[index].Individual, strconv.FormatFloat(inputs[index].GuessAccuracy, 'f', 3, 64))
		}
	}
	return individuals
}

func selfRepair(p *pgs.PGS, cohort solver.Cohort, individual string, indPop string, guessedSnps map[string]uint8,
	guessedRefs map[string]string, guessConfidence map[string]int,
	recoveredSnps map[int]uint8, recoveredRefs map[int]string) [][]uint8 {
	fmt.Println("^^^ No solutions with all the extra loci ^^^")
	highConfidenceSnps := make(map[int]uint8)
	highConfidenceRefs := make(map[int]string)
	for i, snp := range recoveredSnps {
		if guessConfidence[recoveredRefs[i]] > ConfidenceThreshold {
			highConfidenceSnps[i] = snp
			highConfidenceRefs[i] = recoveredRefs[i]
		}
	}
	// Too hard to solve from scratch
	if p.NumVariants-len(highConfidenceSnps) > ScratchSolvingSnpLimit {
		fmt.Printf("Too few high-confidence SNPs (%d/%d), skipping\n", len(highConfidenceSnps), p.NumVariants)
		return nil
	}
	if p.PgsID == "PGS000648" && p.NumVariants-len(highConfidenceSnps) > DeterminismLimit {
		fmt.Printf("Ignoring %s, only %d/%d SNPs are known\n", p.PgsID, len(highConfidenceSnps), p.NumVariants)
		return nil
	}
	fmt.Printf("Solving with high-confidence SNPs: %v\n", highConfidenceRefs)
	solutions := findSolutions(p, cohort, individual, indPop, highConfidenceSnps)
	highConfidence := true
	if len(solutions) == 0 {
		if p.NumVariants > DeterminismLimit {
			return nil
		} else {
			highConfidence = false
			solutions = findSolutions(p, cohort, individual, indPop, make(map[int]uint8))
		}
	}
	if len(solutions) == 0 {
		return nil
	}
	fmt.Printf("Top solution with high-confidence snps/from scratch: accuracy %.2f, likelihood %.2f\n",
		solver.Accuracy(solutions[0], cohort[individual].Genotype),
		solver.CalculateFullSequenceLikelihood(solutions[0], p.PopulationStats[indPop].AF, p.EffectAlleles))

	refMap := make(map[string]struct{})
	refs := make([]string, 0)
	for _, ref := range recoveredRefs {
		if guessConfidence[ref] > ConfidenceThreshold && highConfidence {
			//if guessConfidence[ref] > ConfidenceThreshold {
			fmt.Printf("Already confident about %s\n", ref)
			continue
		}
		if _, ok := refMap[ref]; !ok {
			refMap[ref] = struct{}{}
			refs = append(refs, ref)
			guessConfidence[ref] -= 1
		}
	}

	// Resolving the references to see if there is a valid solution with the new SNP results
	fmt.Println("Resolving mismatching references")
	limit := 30
	if len(solutions) < limit {
		limit = len(solutions)
	}
	var allMatching bool
	var err error
solutionLoop:
	for k, solution := range solutions[:limit] {
		fmt.Printf("Trying a solution %d with accuracy %.3f\n", k, solver.Accuracy(solution, cohort[individual].Genotype))
		refSols := make([][][]uint8, len(refs))
		for j, ref := range refs {
			fmt.Printf("--- %s ---\n", ref)
			allMatching = true
			for i := range p.Loci {
				if recoveredRefs[i] == ref {
					if guess, ok := recoveredSnps[i]; ok && solutions[0][pgs.Ploidy*i]+solutions[0][pgs.Ploidy*i+1] != guess {
						allMatching = false
					}
				}
			}
			if allMatching {
				fmt.Println("No mismatching SNPs")
				continue
			}
			refp := pgs.NewPGS()
			err = refp.LoadCatalogFile(path.Join(params.LocalDataFolder, ref+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file %s: %v\n", ref, err)
				return nil
			}
			err = refp.LoadStats(tools.GG)
			if err != nil {
				log.Printf("Error loading stats for %s: %v\n", ref, err)
				return nil
			}
			refCohort := solver.NewCohort(refp, tools.GG)
			newSnps := make(map[int]uint8)
			for i, locus := range p.Loci {
				if recoveredRefs[i] == ref {
					for j, refLocus := range refp.Loci {
						if refLocus == locus {
							newSnps[j] = solution[pgs.Ploidy*i] + solution[pgs.Ploidy*i+1]
							break
						}
					}
				}
			}
			if refp.NumVariants-len(newSnps) > DeterminismLimit {
				continue
			}
			refSols[j] = findSolutions(refp, refCohort, individual, indPop, newSnps)
			if len(refSols[j]) == 0 {
				fmt.Printf("Could not resolve %s with new SNPs\n", ref)
				continue solutionLoop
			}
			fmt.Printf("New %s accuracy: %.3f\n", ref, solver.Accuracy(refSols[j][0], refCohort[individual].Genotype))
		}
		fmt.Printf("Solution %d is valid\n", k)
		// All refs are solvable, updating old loci
		for j, ref := range refs {
			if len(refSols[j]) == 0 {
				continue
			}
			refp := pgs.NewPGS()
			err = refp.LoadCatalogFile(path.Join(params.LocalDataFolder, ref+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file %s: %v\n", ref, err)
				return nil
			}
			for i, refloc := range refp.Loci {
				if guessConfidence[guessedRefs[refloc]] > ConfidenceThreshold {
					continue
				}
				guessedSnps[refloc] = refSols[j][0][pgs.Ploidy*i] + refSols[j][0][pgs.Ploidy*i+1]
				guessedRefs[refloc] = ref
			}
			guessConfidence[ref] = 0
		}
		solutions[0] = solution
		break solutionLoop
	}
	return solutions
}

func findSolutions(p *pgs.PGS, cht solver.Cohort, idv, ppl string, priorSnps map[int]uint8) [][]uint8 {
	var sm map[string][]uint8
	var precision uint32 = 9
	if len(p.Loci)-len(priorSnps) >= DeterminismLimit {
		precision = 8
	}
	slv := solver.NewDP(cht[idv].Score, p, ppl, precision, priorSnps)
	if len(p.Loci)-len(priorSnps) < DeterminismLimit {
		sm = slv.SolveDeterministic(solver.UseLikelihood)
	} else {
		sm = slv.SolveProbabilistic(solver.UseLikelihood)
	}
	return solver.SortByLikelihoodAndFrequency(sm, p.PopulationStats[ppl], p.EffectAlleles, solver.UseLikelihood)
}

type accuracyOutput struct {
	Individual           string
	Ancestry             string
	Score                string
	TrueLikelihood       string
	ReferenceAccuracy    string
	LikelihoodAccuracies []string
}

func newAccuracyOutput(ind string, pop string, score *apd.Decimal, trueLikelihood, refAcc float32) *accuracyOutput {
	fscore, err := score.Float64()
	if err != nil {
		log.Fatalf("Error converting score to float64: %v", err)
	}
	return &accuracyOutput{
		Individual:        ind,
		Ancestry:          pop,
		Score:             fmt.Sprintf("%.3f", fscore),
		TrueLikelihood:    fmt.Sprintf("%.3f", trueLikelihood),
		ReferenceAccuracy: fmt.Sprintf("%.3f", refAcc),
	}
}

func accuracyParallel() {
	resultFolder := "results/accuracy"
	output := make([]*accuracyOutput, 0)
	catalogFile := "PGS001835_hmPOS_GRCh37.txt"
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(path.Join(params.LocalDataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("======== %s ========\n", p.PgsID)
	p.LoadStats(tools.GG)
	cohort := solver.NewCohort(p, tools.GG)
	populations := tools.LoadAncestry()
	individuals := solver.All1000GenomesSamples()

	chunkNum, chunkSize := getChunkInfo(len(individuals))
	fmt.Println(individuals[chunkNum*chunkSize : (chunkNum+1)*chunkSize])
	filepath := path.Join(resultFolder, fmt.Sprintf("%s-%d.json", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	var individual string
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(individuals) {
			break
		}
		individual = individuals[i]
		fmt.Printf("%s\n", individual)
		indPop := populations[individual]
		if strings.Contains(indPop, ",") {
			indPop = strings.Split(indPop, ",")[0]
		}
		majorReference := solver.AllReferenceAlleleSample(p.PopulationStats[indPop].AF)
		out := newAccuracyOutput(individual, indPop, cohort[individual].Score,
			solver.CalculateFullSequenceLikelihood(cohort[individual].Genotype, p.PopulationStats[indPop].AF, p.EffectAlleles),
			solver.Accuracy(majorReference, cohort[individual].Genotype))

		solutions := findSolutions(p, cohort, individual, indPop, make(map[int]uint8))
		out.LikelihoodAccuracies = make([]string, 0, len(solutions))
		for _, solution := range solutions {
			out.LikelihoodAccuracies = append(out.LikelihoodAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
		}
		output = append(output, out)
	}
	if err = encoder.Encode(output); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func findAllSolutions() {
	//INDIVIDUAL := "NA18595"
	//INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	//INDIVIDUAL := "HG00551" // low 648
	INDIVIDUAL := "HG02024"

	//INDIVIDUAL := "HG01028"
	//INDIVIDUAL := "NA18531"
	//INDIVIDUAL := "NA20872"
	//INDIVIDUAL := "NA20507"
	//INDIVIDUAL := "NA07037"
	//INDIVIDUAL := "HG03015"
	//INDIVIDUAL := "HG01253"
	//INDIVIDUAL := "NA19313"

	p := pgs.NewPGS()
	//catalogFile := "PGS000154_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000753_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000043_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000786_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh37.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000307_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000011_hmPOS_GRCh37.txt"
	//catalogFile := "PGS003436_hmPOS_GRCh37.txt"
	//catalogFile := "PGS002264_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004238_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004222_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004249_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000119_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000880_hmPOS_GRCh37.txt"
	catalogFile := "PGS000083_hmPOS_GRCh37.txt"
	err := p.LoadCatalogFile(path.Join(params.LocalDataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	err = p.LoadStats(tools.GG)
	if err != nil {
		log.Printf("Error loading stats: %v\n", err)
		return
	}
	cohort := solver.NewCohort(p, tools.GG)

	populations := tools.LoadAncestry()
	idvAnc := pgs.GetIndividualAncestry(INDIVIDUAL, populations)
	//recovered := map[int]uint8{0: 0, 1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 1, 8: 2, 10: 1, 11: 0}
	recovered := make(map[int]uint8)
	solutions := findSolutions(p, cohort, INDIVIDUAL, idvAnc, recovered)

	fmt.Printf("\nTrue:\n%s, %.2f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
		solver.CalculateFullSequenceLikelihood(cohort[INDIVIDUAL].Genotype, p.PopulationStats[idvAnc].AF, p.EffectAlleles))
	//solver.CalculateSolutionSpectrumDistance(cohort[INDIVIDUAL].Genotype, p.PopulationStats[idvAnc], p.EffectAlleles))

	fmt.Printf("Guessed %d:\n", len(solutions))
	//Target := solver.ScoreToTarget(cohort[INDIVIDUAL].Score, p)
	target := cohort[INDIVIDUAL].Score
	for _, solution := range solutions {
		diff := new(apd.Decimal)
		p.Context.Sub(diff, target, solver.CalculateDecimalScore(p.Context, solution, p.Weights, p.EffectAlleles))
		fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
			diff.String(), solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[idvAnc].AF, p.EffectAlleles))
		//	solver.CalculateSolutionSpectrumDistance(solution, p.PopulationStats[idvAnc], p.EffectAlleles))
	}
}

type sortingOutput struct {
	Individual                   string
	Ancestry                     string
	Score                        string
	TrueLikelihood               string
	TrueSpectrum                 string
	ReferenceAccuracy            string
	LikelihoodAccuracies         []string
	SpectrumAccuracies           []string
	LikelihoodSpectrumAccuracies []string
}

func newSortingOutput(ind string, pop string, score *apd.Decimal, trueLikelihood, trueSpectrum, refAcc float32) *sortingOutput {
	fscore, err := score.Float64()
	if err != nil {
		log.Fatalf("Error converting score to float64: %v", err)
	}
	return &sortingOutput{
		Individual:        ind,
		Ancestry:          pop,
		Score:             fmt.Sprintf("%.3f", fscore),
		TrueLikelihood:    fmt.Sprintf("%.3f", trueLikelihood),
		TrueSpectrum:      fmt.Sprintf("%.3f", trueSpectrum),
		ReferenceAccuracy: fmt.Sprintf("%.3f", refAcc),
	}
}

//func sortingChoiceParallel() {
//	DeterminismLimit := 35
//	resultFolder := "results/sorting"
//	output := make([]*sortingOutput, 0)
//	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000753_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000154_hmPOS_GRCh37.txt"
//	catalogFile := "PGS003760_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000851_hmPOS_GRCh37.txt"
//	samplesPerPopulation := 10
//	p := pgs.NewPGS()
//	err := p.LoadCatalogFile(path.Join(params.LocalDataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	fmt.Printf("======== %s ========\n", p.PgsID)
//	p.LoadStats()
//	cohort := solver.NewCohort(p)
//	populations := tools.LoadAncestry()
//	sortedSamples := cohort.SortByScore()
//	sortedScores := make([]float64, len(sortedSamples))
//	for i, sample := range sortedSamples {
//		sortedScores[i], _ = cohort[sample].Score.Float64()
//	}
//	minScore, maxScore := sortedScores[0], sortedScores[len(sortedScores)-1]
//	step := (maxScore - minScore) / float64(samplesPerPopulation)
//	individuals := make([]string, 0, samplesPerPopulation*len(pgs.ANCESTRIES))
//	var Count int
//	for _, ppl := range pgs.ANCESTRIES {
//		score := minScore
//		Count = 0
//		for i := 0; i < len(sortedSamples); i++ {
//			if sortedScores[i] < score || populations[sortedSamples[i]] != ppl || sortedScores[i] == 0 {
//				continue
//			}
//			individuals = append(individuals, sortedSamples[i])
//			score += step
//			Count += 1
//			if Count == samplesPerPopulation {
//				break
//			}
//		}
//	}
//	chunkNum, chunkSize := getChunkInfo(len(individuals))
//	fmt.Println(individuals[chunkNum*chunkSize : (chunkNum+1)*chunkSize])
//	filepath := path.Join(resultFolder, fmt.Sprintf("%s-%d.json", p.PgsID, chunkNum))
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	encoder := json.NewEncoder(resFile)
//	var individual string
//	var solmap map[string][]uint8
//	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
//		if i >= len(individuals) {
//			break
//		}
//		individual = individuals[i]
//		fmt.Printf("%s\n", individual)
//		indPop := populations[individual]
//		if strings.Contains(indPop, ",") {
//			indPop = strings.Split(indPop, ",")[0]
//		}
//		majorReference := solver.AllReferenceAlleleSample(p.PopulationStats[indPop].AF)
//		out := newSortingOutput(individual, indPop, cohort[individual].Score,
//			solver.CalculateFullSequenceLikelihood(cohort[individual].Genotype, p.PopulationStats[indPop].AF, p.EffectAlleles),
//			solver.CalculateSolutionSpectrumDistance(cohort[individual].Genotype, p.PopulationStats[indPop], p.EffectAlleles),
//			solver.CountAccuracy(majorReference, cohort[individual].Genotype))
//
//		slv := solver.NewDP(cohort[individual].Score, p, indPop)
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveDeterministic(solver.UseLikelihood)
//		} else {
//			solmap = slv.SolveProbabilistic(solver.UseLikelihood)
//		}
//		solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
//		out.LikelihoodAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodAccuracies = append(out.LikelihoodAccuracies, fmt.Sprintf("%.3f", solver.CountAccuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveDeterministic(solver.UseSpectrum)
//		} else {
//			solmap = slv.SolveProbabilistic(solver.UseSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseSpectrum)
//		out.SpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.SpectrumAccuracies = append(out.SpectrumAccuracies, fmt.Sprintf("%.3f", solver.CountAccuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveDeterministic(solver.UseLikelihoodAndSpectrum)
//		} else {
//			solmap = slv.SolveProbabilistic(solver.UseLikelihoodAndSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihoodAndSpectrum)
//		out.LikelihoodSpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodSpectrumAccuracies = append(out.LikelihoodSpectrumAccuracies, fmt.Sprintf("%.3f", solver.CountAccuracy(solution, cohort[individual].Genotype)))
//		}
//		output = append(output, out)
//	}
//	if err = encoder.Encode(output); err != nil {
//		log.Fatal("Cannot encode json", err)
//	}
//}

func getChunkInfo(totalLen int) (int, int) {
	var err error
	var chunkNum, chunkSize = 0, 0
	if len(os.Args) > 3 {
		chunkNum, err = strconv.Atoi(os.Args[2])
		if err != nil {
			log.Fatalf("Error parsing chunkNum %s: %v", os.Args[2], err)
		}
		chunkSize, err = strconv.Atoi(os.Args[3])
		if err != nil {
			log.Fatalf("Error parsing chunkSize %s: %v", os.Args[3], err)
		}
	} else {
		chunkSize = totalLen
	}
	return chunkNum, chunkSize
}

func getNumThreads() int {
	var num int
	var err error
	numThreadsEnv := os.Getenv(params.NumCpusEnv)
	if numThreadsEnv != "" {
		num, err = strconv.Atoi(numThreadsEnv)
		if err != nil {
			log.Fatalf("Error parsing numCpus %s: %v\n", params.NumCpusEnv, err)
		}
	} else {
		//log.Printf("Could not read numCpus env %s\n", params.NumCpusEnv)
		num = 16
	}
	return num
}
