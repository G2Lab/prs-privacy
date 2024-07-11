package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path"
	"path/filepath"
	"regexp"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
var memprofile = flag.String("memprofile", "", "write cpu profile to file")

const DeterminismLimit = 34

func main() {
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

	//scoreDistribution()
	//likelihoodEffect()
	//scoreToLikelihoodDistribution()
	//scoreToLikelihood()
	//evaluateReferences()
	//accuracyLikelihood()
	//buildDPTables()
	//likelihoodWeight()
	//sortingChoice()
	//sortingChoiceParallel()
	//accuracyParallel()
	findAllSolutions()
	//kinshipExperiment()
	//kingTest()
	//consensusSolving()
	//uniquenessExperiment()
	//sequenceSolving()
	//calculateGenotypeFrequencies()
	//imputeWorkflow()
	//linkingWithImputation()
	//linkingWithGuessed()
	//evaluateImputation()
	//findUnsolvablePRSWithOverlap()
	//predictPRS()
}

func predictPRS() {
	type Result struct {
		Known     int
		Total     int
		Predicted []float64
		Real      []float64
	}
	results := make(map[string]*Result)
	guessed := loadGuessedGenotypes()
	//catalogFile := "PGS000055_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000031_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000515_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000148_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000077_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000788_hmPOS_GRCh37.txt"
	catalogFiles := []string{"PGS000055_hmPOS_GRCh37.txt", "PGS000031_hmPOS_GRCh37.txt", "PGS000515_hmPOS_GRCh37.txt",
		"PGS000148_hmPOS_GRCh37.txt", "PGS000720_hmPOS_GRCh37.txt"}
	for _, catalogFile := range catalogFiles {
		p := pgs.NewPGS()
		err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		p.LoadStats()
		results[p.PgsID] = &Result{
			Known:     0,
			Total:     len(p.Loci),
			Predicted: make([]float64, 0),
			Real:      make([]float64, 0),
		}
		cohort := solver.NewCohort(p)
		scores := make([]float64, 0)
		var scf float64
		for idv := range cohort {
			scf, err = cohort[idv].Score.Float64()
			if err != nil {
				log.Printf("Error converting predictedScore to float64: %v\n", err)
				return
			}
			scores = append(scores, scf)
		}
		sort.Slice(scores, func(i, j int) bool {
			return scores[i] < scores[j]
		})
		weights := make([]float64, len(p.Weights))
		for i, w := range p.Weights {
			weights[i], err = w.Float64()
			if err != nil {
				log.Printf("Error converting weight to float64: %v\n", err)
				return
			}
		}
		minScore, maxScore := scores[0], scores[len(scores)-1]
		populations := tools.LoadAncestry()
		individuals := solver.AllRelativeSamples()
		var predictedScore, trueScore float64
		var majorFreq float32
		var guess, majorGtp uint8
		var ok bool
		for k, idv := range individuals {
			predictedScore = 0
			trueScore, err = cohort[idv].Score.Float64()
			if err != nil {
				log.Printf("Error converting trueScore to float64: %v\n", err)
				return
			}
			for i, locus := range p.Loci {
				chr, pos := tools.SplitLocus(locus)
				guess, ok = guessed[idv][chr][pos]
				if !ok {
					// we do not have a guess for the locus, so we assume the most common genotype
					majorGtp = 0
					majorFreq = p.PopulationStats[pgs.GetIndividualAncestry(idv, populations)].GF[i][0]
					for j, freq := range p.PopulationStats[pgs.GetIndividualAncestry(idv, populations)].GF[i] {
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
					predictedScore += 2 * weights[i]
				case guess == 1:
					predictedScore += weights[i]
				default:
					if guess != 0 && guess != 1 && guess != 2 {
						log.Printf("Unknown guess genotype at %s: %d\n", locus, guess)
					}
				}
			}
			results[p.PgsID].Predicted = append(results[p.PgsID].Predicted, percentile(predictedScore, minScore, maxScore))
			results[p.PgsID].Real = append(results[p.PgsID].Real, percentile(trueScore, minScore, maxScore))
			fmt.Printf("%s: %.3f%% -> %.3f%%\n", idv, percentile(trueScore, minScore, maxScore),
				percentile(predictedScore, minScore, maxScore))
		}
		guessedSumPos, guessedSumNeg := 0.0, 0.0
		restSumPos, restSumNeg := 0.0, 0.0
		for _, idv := range individuals {
			for i, locus := range p.Loci {
				chr, pos := tools.SplitLocus(locus)
				if _, ok = guessed[idv][chr][pos]; ok {
					if weights[i] > 0 {
						guessedSumPos += weights[i]
					} else {
						guessedSumNeg += weights[i]
					}
				} else {
					if weights[i] > 0 {
						restSumPos += weights[i]
					} else {
						restSumNeg += weights[i]
					}
				}
			}
			break
		}
		fmt.Printf("Guessed: %.3f, %.3f\n", guessedSumPos, guessedSumNeg)
		fmt.Printf("Rest: %.3f, %.3f\n", restSumPos, restSumNeg)
	}
	//filePath := path.Join("results/predict", "prediction.json")
	//resFile, err := os.OpenFile(filePath, os.O_CREATE|os.O_WRONLY, 0644)
	//if err != nil {
	//	log.Fatalf("Error opening result file: %v", err)
	//}
	//defer resFile.Close()
	//encoder := json.NewEncoder(resFile)
	//if err = encoder.Encode(results); err != nil {
	//	log.Fatal("Cannot encode json", err)
	//}
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

func findUnsolvablePRSWithOverlap() {
	guessedPerIndividual := loadGuessedGenotypes()
	guessedLoci := make(map[string]struct{})
	for idv := range guessedPerIndividual {
		for chr := range guessedPerIndividual[idv] {
			for pos := range guessedPerIndividual[idv][chr] {
				guessedLoci[fmt.Sprintf("%s:%s", chr, pos)] = struct{}{}
			}
		}
	}
	ids, err := fewerVariantsPGS(50, 100)
	if err != nil {
		log.Println("Error:", err)
		return
	}
	catalogFolder := "catalog"
	unknownPgs := make(map[string][]int)
	for _, id := range ids {
		fmt.Printf("==== %s ====\n", id)
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(filepath.Join(catalogFolder, id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Println("Error:", err)
			continue
		}
		for _, locus := range p.Loci {
			if strings.HasPrefix(locus, "X:") || strings.HasPrefix(locus, "Y:") {
				continue
			}
			if !isNumber(strings.Split(locus, ":")[1]) || !isNumberInRange(strings.Split(locus, ":")[0], 1, 22) {
				continue
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
	pgs := make([]string, 0)
	for id, counts := range unknownPgs {
		if counts[0] > 0 {
			pgs = append(pgs, id)
		}
	}
	sort.Slice(pgs, func(i, j int) bool {
		return float64(unknownPgs[pgs[i]][0])/float64(unknownPgs[pgs[i]][1]) >
			float64(unknownPgs[pgs[j]][0])/float64(unknownPgs[pgs[j]][1])
	})
	for _, id := range pgs {
		fmt.Printf("%s: %d/%d\n", id, unknownPgs[id][0], unknownPgs[id][1])
	}
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

type Result struct {
	Individual  string
	Score       string
	Accuracies  []string
	Likelihoods []string
	Spectrums   []string
}

func NewResult(ind string, score *apd.Decimal) *Result {
	return &Result{
		Individual:  ind,
		Score:       score.String(),
		Accuracies:  make([]string, 0),
		Likelihoods: make([]string, 0),
		Spectrums:   make([]string, 0),
	}
}

type Reference struct {
	ID         string
	Score      string
	Likelihood string
	Accuracies []string
}

func NewReference(id string, score *apd.Decimal, likelihood float64) *Reference {
	return &Reference{
		ID:         id,
		Score:      score.String(),
		Likelihood: fmt.Sprintf("%.3f", likelihood),
		Accuracies: make([]string, 0),
	}
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
		err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		err = p.LoadStats()
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return
		}
		cohort := solver.NewCohort(p)
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
		var solmap map[string][]uint8
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
			err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats()
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
			cohort := solver.NewCohort(p)
			for _, locus := range p.Loci {
				for idv := range db {
					if _, ok := db[idv][locus]; ok {
						continue
					}
					db[idv][locus] = cohort[idv].Genotype[0] + cohort[idv].Genotype[1]
				}
			}
			slv := solver.NewDP(cohort[relative].Score, p, indPop, recoveredSnps)
			if len(p.Loci)-len(recoveredSnps) < DeterminismLimit {
				solmap = slv.SolveDeterministic(solver.UseLikelihood)
			} else {
				solmap = slv.SolveProbabilistic(solver.UseLikelihood)
			}
			solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				fmt.Println("^^^ No solutions with extra loci, calculating from scratch ^^^")
				slv := solver.NewDP(cohort[relative].Score, p, indPop, make(map[int]uint8))
				if len(p.Loci) < DeterminismLimit {
					solmap = slv.SolveDeterministic(solver.UseLikelihood)
				} else {
					solmap = slv.SolveProbabilistic(solver.UseLikelihood)
				}
				solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
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

func consensusSolving() {
	threshold := 2000
	fmt.Printf("Consensus solving for %d\n", threshold)
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

	individuals := []string{"NA19678"}
	//individuals := []string{"HG00119", "HG00524", "HG00581", "HG00656", "HG00731", "HG01936", "HG02025", "HG02026",
	//	"HG02067", "HG02353", "HG02371", "HG02250", "HG02373", "HG02386", "HG02375", "HG03713", "HG03673", "NA19238",
	//	"NA19239", "NA19027", "NA19334", "NA19331", "NA19664", "NA19678", "NA19661", "NA19713", "NA20321", "NA20334",
	//	"NA20289", "NA20792", "NA20868", "NA20895"}

	chunkNum, chunkSize := getChunkInfo(len(individuals))
	fmt.Println(individuals[chunkNum*chunkSize : (chunkNum+1)*chunkSize])

	type Result struct {
		Individual string
		Accuracies map[int]float32
	}
	//results := make([]*Result, 0)
	var individual string
	var ok bool
	var guessedGtp uint8
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
		populations := tools.LoadAncestry()
		indPop := populations[individual]
		if strings.Contains(indPop, ",") {
			indPop = strings.Split(indPop, ",")[0]
		}
		var pgsID string
		var solmap map[string][]uint8
		recoveredSnps := make(map[int]uint8)
		guessedSnps := make(map[string]map[uint8]map[string]struct{})
		trueSnps := make(map[string]uint8) // Count all the snps
		var solutions [][]uint8
		for {
			if len(allPgs) == 0 {
				break
			}
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats()
			if err != nil {
				log.Printf("Error loading stats: %v\n", err)
				return
			}
			lociToPgs[p.PgsID] = p.Loci
			cohort := solver.NewCohort(p)
			slv := solver.NewDP(cohort[individual].Score, p, indPop, recoveredSnps)
			if len(p.Loci) > DeterminismLimit {
				break
			}
			solmap = slv.SolveDeterministic(solver.UseLikelihood)
			solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[individual].Genotype))
				for i, locus := range p.Loci {
					if _, ok = guessedSnps[locus]; !ok {
						guessedSnps[locus] = make(map[uint8]map[string]struct{})
					}
					guessedGtp = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
					if _, ok = guessedSnps[locus][guessedGtp]; !ok {
						guessedSnps[locus][guessedGtp] = make(map[string]struct{})
					}
					guessedSnps[locus][guessedGtp][p.PgsID] = struct{}{}
				}
			}
			for i, locus := range p.Loci {
				if _, ok := trueSnps[locus]; ok {
					continue
				}
				trueSnps[locus] = cohort[individual].Genotype[pgs.Ploidy*i] + cohort[individual].Genotype[pgs.Ploidy*i+1]
			}
			fmt.Printf("Guessed %d\n", len(guessedSnps))
			if len(guessedSnps) >= threshold {
				break
			}
			allPgs = allPgs[1:]
			//for _, locus := range p.Loci {
			//	for _, id := range lociToPgs[locus] {
			//		if id == pgsID {
			//			continue
			//		}
			//		if _, ok := pgsToLoci[id][locus]; ok {
			//			delete(pgsToLoci[id], locus)
			//		}
			//	}
			//}
		}
		//results = append(results, &Result{Individual: individual, Accuracies: accuracies})
		confirmedSnps := make(map[string]uint8)
		for locus, options := range guessedSnps {
			var winGtp uint8
			var winNum = 0
			for gtp, pgsIDs := range options {
				if len(pgsIDs) > winNum {
					winNum = len(pgsIDs)
					winGtp = gtp
				}
			}
			for gtp, pgsIDs := range options {
				if gtp == winGtp {
					continue
				}
				if winNum == len(pgsIDs) {
					// Tie breaking
					winTotalSnps, cndTotalSnps := 0, 0
					for pgs := range options[winGtp] {
						winTotalSnps += idsToNumVariants[pgs]
					}
					for pgsID := range pgsIDs {
						cndTotalSnps += idsToNumVariants[pgsID]
					}
					if cndTotalSnps < winTotalSnps {
						winGtp = gtp
					}
				}
			}
			confirmedSnps[locus] = winGtp
		}
		var accuracy float32 = 0.0
		for locus, snp := range trueSnps {
			if guessed, ok := confirmedSnps[locus]; ok && snp == guessed {
				accuracy += 1
			}
		}
		fmt.Printf("### CountAccuracy after determinism (%d SNPs): %.3f\n", len(confirmedSnps), accuracy/float32(len(trueSnps)))
	}
	//resultFolder := "results/sequential"
	//filepath := path.Join(resultFolder, fmt.Sprintf("chunk%d.json", chunkNum))
	//resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	//if err != nil {
	//	log.Fatalf("Error opening result file: %v", err)
	//}
	//defer resFile.Close()
	//encoder := json.NewEncoder(resFile)
	//if err = encoder.Encode(results); err != nil {
	//	log.Fatal("Cannot encode json", err)
	//}
}

func uniquenessExperiment() {
	fmt.Printf("--- Uniqueness experiment ---\n")
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

	//allPgs = []string{"PGS004164"}
	populations := tools.LoadAncestry()
	type Result struct {
		PgsID           string
		NumVariants     int
		NumUniqueScores int
		AnonymitySets   []int
	}
	results := make([]*Result, 0)
	var sc, idvPop string
	var numUnique int
	var scores map[string][]string
	var scoreList []string
	for _, pgsID := range allPgs {
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		err = p.LoadStats()
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return
		}
		cohort := solver.NewCohort(p)
		scores = make(map[string][]string)
		for idv := range cohort {
			sc = cohort[idv].Score.String()
			if _, ok := scores[sc]; !ok {
				scores[sc] = make([]string, 0)
			}
			scores[sc] = append(scores[sc], idv)
			scoreList = append(scoreList, sc)
			// Check for impossible SNPs
			idvPop = populations[idv]
			if strings.Contains(idvPop, ",") {
				idvPop = strings.Split(idvPop, ",")[0]
			}
			for j, af := range p.PopulationStats[idvPop].AF {
				if (af[0] == 0 && cohort[idv].Genotype[pgs.Ploidy*j]+cohort[idv].Genotype[pgs.Ploidy*j+1] == 0) ||
					(af[1] == 0 && cohort[idv].Genotype[pgs.Ploidy*j]+cohort[idv].Genotype[pgs.Ploidy*j+1] > 0) {
					fmt.Printf("##### %s: Impossible SNP -- idv %s, locus %s (%d), AF %v, genotype %v\n", pgsID, idv,
						p.Loci[j], j, af, cohort[idv].Genotype[pgs.Ploidy*j:pgs.Ploidy*j+2])

				}
			}
		}
		numUnique = 0
		anonsets := make([]int, 0)
		for _, idvs := range scores {
			if len(idvs) == 1 {
				numUnique++
			}
			anonsets = append(anonsets, len(idvs))
		}
		results = append(results, &Result{
			PgsID:           pgsID,
			NumVariants:     pgsToNumVariants[pgsID],
			NumUniqueScores: numUnique,
			AnonymitySets:   anonsets,
		})
	}
	resultFolder := "results/uniqueness"
	filepath := path.Join(resultFolder, "scores.json")
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

type Contender struct {
	id         string
	solution   []uint8
	likelihood float32
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

	//individuals := []string{"HG00119", "HG00524", "HG00581", "HG00656", "HG00731", "HG01936", "HG02025", "HG02026",
	//	"HG02067", "HG02353", "HG02371", "HG02250", "HG02373", "HG02386", "HG02375", "HG03713", "HG03673", "NA19238",
	//	"NA19239", "NA19027", "NA19334", "NA19331", "NA19664", "NA19678", "NA19661", "NA19713", "NA20321", "NA20334",
	//	"NA20289", "NA20792", "NA20868", "NA20895"}
	related := solver.ReadRelatedIndividuals()
	individuals := make([]string, 0)
	for idv := range related {
		individuals = append(individuals, idv)
	}
	sort.Strings(individuals)

	chunkNum, chunkSize := getChunkInfo(len(individuals))
	fmt.Println(individuals[chunkNum*chunkSize : (chunkNum+1)*chunkSize])

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
		indPop := pgs.GetIndividualAncestry(individual, populations)
		guessedSnps := make(map[string]uint8)
		guessConfidence := make(map[string]int)
		recoveredPgsReferences := make(map[string]string)
		trueSnps := make(map[string]uint8) // Count only the snps that were "guessed", i.e., no solutions does not Count
		accuracies := make(map[int]float32)
		var numMatches float32
		thresholdPos := 0
		for {
			if len(allPgs) == 0 {
				break
			}
			//sort.Slice(allPgs, func(i, j int) bool {
			//	if len(pgsToLoci[allPgs[i]]) == len(pgsToLoci[allPgs[j]]) {
			//		return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
			//	}
			//	return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
			//})
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats()
			if err != nil {
				log.Printf("Error loading stats: %v\n", err)
				return
			}
			recoveredSnps := make(map[int]uint8)
			recoveredRefs := make(map[int]string)
			guessedRefs := make(map[string]struct{})
			for l, locus := range p.Loci {
				if guess, ok := guessedSnps[locus]; ok {
					recoveredSnps[l] = guess
					recoveredRefs[l] = recoveredPgsReferences[locus]
					if _, ok = guessedRefs[recoveredPgsReferences[locus]]; !ok {
						guessedRefs[recoveredPgsReferences[locus]] = struct{}{}
					}
				}
			}
			fmt.Printf("Total SNPs %d, unknown %d\n", len(p.Loci), len(p.Loci)-len(recoveredSnps))
			fmt.Printf("Recovered snps: %v\n", recoveredSnps)
			fmt.Printf("Recovered references: %v\n", recoveredRefs)
			cohort := solver.NewCohort(p)
			solutions := findSolutions(p, cohort, individual, indPop, recoveredSnps)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				solutions = selfRepair(p, cohort, individual, indPop, guessedSnps, guessedRefs, guessConfidence, recoveredSnps, recoveredPgsReferences)
			}
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[individual].Genotype))
				for i, locus := range p.Loci {
					if _, ok := guessedSnps[locus]; ok {
						continue
					}
					guessedSnps[locus] = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
					recoveredPgsReferences[locus] = pgsID
				}
				//
				guessConfidence[pgsID] = 1
				fmt.Printf("Increasing guess confidence for ")
				for ref := range guessedRefs {
					guessConfidence[ref]++
					// If all or all but one snps have been guessed prior, and they work,
					// we have higher confidence that they are correct
					if len(p.Loci)-len(recoveredSnps) < 2 {
						guessConfidence[ref] += 5
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

func selfRepair(p *pgs.PGS, cohort solver.Cohort, individual string, indPop string, guessedSnps map[string]uint8,
	guessedRefs map[string]struct{}, guessConfidence map[string]int,
	recoveredSnps map[int]uint8, recoveredPgsReferences map[string]string) [][]uint8 {
	fmt.Println("^^^ No solutions with extra loci, trying with a subset")
	var err error
	ids := make([]string, 0)
	for ref := range guessedRefs {
		if guessConfidence[ref] > 10 {
			continue
		}
		ids = append(ids, ref)
	}
	contenders := make([]Contender, 0)
	var solutions [][]uint8
	for _, id := range ids {
		recoveredSnps = make(map[int]uint8)
		for l, locus := range p.Loci {
			if guess, ok := guessedSnps[locus]; ok {
				if recoveredPgsReferences[locus] != id {
					recoveredSnps[l] = guess
				}
			}
		}
		fmt.Printf("Solving w/o %s\n", id)
		// Decrease the confidence now so that it does not increase at the next step
		guessConfidence[id]--
		solutions = findSolutions(p, cohort, individual, indPop, recoveredSnps)
		if len(solutions) > 0 {
			contenders = append(contenders, Contender{id, solutions[0],
				solver.CalculateFullSequenceLikelihood(solutions[0], p.PopulationStats[indPop].AF, p.EffectAlleles)})
			fmt.Printf("Got %s: %.3f\n", id, contenders[len(contenders)-1].likelihood)
		} else {
			fmt.Printf("No solution w/o %s\n", id)
		}
	}
	sort.Slice(contenders, func(i, j int) bool {
		return contenders[i].likelihood < contenders[j].likelihood
	})
	// Resolving the references to see if there is a valid solution with the new SNP results
	for _, cnt := range contenders {
		pcnt := pgs.NewPGS()
		err = pcnt.LoadCatalogFile(path.Join(params.DataFolder, cnt.id+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return nil
		}
		err = pcnt.LoadStats()
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return nil
		}
		fmt.Printf("Re-solving %s with ", cnt.id)
		recoveredSnps = make(map[int]uint8)
		for l, prevLocus := range pcnt.Loci {
			for r, resolvedLocus := range p.Loci {
				if resolvedLocus != prevLocus {
					continue
				}
				recoveredSnps[l] = cnt.solution[pgs.Ploidy*r] + cnt.solution[pgs.Ploidy*r+1]
				break
			}
			fmt.Printf("%d:%d ", l, recoveredSnps[l])
		}
		fmt.Println()
		cohortcnt := solver.NewCohort(pcnt)
		solcnt := findSolutions(pcnt, cohortcnt, individual, indPop, recoveredSnps)
		if len(solcnt) > 0 {
			for i, lcnt := range pcnt.Loci {
				guessedSnps[lcnt] = solcnt[0][pgs.Ploidy*i] + solcnt[0][pgs.Ploidy*i+1]
			}
			fmt.Printf("New %s accuracy: %.3f\n", cnt.id, solver.Accuracy(solcnt[0], cohortcnt[individual].Genotype))
			solutions = [][]uint8{cnt.solution}
			guessConfidence[cnt.id] = 0
			break
		}
	}
	if len(solutions) == 0 {
		// Solve from scratch
		fmt.Println("No solutions with any subset, calculating from scratch")
		solutions = findSolutions(p, cohort, individual, indPop, make(map[int]uint8))
	}
	return solutions
}

func findSolutions(p *pgs.PGS, cht solver.Cohort, idv, ppl string, priorSnps map[int]uint8) [][]uint8 {
	var sm map[string][]uint8
	slv := solver.NewDP(cht[idv].Score, p, ppl, priorSnps)
	if len(p.Loci) < DeterminismLimit {
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
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("======== %s ========\n", p.PgsID)
	p.LoadStats()
	cohort := solver.NewCohort(p)
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
	var solmap map[string][]uint8
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

		slv := solver.NewDP(cohort[individual].Score, p, indPop, make(map[int]uint8))
		//
		if len(p.Loci) < DeterminismLimit {
			solmap = slv.SolveDeterministic(solver.UseLikelihood)
		} else {
			solmap = slv.SolveProbabilistic(solver.UseLikelihood)
		}
		solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
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
	INDIVIDUAL := "HG00124"
	//
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
	//catalogFile := "PGS000648_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh37.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
	catalogFile := "PGS002302_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000307_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000011_hmPOS_GRCh37.txt"
	//catalogFile := "PGS003436_hmPOS_GRCh37.txt"
	//catalogFile := "PGS002264_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004238_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004222_hmPOS_GRCh37.txt"
	//catalogFile := "PGS004249_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000589_hmPOS_GRCh37.txt"
	//catalogFile := "PGS001828_hmPOS_GRCh37.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	err = p.LoadStats()
	if err != nil {
		log.Printf("Error loading stats: %v\n", err)
		return
	}
	cohort := solver.NewCohort(p)

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
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
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
//	individuals := make([]string, 0, samplesPerPopulation*len(pgs.POPULATIONS))
//	var Count int
//	for _, ppl := range pgs.POPULATIONS {
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

//func sortingChoice() {
//	DeterminismLimit := 35
//	resultFolder := "results/sorting"
//	output := make([]*sortingOutput, 0)
//	//catalogFile := "PGS003760_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000753_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000154_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS003760_hmPOS_GRCh37.txt"
//	catalogFile := "PGS000851_hmPOS_GRCh37.txt"
//	samplesPerPopulation := 10
//	p := pgs.NewPGS()
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
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
//	individuals := make([]string, 0, samplesPerPopulation*len(pgs.POPULATIONS))
//	var Count int
//	for _, ppl := range pgs.POPULATIONS {
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
//
//	fmt.Println(individuals)
//	var solmap map[string][]uint8
//	for _, individual := range individuals {
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
//	outputPath := path.Join(resultFolder, fmt.Sprintf("%s.json", p.PgsID))
//	resFile, err := os.OpenFile(outputPath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	encoder := json.NewEncoder(resFile)
//	if err = encoder.Encode(output); err != nil {
//		log.Fatal("Cannot encode json", err)
//	}
//}

func accuracyLikelihood() {
	resFolder := "results/accuracyLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
	catalogFile := "PGS000037_hmPOS_GRCh37.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh37.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	populations := tools.LoadAncestry()
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	//samples := cohort.SortByScore()[len(cohort)-100:]
	samples := solver.All1000GenomesSamples()
	// Create a result file
	chunkNum, chunkSize := getChunkInfo(len(samples))
	//chunkNum := 2
	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.json", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	var acc float32
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d\n", i)
		indPop := populations[samples[i]]
		if strings.Contains(indPop, ",") {
			indPop = strings.Split(indPop, ",")[0]
		}
		slv := solver.NewDP(cohort[samples[i]].Score, p, indPop, make(map[int]uint8))
		solmap := slv.SolveDeterministic(solver.UseLikelihood)
		solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
		result := NewResult(samples[i], cohort[samples[i]].Score)
		for _, solution := range solutions {
			acc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
			result.Accuracies = append(result.Accuracies, fmt.Sprintf("%.3f", acc))
			result.Likelihoods = append(result.Likelihoods, fmt.Sprintf("%.3f",
				solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles)))
			//result.Spectrums = append(result.Spectrums, fmt.Sprintf("%.3f",
			//	solver.CalculateSolutionSpectrumDistance(solution, p.PopulationStats[indPop], p.EffectAlleles)))
		}
		if err = encoder.Encode(result); err != nil {
			log.Fatal("Cannot encode json", err)
		}
	}
	fmt.Println("Finished")
}

//func scoreToLikelihood() {
//	resFolder := "results/scoreLikelihood"
//	p := pgs.NewPGS()
//	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	populations := tools.LoadAncestry()
//	fmt.Printf("%s\n", p.PgsID)
//	cohort := solver.NewCohort(p)
//	samples := solver.All1000GenomesSamples()
//	// Create csv result file
//	chunkNum, chunkSize := getChunkInfo(len(samples))
//	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum))
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	writer := csv.NewWriter(resFile)
//	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
//		if i >= len(samples) {
//			break
//		}
//		fmt.Printf("%d ", i)
//		output := []string{samples[i], fmt.Sprintf("%.3f", cohort[samples[i]].Score)}
//		indPop := populations[samples[i]]
//		slv := solver.NewDP(cohort[samples[i]].Score, p, indPop, make(map[int]uint8))
//		solmap := slv.SolveFromSavedProbabilistic(solver.UseLikelihood)
//		solutions := solver.SortByAccuracy(solmap, cohort[samples[i]].Genotype)
//		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[samples[i]].Genotype) != 1.0 {
//			fmt.Printf("No solution for %s: %s, %.3f\n", samples[i], solver.ArrayToString(solutions[0]),
//				solver.Accuracy(solutions[0], cohort[samples[i]].Genotype))
//			continue
//		}
//		for _, solution := range solutions {
//			output = append(output, fmt.Sprintf("%.3f",
//				solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles)))
//		}
//		writer.Write(output)
//		writer.Flush()
//	}
//}

//func scoreToLikelihoodDistribution() {
//	resFolder := "results/scoreLikelihood"
//	p := pgs.NewPGS()
//	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	//fmt.Printf("%s\n", p.PgsID)
//	cohort := solver.NewCohort(p)
//	sorted := cohort.SortByScore()
//	numSamples := 5 * 2
//	step := (cohort[sorted[len(sorted)-1]].Score - cohort[sorted[0]].Score) / (float64(numSamples) * 10 / 6)
//	samples := make([]string, 0, numSamples)
//	for i := 0; i < len(sorted); i++ {
//		if len(samples) == numSamples {
//			break
//		}
//		if cohort[sorted[i]].Score >= cohort[sorted[0]].Score+float64(len(samples))*step {
//			samples = append(samples, sorted[i])
//			score := cohort[sorted[i]].Score
//			for {
//				i++
//				if i >= len(sorted) {
//					break
//				}
//				if cohort[sorted[i]].Score != score {
//					samples = append(samples, sorted[i])
//					break
//				}
//			}
//		}
//	}
//	for _, sample := range samples {
//		fmt.Printf("%s(%.3f) ", sample, cohort[sample].Score)
//	}
//	// Create csv result file
//	filepath := path.Join(resFolder, fmt.Sprintf("%s-distro.csv", p.PgsID))
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	writer := csv.NewWriter(resFile)
//	defer writer.Flush()
//	ctx, cancel := context.WithCancel(context.Background())
//	numThreads := getNumThreads()
//	for _, sample := range samples {
//		fmt.Printf("%s\n", sample)
//		output := []string{sample, fmt.Sprintf("%.3f", cohort[sample].Score)}
//		slv := solver.NewTwoSplitDP(ctx, cohort[sample].Score, p, numThreads)
//		solmap := slv.SolveFromSavedProbabilistic(numThreads)
//		solutions := solver.SortByAccuracy(solmap, cohort[sample].Genotype)
//		if len(solutions) == 0 || solver.CountAccuracy(solutions[0], cohort[sample].Genotype) != 1.0 {
//			fmt.Printf("No solution for %s: %s, %.3f\n", sample, solver.ArrayToString(solutions[0]),
//				solver.CountAccuracy(solutions[0], cohort[sample].Genotype))
//			continue
//		}
//		for _, solution := range solutions {
//			output = append(output, fmt.Sprintf("%.3f", p.CalculateFullSequenceLikelihood(solution)))
//		}
//		writer.Write(output)
//	}
//	cancel()
//}

//func likelihoodEffect() {
//	resFolder := "results"
//	p := pgs.NewPGS()
//	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	fmt.Printf("%s\n", p.PgsID)
//	cohort := solver.NewCohort(p)
//	samples := All1000GenomesSamples()
//	chunkNum, chunkSize := getChunkInfo(len(samples))
//	// Create csv result file
//	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum))
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	writer := csv.NewWriter(resFile)
//	defer writer.Flush()
//	writer.Write([]string{"individual",
//		"score", "major score", "median score",
//		"first accuracy", "major accuracy",
//		"true position", "total solutions",
//		"first likelihood", "true likelihood", "major likelihood"})
//	majorReference := p.AllReferenceAlleleSample()
//	majorScore := solver.CalculateDecimalSum(p.Context, majorReference, p.Weights)
//	majorLikelihood := p.CalculateFullSequenceLikelihood(majorReference)
//	medianScore := cohort[cohort.SortByScore()[len(cohort)/2]].Score
//	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
//		if i >= len(samples) {
//			break
//		}
//		fmt.Printf("%d ", i)
//		slv := solver.NewDP(cohort[samples[i]].Score, p, nil)
//		solmap := slv.SolveFromSavedProbabilistic()
//		solutions := solver.SortByLikelihood(solmap, p)
//		firstAcc := 0.0
//		firstLikelihood := 0.0
//		acc := 0.0
//	solLoop:
//		for j, solution := range solutions {
//			if j == 0 {
//				firstAcc = solver.CountAccuracy(solution, cohort[samples[i]].Genotype)
//				firstLikelihood = solver.CalculateFullSequenceLikelihood(solution)
//			}
//			acc = solver.CountAccuracy(solution, cohort[samples[i]].Genotype)
//			if acc == 1.0 {
//				writer.Write([]string{samples[i], fmt.Sprintf("%.3f", solver.CalculateDecimalSum(p.Context, solution, p.Weights)),
//					fmt.Sprintf("%.3f", majorScore), fmt.Sprintf("%.3f", medianScore),
//					fmt.Sprintf("%.3f", firstAcc),
//					fmt.Sprintf("%.3f", solver.CountAccuracy(majorReference, cohort[samples[i]].Genotype)),
//					fmt.Sprintf("%d", j), fmt.Sprintf("%d", len(solutions)), fmt.Sprintf("%.3f", firstLikelihood),
//					fmt.Sprintf("%.3f", solver.CalculateFullSequenceLikelihood(solution)), fmt.Sprintf("%.3f", majorLikelihood)})
//				break solLoop
//			}
//		}
//		if acc != 1.0 {
//			diff := new(apd.Decimal)
//			p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalSum(p.Context, cohort[samples[i]].Genotype, p.Weights))
//			fmt.Printf("\nNo solution for %s\n", samples[i])
//			fmt.Printf("\nTrue:\n%s -- %s, %f\n", solver.ArrayToString(cohort[samples[i]].Genotype),
//				diff.String(),
//				solver.CalculateFullSequenceLikelihood(cohort[samples[i]].Genotype))
//			for _, sol := range solutions {
//				p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalSum(p.Context, sol, p.Weights))
//				fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(sol), solver.CountAccuracy(sol, cohort[samples[i]].Genotype),
//					diff.String(), solver.CalculateFullSequenceLikelihood(sol))
//			}
//		}
//	}
//}

func getChunkInfo(totalLen int) (int, int) {
	var err error
	var chunkNum, chunkSize = 0, 0
	if len(os.Args) > 2 {
		chunkNum, err = strconv.Atoi(os.Args[1])
		if err != nil {
			log.Fatalf("Error parsing chunkNum %s: %v", os.Args[1], err)
		}
		chunkSize, err = strconv.Atoi(os.Args[2])
		if err != nil {
			log.Fatalf("Error parsing chunkSize %s: %v", os.Args[2], err)
		}
	} else {
		chunkSize = totalLen
	}
	return chunkNum, chunkSize
}

func likelihoodWeight() {
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000043_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh37.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
	catalogFile := "PGS002302_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000307_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000534_hmPOS_GRCh37.txt"
	//catalogFile := "PGS000011_hmPOS_GRCh37.txt"
	//catalogFile := "PGS003436_hmPOS_GRCh37.txt"
	//catalogFile := "PGS002264_hmPOS_GRCh37.txt"
	//catalogFile := "PGS003760_hmPOS_GRCh37.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	fmt.Printf("%s\n", p.PgsID)
	type Relation struct {
		Weights []float64
		AF      map[string][]float32
	}
	r := new(Relation)
	r.Weights = make([]float64, len(p.Weights))
	for i := 0; i < len(p.Weights); i++ {
		r.Weights[i], err = p.Weights[i].Float64()
		if err != nil {
			log.Fatalf("Error converting weight %d: %v\n", i, err)
		}
	}
	r.AF = make(map[string][]float32)
	for _, pop := range pgs.POPULATIONS {
		r.AF[pop] = make([]float32, len(p.PopulationStats[pop].AF))
		for i := 0; i < len(p.PopulationStats[pop].AF); i++ {
			r.AF[pop][i] = p.PopulationStats[pop].AF[i][p.EffectAlleles[i]]
		}
	}

	resFolder := "results/weights"
	filepath := path.Join(resFolder, fmt.Sprintf("%s.json", p.PgsID))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(r); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

//func scoreDistribution() {
//	p := pgs.NewPGS()
//	catalogFile := "PGS000073_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000040_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000648_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000891_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS001827_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh37.txt"
//	//catalogFile := "PGS000066_hmPOS_GRCh37.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	cohort := solver.NewCohort(p)
//	samples := All1000GenomesSamples()
//	indices := []int{0, len(p.Weights) / 2, len(p.Weights)}
//
//	resFolder := "results"
//	filepath := path.Join(resFolder, "distro.json")
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	encoder := json.NewEncoder(resFile)
//
//	multiplier := apd.New(1, int32(p.WeightPrecision))
//	if p.WeightPrecision > params.PrecisionsLimit {
//		multiplier.SetFinite(1, params.PrecisionsLimit)
//	}
//	weights := solver.DecimalsToInts(p.Context, p.Weights, multiplier)
//	maxTotalPositive, maxTotalNegative := solver.GetMaxTotal(weights)
//	mf, _ := multiplier.Float64()
//	maxTotalPositiveF, maxTotalNegativeF := float64(maxTotalPositive)/mf, float64(maxTotalNegative)/mf
//	if err = encoder.Encode([]float64{maxTotalPositiveF, maxTotalNegativeF}); err != nil {
//		log.Fatal("Cannot encode json for major-minor", err)
//	}
//
//	betas := make(map[uint16]int64, len(p.Weights))
//	for i := 0; i < len(weights); i++ {
//		betas[uint16(i)] = weights[i]
//	}
//	sampledMax, sampledMin := solver.SampleMaxMinScores(0, len(p.Weights), 100*len(p.Weights), betas, p)
//	sampledMaxF, sampledMinF := float64(sampledMax)/mf, float64(sampledMin)/mf
//
//	if err = encoder.Encode([]float64{sampledMaxF, sampledMinF}); err != nil {
//		log.Fatal("Cannot encode json for major-minor", err)
//	}
//
//	fmt.Println(sampledMaxF, sampledMinF)
//	scores := make([][]float64, 2)
//	for r := 0; r < len(scores); r++ {
//		scores[r] = make([]float64, 0, len(samples))
//	}
//	for _, sample := range samples {
//		for i := 0; i < len(indices)-1; i++ {
//			score, _ := solver.CalculateDecimalSum(p.Context, cohort[sample].Genotype[indices[i]*pgs.Ploidy:indices[i+1]*pgs.Ploidy],
//				p.Weights[indices[i]:indices[i+1]]).Float64()
//			scores[i] = append(scores[i], score)
//		}
//		//fullscore, _ := solver.CalculateDecimalSum(p.Context, cohort[sample].Genotype, p.Weights).Float64()
//		//fmt.Println(sample, fullscore, scores[0][len(scores[0])-1], scores[1][len(scores[1])-1])
//	}
//	for r := 0; r < len(scores); r++ {
//		if err = encoder.Encode(scores[r]); err != nil {
//			log.Fatal("Cannot encode scores json", err)
//		}
//	}
//}

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

func allScores(catalogId string) []float64 {
	filename := path.Join(params.DataFolder, fmt.Sprintf("%s.scores", catalogId))
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Cannot open file %s: %v", filename, err)
	}
	defer file.Close()
	scores := make([]float64, 0, 2504)
	var score float64
	for {
		_, err := fmt.Fscanf(file, "%f\n", &score)
		if err != nil {
			break
		}
		scores = append(scores, score)
	}
	return scores
}

func numSNPs(seq []uint8) int {
	num := 0
	for _, val := range seq {
		if val != 0 {
			num++
		}
	}
	return num
}

func numHeterozygous(seq []uint8) int {
	num := 0
	for i := 0; i < len(seq); i += pgs.Ploidy {
		if seq[i] == 1 && seq[i+1] == 1 {
			num++
		}
	}
	return num
}

//func selectPGS() {
//	var pgsToNumVariants map[string]int
//	var pgsToPgp map[string]string
//	numVariantsFile, err := os.Open("results/filtered_pgs.json")
//	if err != nil {
//		log.Fatalf("Error opening file: %v", err)
//	}
//	defer numVariantsFile.Close()
//	decoder := json.NewDecoder(numVariantsFile)
//	err = decoder.Decode(&pgsToNumVariants)
//	if err != nil {
//		log.Fatalf("Error decoding num variants json: %v", err)
//	}
//	pgpFile, err := os.Open("results/pgs_pgp.json")
//	if err != nil {
//		log.Fatalf("Error opening file: %v", err)
//	}
//	defer pgpFile.Close()
//	decoder = json.NewDecoder(pgpFile)
//	err = decoder.Decode(&pgsToPgp)
//	if err != nil {
//		log.Fatalf("Error decoding pgp json: %v", err)
//	}
//
//	baggedPgp := make(map[string]struct{})
//	pgsToSelect := make(map[string]int)
//
//	for i := 5; i <= 50; i += 5 {
//		if i == 50 {
//			i = 49
//		}
//		for pgs, num := range pgsToNumVariants {
//			if num == i {
//				if _, ok := baggedPgp[pgsToPgp[pgs]]; ok {
//					continue
//				}
//				pgsToSelect[pgs] = num
//				baggedPgp[pgsToPgp[pgs]] = struct{}{}
//				break
//			}
//		}
//	}
//	fmt.Println(pgsToSelect)
//}
