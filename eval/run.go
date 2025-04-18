package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/data"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"gonum.org/v1/gonum/stat/distuv"
)

const (
	DeterminismLimit       = 34
	ScratchSolvingSnpLimit = 44
	ConfidenceThreshold    = 5
	perAncestryAf          = "ancestry"
	globalAf               = "global"
)

func main() {
	expr := flag.String("e", "", "Experiment type")
	flag.Parse()
	switch *expr {
	case "demo":
		demo()
	case "solve":
		chainSolving()
	case "solveaf":
		chainSolvingOneAF()
	case "accuracy":
		recoveryAccuracy(perAncestryAf)
	case "gaccuracy":
		recoveryAccuracy(globalAf)
	case "linking":
		linking()
	case "uniqueness_gg":
		uniquenessExperiment(data.GG)
	case "uniqueness_uk":
		uniquenessExperiment(data.UKB)
	case "scores":
		calculateScoresAndStats()
	case "rounding":
		precisionToUniqueness()
	case "snps":
		numberSnpsToPossibleScores()
	case "percentiles":
		percentilesToScores()
	case "genfreq":
		calculateGenotypeFrequenciesForGuessed()
	}
}

func demo() {
	idsToNumVariants, lociToPgs, err := loadValidatedPgsAndLoci()
	if err != nil {
		log.Println(err)
		return
	}
	individuals := solver.All1000GenomesAndRelativeSamples()
	individual := individuals[rand.Intn(len(individuals))]
	populations := data.LoadAncestry()
	idvPop := pgs.GetIndividualAncestry(individual, populations)
	fmt.Printf("--- Recovering the genotypes of %s (%s) from the 1000 Genomes dataset ---\n", individual, idvPop)

	var pgsID string
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
	guessedSnps := make(map[string]uint8)
	guessedRefs := make(map[string]string)
	guessConfidence := make(map[string]int)
	trueSnps := make(map[string]uint8)
	for {
		if len(allPgs) == 0 || idsToNumVariants[allPgs[0]] > 30 {
			break
		}
		pgsID = allPgs[0]
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		fmt.Printf("======== %s ========\n", p.PgsID)
		err = p.LoadStats(data.GG)
		if err != nil {
			log.Printf("Error loading stats: %v\n", err)
			return
		}
		knownSnps := make(map[int]uint8)
		knownRefs := make(map[int]string)
		includedRefs := make(map[string]struct{})
		for l, locus := range p.Loci {
			if guess, ok := guessedSnps[locus]; ok {
				knownSnps[l] = guess
				knownRefs[l] = guessedRefs[locus]
				if _, ok = includedRefs[guessedRefs[locus]]; !ok {
					includedRefs[guessedRefs[locus]] = struct{}{}
				}
			}
		}
		fmt.Printf("Total SNPs %d, unknown %d\n", len(p.Loci), len(p.Loci)-len(knownSnps))
		cohort := solver.NewCohort(p, data.GG)
		solutions := findSolutions(p, cohort, individual, idvPop, knownSnps)
		if len(solutions) == 0 && len(knownSnps) > 0 {
			solutions = selfRepair(p, cohort, individual, idvPop, guessedSnps, guessedRefs, guessConfidence,
				knownSnps, knownRefs)
			knownSnps = nil
		}
		if len(solutions) > 0 {
			for i, locus := range p.Loci {
				if _, ok := guessedSnps[locus]; ok {
					continue
				}
				guessedSnps[locus] = solutions[0][pgs.Ploidy*i] + solutions[0][pgs.Ploidy*i+1]
				guessedRefs[locus] = pgsID
				trueSnps[locus] = cohort[individual].Genotype[pgs.Ploidy*i] + cohort[individual].Genotype[pgs.Ploidy*i+1]
			}
			fmt.Printf("Recovered genotypes:\t%s\n", solver.ArrayToString(solutions[0]))
			fmt.Printf("True genotypes:\t\t%s\n", solver.ArrayToString(cohort[individual].Genotype))
			fmt.Printf("Accuracy: %.0f%%\n", 100*solver.Accuracy(cohort[individual].Genotype, solutions[0]))
			//
			guessConfidence[pgsID] = 1
			for ref := range includedRefs {
				guessConfidence[ref]++
				if len(p.Loci)-len(knownSnps) <= 1 {
					guessConfidence[ref]++
				}
			}
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
	accuracy := 0.0
	for locus, guess := range guessedSnps {
		if guess == trueSnps[locus] {
			accuracy++
		}
	}
	accuracy /= float64(len(guessedSnps))
	fmt.Println("==============================")
	fmt.Printf("Recovered %d genotypes for %s (%s) with %.1f%% accuracy\n",
		len(guessedSnps), individual, idvPop, accuracy*100)
}

func numberSnpsToPossibleScores() {
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
	type Result struct {
		NumVariants         int
		TotalPossibleScores uint64
		Precision           uint32
	}
	results := make([]Result, 0)
	for _, pgsID := range allPgs {
		fmt.Printf("======== %s (%d) ========\n", pgsID, pgsToNumVariants[pgsID])
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		results = append(results, Result{
			NumVariants:         len(p.Loci),
			TotalPossibleScores: calculateNumPossibleScores(p),
			Precision:           p.WeightPrecision,
		})
	}
	resultFolder := "results/uniqueness"
	filepath := path.Join(resultFolder, "num_snps.json")
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

func precisionToUniqueness() {
	//file, _ := os.Open("results/validated_pgs.json")
	//defer file.Close()
	//decoder := json.NewDecoder(file)
	//var pgsToNumVariants map[string]int
	//decoder.Decode(&pgsToNumVariants)
	//allPgs := make([]string, 0, len(pgsToNumVariants))
	//for id, _ := range pgsToNumVariants {
	//	allPgs = append(allPgs, id)
	//}
	//sort.Slice(allPgs, func(i, j int) bool {
	//	if pgsToNumVariants[allPgs[i]] == pgsToNumVariants[allPgs[j]] {
	//		return allPgs[i] < allPgs[j]
	//	}
	//	return pgsToNumVariants[allPgs[i]] < pgsToNumVariants[allPgs[j]]
	//})
	//for _, pgsID := range allPgs {
	//	p := pgs.NewPGS()
	//	p.LoadCatalogFile(path.Join(params.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
	//	if p.WeightPrecision >= 17 {
	//		maxw := p.FindMaxAbsoluteWeight() * math.Pow(10, float64(p.WeightPrecision))
	//		density := float64(p.NumVariants) / Log3(maxw)
	//		fmt.Printf("%s (precision %d, variants %d, density %.3f)\n", pgsID, p.WeightPrecision,
	//			pgsToNumVariants[pgsID], density)
	//	}
	//}
	//return
	type Result struct {
		Precision              uint32
		Density                float64
		TotalPossibleScores    uint64
		PresentScores          int
		PercentageUnique       float32
		MedianAnonymitySetSize int
	}
	dataset := data.UKB
	pgsID := "PGS000869"
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats(dataset)
	exp := apd.New(1, int32(p.WeightPrecision))
	weights := solver.ScaleWeights(p.Context, p.Weights, exp)
	for i := range p.Weights {
		fmt.Printf("%s\t->\t%s\n", p.Weights[i].String(), weights[i].String())
	}
	cohort := solver.NewCohort(p, dataset)
	results := make([]Result, 0)
	ten := new(apd.BigInt).SetInt64(10)
	divisor := new(apd.Decimal).SetInt64(int64(len(weights) * pgs.Ploidy))
	individuals := solver.AllUKBiobankSamples()
	allScores := make(map[uint32][]string)
	for wp := p.WeightPrecision; wp > 0; wp-- {
		fmt.Printf("=== %d ===\n", wp)
		allScores[wp] = make([]string, 0)
		sums := make(map[string]int)
		preScores := make([]string, 0)
		var sf string
		for _, idv := range individuals {
			sf = solver.CalculateBigIntSum(cohort[idv].Genotype, weights, p.EffectAlleles).String()
			preScores = append(preScores, sf)
			sums[sf]++
		}
		pow := apd.New(1, -int32(wp))
		fmt.Printf("Pow: %s\n", pow.String())
		for _, score := range preScores {
			sbi := new(apd.Decimal)
			_, _, err := sbi.SetString(score)
			if err != nil {
				log.Fatalf("Error setting string: %v", err)
			}
			p.Context.Mul(sbi, sbi, pow)
			p.Context.Quo(sbi, sbi, divisor)
			allScores[wp] = append(allScores[wp], sbi.String())
		}
		realNumUnique := 0
		anonsets := make([]int, 0)
		for _, count := range sums {
			if count == 1 {
				realNumUnique++
			}
			anonsets = append(anonsets, count)
		}
		density := float64(p.NumVariants) / Log3(float64(solver.FindMaxAbsoluteBigInt(weights)))
		sort.Ints(anonsets)
		results = append(results, Result{
			Precision:              wp,
			Density:                density,
			TotalPossibleScores:    getTotalNumberPossibleScores(weights),
			PresentScores:          len(sums),
			PercentageUnique:       float32(realNumUnique) * 100 / float32(len(cohort)),
			MedianAnonymitySetSize: anonsets[len(anonsets)/2],
		})
		for i := range weights {
			fmt.Printf("%s ", weights[i].String())
		}
		fmt.Printf("\nPrecision: %d, Density: %.3f, Total possible sums: %d, Present sums: %d, Unique: %.2f%%, Median set size: %d\n",
			wp, density, results[len(results)-1].TotalPossibleScores, len(sums),
			float32(realNumUnique)*100/float32(len(cohort)), anonsets[len(anonsets)/2])
		for i := range weights {
			weights[i] = divideAndRound(weights[i], ten)
		}
	}
	resultFolder := "results/rounding"
	filepath := path.Join(resultFolder, fmt.Sprintf("%s_stats.json", pgsID))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	filepath = path.Join(resultFolder, fmt.Sprintf("%s_scores.json", pgsID))
	distroFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer distroFile.Close()
	encoder = json.NewEncoder(distroFile)
	if err = encoder.Encode(allScores); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func percentilesToScores() {
	//type Result struct {
	//	Precision              uint32
	//	Density                float64
	//	TotalPossibleScores    uint64
	//	PresentScores          int
	//	PercentageUnique       float32
	//	MedianAnonymitySetSize int
	//}
	dataset := data.UKB
	pgsID := "PGS000869"
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats(dataset)
	cohort := solver.NewCohort(p, dataset)
	//results := make([]Result, 0)
	individuals := solver.AllUKBiobankSamples()
	scores := make([]*apd.Decimal, 0, len(individuals))
	for _, idv := range individuals {
		scores = append(scores, cohort[idv].Score)
	}
	// Sort scores in the same order as individuals
	sort.SliceStable(scores, func(i, j int) bool {
		return scores[i].Cmp(scores[j]) < 0
	})

	idvPerPercentile := len(scores) / 100
	scoresPerPercentile := make(map[int]int, 100)
	var start, end int
	var uniqueScoresInPercentile map[string]bool
	for percentile := 1; percentile <= 100; percentile++ {
		start = idvPerPercentile * (percentile - 1)
		end = idvPerPercentile * percentile
		if percentile == 100 {
			end = len(scores)
		}
		uniqueScoresInPercentile = make(map[string]bool)
		for i := start; i < end; i++ {
			uniqueScoresInPercentile[scores[i].String()] = true
		}
		scoresPerPercentile[percentile] = len(uniqueScoresInPercentile)
	}

	for percentile := 1; percentile <= 100; percentile++ {
		fmt.Printf("Percentile %d: %d\n", percentile, scoresPerPercentile[percentile])
	}

	//resultFolder := "results/percentiles"
	//filepath := path.Join(resultFolder, fmt.Sprintf("%s.json", pgsID))
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

func divideAndRound(n, divisor *apd.BigInt) *apd.BigInt {
	quotient := new(apd.BigInt).Quo(n, divisor)
	remainder := new(apd.BigInt).Rem(n, divisor)

	if remainder.Cmp(apd.NewBigInt(5)) >= 0 {
		quotient.Add(quotient, apd.NewBigInt(1))
	}

	return quotient
}

func getTotalNumberPossibleScores(weights []*apd.BigInt) uint64 {
	numPossibleSubsets := calculateNumPossibleSubsets(weights)
	minScoreBig, maxScoreBig, secondMinBig, secondMaxBig := FindMinAndMaxBigIntSums(weights)
	minVal := minScoreBig.Int64()
	maxVal := maxScoreBig.Int64()
	secondMin := secondMinBig.Int64()
	secondMax := secondMaxBig.Int64()
	fmt.Printf("Min: %d, 2ndMin: %d, Max: %d, 2ndMax: %d\n", minVal, secondMin,
		maxVal, secondMax)
	// Adding 2 to account for the minimum and maximum scores
	numPossibleScoresDuePrecision := uint64(secondMax - secondMin + 1 + 2)
	numPossibleScores := numPossibleScoresDuePrecision
	if numPossibleSubsets > 0 && numPossibleSubsets < numPossibleScoresDuePrecision {
		fmt.Printf("Number of subsets: %d < In-between scores: %d\n", numPossibleSubsets, numPossibleScoresDuePrecision)
		numPossibleScores = numPossibleSubsets
	} else {
		fmt.Printf("Num subsums: %d > In-between scores: %d\n", numPossibleSubsets, numPossibleScoresDuePrecision)
	}
	return numPossibleScores
}

func FindMinAndMaxBigIntSums(weights []*apd.BigInt) (*apd.BigInt, *apd.BigInt, *apd.BigInt, *apd.BigInt) {
	maxScore, minScore := apd.NewBigInt(0), apd.NewBigInt(0)
	secondMaxScore, secondMinScore := apd.NewBigInt(0), apd.NewBigInt(0)
	var smallestPositiveWeight, smallestNegativeWeight *apd.BigInt
	allPositive, allNegative := true, true
	for _, weight := range weights {
		if weight.Sign() == -1 {
			minScore.Add(minScore, weight)
			minScore.Add(minScore, weight)
			allPositive = false
			if smallestNegativeWeight == nil || weight.Cmp(smallestNegativeWeight) > 0 {
				smallestNegativeWeight = weight
			}
			continue
		}
		maxScore.Add(maxScore, weight)
		maxScore.Add(maxScore, weight)
		allNegative = false
		if smallestPositiveWeight == nil || weight.Cmp(smallestPositiveWeight) < 0 {
			smallestPositiveWeight = weight
		}
	}
	switch {
	case allNegative:
		secondMaxScore.Add(maxScore, smallestNegativeWeight)
		secondMinScore.Sub(minScore, smallestNegativeWeight)
	case allPositive:
		secondMinScore.Set(smallestPositiveWeight)
		secondMaxScore.Sub(maxScore, smallestPositiveWeight)
	default:
		secondMaxScore.Sub(maxScore, smallestPositiveWeight)
		secondMinScore.Sub(minScore, smallestNegativeWeight)
	}
	return minScore, maxScore, secondMinScore, secondMaxScore
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
	err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, allPgs[pgsNum]+"_hmPOS_GRCh37.txt"))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats(data.UKB)
	solver.NewCohort(p, data.UKB)
	fmt.Printf("%s complete\n", allPgs[pgsNum])
}

func uniquenessExperiment(dataset string) {
	fmt.Printf("--- Uniqueness experiment ---\n")
	type Result struct {
		PgsID                       string
		NumVariants                 int
		TotalPresentScores          int
		TotalPossibleScores         uint64
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
		if pgsToNumVariants[allPgs[i]] == pgsToNumVariants[allPgs[j]] {
			return allPgs[i] < allPgs[j]
		}
		return pgsToNumVariants[allPgs[i]] < pgsToNumVariants[allPgs[j]]
	})

	chunkNum, totalChunks := 0, 0
	if len(os.Args) > 3 {
		chunkNum, err = strconv.Atoi(os.Args[2])
		if err != nil {
			log.Fatalf("Error parsing chunkNum %s: %v", os.Args[2], err)
		}
		totalChunks, err = strconv.Atoi(os.Args[3])
		if err != nil {
			log.Fatalf("Error parsing totalChunks %s: %v", os.Args[3], err)
		}
	}

	precisionEstimationLimit := 13
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
				var sf string
				var realNumUnique int
				p := pgs.NewPGS()
				err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
				if err != nil {
					log.Printf("Error loading catalog file: %v\n", err)
					return
				}
				p.LoadStats(dataset)
				fmt.Printf("======== %s (%d) ========\n", p.PgsID, len(p.Loci))
				cohort := solver.NewCohort(p, dataset)
				fmt.Printf("%s: cohort size %d\n", p.PgsID, len(cohort))
				scores := make(map[string]int)
				for idv := range cohort {
					sf = cohort[idv].Score.String()
					scores[sf]++
				}
				realNumUnique = 0
				anonsets := make([]int, 0)
				for _, count := range scores {
					if count == 1 {
						realNumUnique++
					}
					anonsets = append(anonsets, count)
				}
				predictedNumUnique, predictedMedianAnonSize, numPossibleScores :=
					estimateNumUniqueScores(p, len(cohort), precisionEstimationLimit)
				if math.IsNaN(predictedNumUnique) || predictedNumUnique == -1 {
					predictedNumUnique = -1
					fmt.Printf("Too many possible scores for %s\n", pgsID)
				}
				result := &Result{
					PgsID:                       pgsID,
					NumVariants:                 pgsToNumVariants[pgsID],
					TotalPresentScores:          len(scores),
					TotalPossibleScores:         numPossibleScores,
					RealPercentageUnique:        float64(realNumUnique) * 100 / float64(len(cohort)),
					PredictedPercentageUnique:   predictedNumUnique * 100 / float64(len(cohort)),
					RealMedianAnonymitySet:      median(anonsets),
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
	for i, pgsID := range allPgs {
		if i%totalChunks != chunkNum {
			continue
		}
		tasks <- pgsID
	}
	fmt.Println("---- All tasks sent ----")
	close(tasks)
	wg.Wait()
	for _, res := range results {
		fmt.Println(*res)
	}
	resultFolder := "results/uniqueness"
	filepath := path.Join(resultFolder, fmt.Sprintf("scores_%s_%d.json", dataset, chunkNum))
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

func calculateNumPossibleScores(p *pgs.PGS) uint64 {
	numPossibleSubsets := calculateNumPossibleSubsets(p.Weights)
	_, _, secondMinDecimal, secondMaxDecimal := p.FindMinAndMaxScores()
	secondMin, _ := secondMinDecimal.Float64()
	secondMax, _ := secondMaxDecimal.Float64()
	step := math.Pow10(-int(p.WeightPrecision))
	// Adding 2 to account for the minimum and maximum scores
	numPossiblePrecisionScores := (uint64((secondMax-secondMin)/step) + 1) + 2
	if numPossibleSubsets < numPossiblePrecisionScores {
		return numPossibleSubsets
	}
	return numPossiblePrecisionScores
}

func estimateNumUniqueScores(p *pgs.PGS, M int, powerUpperBound int) (float64, int, uint64) {
	numPossibleSubsets := calculateNumPossibleSubsets(p.Weights)
	minScoreDecimal, maxScoreDecimal, secondMinDecimal, secondMaxDecimal := p.FindMinAndMaxScores()
	minVal, _ := minScoreDecimal.Float64()
	maxVal, _ := maxScoreDecimal.Float64()
	secondMin, _ := secondMinDecimal.Float64()
	secondMax, _ := secondMaxDecimal.Float64()
	step := math.Pow10(-int(p.WeightPrecision))
	// Adding 2 to account for the minimum and maximum scores
	numPossibleScoresDuePrecision := (uint64((secondMax-secondMin)/step) + 1) + 2
	numPossibleScores := numPossibleScoresDuePrecision
	if numPossibleSubsets > 0 && numPossibleSubsets < numPossibleScoresDuePrecision {
		fmt.Printf("Number of subsets: %d < In-between scores: %d\n", numPossibleSubsets, numPossibleScoresDuePrecision)
		step = (secondMax - secondMin) / float64(numPossibleSubsets-2)
		numPossibleScores = numPossibleSubsets
	} else {
		fmt.Printf("Num subsums: %d > In-between scores: %d\n", numPossibleSubsets, numPossibleScoresDuePrecision)
	}
	if numPossibleScores > uint64(math.Pow10(powerUpperBound)) {
		return -1, -1, numPossibleScores
	}

	predictedMeans, predictedStds := p.EstimateMeanAndStd()
	normalDist := distuv.Normal{Mu: predictedMeans["ALL"], Sigma: predictedStds["ALL"]}

	type result struct {
		unique float64
		sets   map[uint16]uint16
	}
	numWorkers := 8
	segmentSize := (secondMax - secondMin) / float64(numWorkers)
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
			for x := secondMin + float64(i)*segmentSize; x < minVal+float64(i+1)*segmentSize; x += step {
				lowerBound = x - step/2
				upperBound = x + step/2
				prob = normalDist.CDF(upperBound) - normalDist.CDF(lowerBound)
				averageCount = prob * float64(M)
				segmentTotalUnique += averageCount * math.Pow(1-prob, float64(M-1))
				roundedCount = probabilisticallyRound(averageCount)
				if roundedCount == 0 {
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
	// Calculating the probability of the minimum and maximum scores (we have excluded them before to skip
	// the impossible scores [min, 2ndMin] and [2ndMax, max])
	totalUnique += estimateScore(minVal, normalDist, M, step, allSets)
	totalUnique += estimateScore(maxVal, normalDist, M, step, allSets)
	//
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

func estimateScore(score float64, dist distuv.Normal, M int, step float64, sets map[uint16]uint16) float64 {
	prob := dist.CDF(score+step/2) - dist.CDF(score-step/2)
	averageCount := prob * float64(M)
	roundedCount := probabilisticallyRound(averageCount)
	if roundedCount > 0 {
		if _, exists := sets[roundedCount]; exists {
			sets[roundedCount]++
		} else {
			sets[roundedCount] = 1
		}
	}
	return averageCount * math.Pow(1-prob, float64(M-1))
}

func probabilisticallyRound(value float64) uint16 {
	rounded := uint16(math.Round(value))
	if value < 1 {
		randomValue := rand.Float64()
		if randomValue >= value {
			return 0
		}
		return 1
	}
	return rounded
}

func median(values []int) int {
	sort.Ints(values)
	n := len(values)
	if n%2 == 0 {
		return (values[n/2-1] + values[n/2]) / 2
	}
	return values[n/2]
}

type Stringer interface {
	String() string
}

func calculateNumPossibleSubsets[T Stringer](weights []T) uint64 {
	weightRepetitions := make(map[string]int)
	for i := range weights {
		w := weights[i].String()
		weightRepetitions[w]++
	}
	numPossibleSubsets := uint64(1)
	for _, count := range weightRepetitions {
		numPossibleSubsets *= 2*uint64(count) + 1
	}
	return numPossibleSubsets
}

func chainSolving() {
	fmt.Printf("Chain solving\n")
	idsToNumVariants, lociToPgs, err := loadValidatedPgsAndLoci()
	if err != nil {
		log.Println(err)
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

	type Guess struct {
		Individual string
		Ancestry   string
		SNPs       map[string]uint8
	}
	guesses := make([]*Guess, 0)
	populations := data.LoadAncestry()
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
		for {
			if len(allPgs) == 0 {
				break
			}
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats(data.GG)
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
			cohort := solver.NewCohort(p, data.GG)
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
				fmt.Println("!!! Still could not find a solution !!!")
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
			_ = p
		}
		guesses = append(guesses, &Guess{Individual: individual, Ancestry: indPop, SNPs: guessedSnps})
	}

	resultFolder := "results/recoveryOutput"
	filepath := path.Join(resultFolder, fmt.Sprintf("guesses%d.json", chunkNum))
	guessFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer guessFile.Close()
	encoder := json.NewEncoder(guessFile)
	if err = encoder.Encode(guesses); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func chainSolvingOneAF() {
	fmt.Printf("Chain solving\n")
	idsToNumVariants, lociToPgs, err := loadValidatedPgsAndLoci()
	if err != nil {
		log.Println(err)
		return
	}
	individuals := getIndividualsSample()
	ppl := os.Args[2]
	//fmt.Printf("%s: %s\n", ppl, individuals[ppl])

	type Guess struct {
		Individual string
		Ancestry   string
		SNPs       map[string]uint8
	}
	guesses := make([]*Guess, 0)
	populations := data.LoadAncestry()
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
		indPop := "ALL"
		guessedSnps := make(map[string]uint8)
		guessedRefs := make(map[string]string)
		guessConfidence := make(map[string]int)
		for {
			if len(allPgs) == 0 {
				break
			}
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file: %v\n", err)
				return
			}
			fmt.Printf("======== %s ========\n", p.PgsID)
			err = p.LoadStats(data.GG)
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
			cohort := solver.NewCohort(p, data.GG)
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
		guesses = append(guesses, &Guess{Individual: individual,
			Ancestry: pgs.GetIndividualAncestry(individual, populations), SNPs: guessedSnps})
	}

	resultFolder := "results/recoveryOutput"
	filepath := path.Join(resultFolder, fmt.Sprintf("af-guesses-%s.json", ppl))
	guessFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer guessFile.Close()
	encoder := json.NewEncoder(guessFile)
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
	file, err := os.Open("results/recoveryAccuracy/individuals.json")
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
	fmt.Println("--- No solutions with all the known genotypes ---")
	fmt.Println("Running self-repair")
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
	//fmt.Printf("Solving with high-confidence SNPs: %v\n", highConfidenceRefs)
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
	fmt.Printf("Top solution from scratch: accuracy %.2f\n",
		solver.Accuracy(solutions[0], cohort[individual].Genotype))

	refMap := make(map[string]struct{})
	refs := make([]string, 0)
	for _, ref := range recoveredRefs {
		if guessConfidence[ref] > ConfidenceThreshold && highConfidence {
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
			err = refp.LoadCatalogFile(path.Join(data.LocalInputFolder, ref+"_hmPOS_GRCh37.txt"))
			if err != nil {
				log.Printf("Error loading catalog file %s: %v\n", ref, err)
				return nil
			}
			err = refp.LoadStats(data.GG)
			if err != nil {
				log.Printf("Error loading stats for %s: %v\n", ref, err)
				return nil
			}
			refCohort := solver.NewCohort(refp, data.GG)
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
		// All refs are solvable, updating old loci
		for j, ref := range refs {
			if len(refSols[j]) == 0 {
				continue
			}
			refp := pgs.NewPGS()
			err = refp.LoadCatalogFile(path.Join(data.LocalInputFolder, ref+"_hmPOS_GRCh37.txt"))
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
		sm = slv.SolveDeterministic()
	} else {
		sm = slv.SolveProbabilistic()
	}
	return solver.SortByLikelihood(sm, p.PopulationStats[ppl], p.EffectAlleles)
}

func linking() {
	guessedSnps := loadGuessedGenotypes(perAncestryAf)
	guessedIndividuals := make([]string, 0)
	for idv := range guessedSnps {
		guessedIndividuals = append(guessedIndividuals, idv)
	}
	sort.Strings(guessedIndividuals)
	vcfFile := "results/linking/plink.vcf"
	compressedVcfFile := vcfFile + ".gz"
	preparePlinkInput(guessedIndividuals, guessedSnps, vcfFile)
	compressVcfFile(vcfFile, compressedVcfFile)
	indexVcfFile(compressedVcfFile)
	plinkGCTA(compressedVcfFile)
	plinkKing(compressedVcfFile)
	removeFiles(vcfFile, compressedVcfFile)
}

func preparePlinkInput(guessedIndividuals []string, guessed map[string]map[string]map[string]uint8, output string) {
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
		retrieved = getIndividualPositionStrings(chr, strPos, mainIndividuals, data.GG)
		transferSnps(trueSNPs, retrieved, chr)
		retrieved = getIndividualPositionStrings(chr, strPos, relativeIndividuals, data.RL)
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
	vcfFile, err := os.OpenFile(output, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Cannot create plink ibd file: %v", err)
	}
	defer vcfFile.Close()
	_, err = vcfFile.WriteString(formatHeader)
	if err != nil {
		log.Fatalf("Cannot write format header to the plink file: %v", err)
	}
	_, err = vcfFile.WriteString(samplesHeader)
	if err != nil {
		log.Fatalf("Cannot write samples header to the plink file: %v", err)
	}
	var snp uint8
	for chrId := 1; chrId <= 22; chrId++ {
		chr = strconv.Itoa(chrId)
	positionLoop:
		for _, posId := range positions[chr] {
			pos := strconv.Itoa(posId)
			if _, ok = alleles[chr]; !ok || alleles[chr][pos][0] == "." || alleles[chr][pos][1] == "." {
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
	runCommand("plink", []string{"--vcf", path, "--make-rel", "square", "--out", "results/linking/plink"})
}

func plinkKing(path string) {
	runCommand("plink", []string{"--vcf", path, "--make-king", "square", "--out", "results/linking/plink"})
}

func retrievePositionAlleles(chr, pos string) (string, string) {
	ref, alt := retrievePositionAllelesDataset(chr, pos, data.GG)
	if ref == "." || alt == "." {
		ref, alt = retrievePositionAllelesDataset(chr, pos, data.RL)
	}
	return ref, alt
}

func retrievePositionAllelesDataset(chr, pos, dataset string) (string, string) {
	prg := "bcftools"
	args := []string{
		"query",
		"-f", "%POS\t%REF\t%ALT\n",
		"-r", fmt.Sprintf("%s:%s-%s", chr, pos, pos),
		data.GetChromosomeFilepath(chr, dataset),
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

func removeFiles(filepaths ...string) {
	for _, path := range filepaths {
		path = strings.TrimSpace(path)
		if _, err := os.Stat(path); os.IsNotExist(err) {
			continue
		}
		err := os.Remove(path)
		if err != nil {
			log.Printf("Error removing file %s: %v\n", path, err)
		}
	}
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

func recoveryAccuracy(afType string) {
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
		retrieved := getIndividualPositionSNPs(chr, positions, mainSamples, data.GG)
		for idv = range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
		retrieved = getIndividualPositionSNPs(chr, positions, relatives, data.RL)
		for idv = range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
	}
	//
	references := getMajorGenotypesForGuessedLoci()
	ancestries := data.LoadAncestry()
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
				locus = data.MergeLocus(chr, pos)
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
	resultFolder := "results/recoveryAccuracy"
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
		err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		err = p.LoadStats(data.GG)
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
				float32(Log3(p.FindMaxAbsoluteWeight()*math.Pow(10, float64(p.WeightPrecision))))
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
	folder := "results/recoveryOutput"
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
					chr, pos := data.SplitLocus(locus)
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
	samples := solver.All1000GenomesSamples()
	var filepath string
	var outFile *os.File
	var err error
	var ok bool
	var positions []string
	var chrStr, pos, idv, anc string
	var positionMap map[string]struct{}
	frequencies := make(map[string]map[string]map[string][]float32)
	ancestries := data.LoadAncestry()
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
		retrieved := getIndividualPositionSNPs(chrStr, positions, samples, data.GG)
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
	filepath = "info/genotype_frequencies.json"
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
	folder := "info"
	filepath := path.Join(folder, "genotype_frequencies.json")
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
		data.GetChromosomeFilepath(chr, dataset),
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
				gt, err = data.NormalizeSnp(fields[1])
				if err != nil {
					log.Printf("Error normalizing SNP: %s -- %v\n", line, err)
					continue
				}
				snp, err = data.SnpToSum(gt)
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
		data.GetChromosomeFilepath(chr, dataset),
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

func loadValidatedPgsAndLoci() (map[string]int, map[string][]string, error) {
	file, err := os.Open("results/validated_pgs.json")
	if err != nil {
		log.Println("Error opening validated ids file")
		return nil, nil, err
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	var idsToNumVariants map[string]int
	err = decoder.Decode(&idsToNumVariants)
	if err != nil {
		log.Println("Error decoding validated ids")
		return nil, nil, err
	}
	file, err = os.Open("results/validated_loci.json")
	if err != nil {
		log.Println("Error opening validated ids file")
		return nil, nil, err
	}
	defer file.Close()
	decoder = json.NewDecoder(file)
	var lociToPgs map[string][]string
	err = decoder.Decode(&lociToPgs)
	if err != nil {
		log.Println("Error decoding validated loci")
		return nil, nil, err
	}
	return idsToNumVariants, lociToPgs, nil
}

func Log3(x float64) float64 {
	return math.Log2(x) / math.Log2(3)
}
