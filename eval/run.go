package main

import (
	"bufio"
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
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"syscall"
	"time"

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
	case "density":
		measureDensityAccuracyEffect()
	case "resources":
		measureResources()
	case "confusion":
		genotypeConfusionMatrix()
	case "ukblinking":
		ukbLinking()
	case "ukblinkingstats":
		ukbLinkingStats()
	case "ukballeles":
		retrieveUkbAlleles()
	case "ukbking":
		ukbKingRelatedness()
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
	for id := range idsToNumVariants {
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
			fmt.Printf("GuessAccuracy: %.0f%%\n", 100*solver.Accuracy(cohort[individual].Genotype, solutions[0]))
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
	for id := range pgsToNumVariants {
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
	for id := range pgsToNumVariants {
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
	for id := range pgsToNumVariants {
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
		for id := range idsToNumVariants {
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
		for id := range idsToNumVariants {
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

// It samples one individual from each decile of accuracy for each ancestry group
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
		//fmt.Printf("--- %s ---\n", ancestry)
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
			//fmt.Printf("%s: %s\n", inputs[index].Individual, strconv.FormatFloat(inputs[index].GuessAccuracy, 'f', 3, 64))
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
	Individual         string
	Ancestry           string
	SNPs               int
	GuessAccuracy      float32
	ReferenceAccuracy  float32
	StochasticAccuracy float32
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

	// Per-locus AAF map for stochastic prediction (the probability of 1)
	locusAAF := make(map[string]map[string]float32) // locus -> ancestry -> AAF
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
		for i, locus := range p.Loci {
			if _, ok = observedLoci[locus]; ok {
				continue
			}
			observedLoci[locus] = struct{}{}
			locusAAF[locus] = make(map[string]float32)
			for ppl := range p.PopulationStats {
				locusAAF[locus][ppl] = p.PopulationStats[ppl].AF[i][1]
			}
		}
	}

	rng := rand.New(rand.NewSource(int64('p')<<16 | int64('r')<<8 | int64('s')))
	idvResults := make([]*IdvResult, 0)
	locusResults := make(map[string]*LocusResult)
	var locus string
	for idv = range guessed {
		res := IdvResult{Individual: idv, Ancestry: pgs.GetIndividualAncestry(idv, ancestries),
			SNPs: 0, GuessAccuracy: 0, ReferenceAccuracy: 0, StochasticAccuracy: 0}
		for chr = range guessed[idv] {
			for pos = range guessed[idv][chr] {
				// Individual accuracy
				if guessed[idv][chr][pos] == trueSnps[idv][chr][pos] {
					res.GuessAccuracy++
				}
				if guessed[idv][chr][pos] == references[res.Ancestry][chr][pos] {
					res.ReferenceAccuracy++
				}
				// Stochastic prediction
				locus = data.MergeLocus(chr, pos)
				if aafMap, aafOk := locusAAF[locus]; aafOk {
					pAlt := float64(aafMap[res.Ancestry])
					if pAlt <= 0 {
						pAlt = 0.001
					} else if pAlt >= 1 {
						pAlt = 0.999
					}
					// Two independent allele draws, each 1 with probability pAlt
					var predicted uint8
					if rng.Float64() < pAlt {
						predicted++
					}
					if rng.Float64() < pAlt {
						predicted++
					}
					if predicted == trueSnps[idv][chr][pos] {
						res.StochasticAccuracy++
					}
				}
				res.SNPs++
				// Locus accuracy
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
		res.StochasticAccuracy /= float32(res.SNPs)
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

	observedLoci = make(map[string]struct{})
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
	var args []string
	if dataset == data.UKB {
		args = []string{
			"query",
			"-S", data.UKBBSamplesFile,
			"-f", "%POS~[%SAMPLE=%GT\t]\n",
			data.GetChromosomeFilepath(chr, dataset),
			"-r", "",
		}
	} else {
		args = []string{
			"query",
			"-s", strings.Join(individuals, ","),
			"-f", "%POS~[%SAMPLE=%GT\t]\n",
			data.GetChromosomeFilepath(chr, dataset),
			"-r", "",
		}
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
			fields := strings.SplitN(line, "~", 2)
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

func getIndividualPositionSNPsBatched(chr string, positions []string, dataset string) map[string]map[string]uint8 {
	// Build comma-separated regions for a single bcftools call
	regions := make([]string, len(positions))
	for i, pos := range positions {
		regions[i] = fmt.Sprintf("%s:%s-%s", chr, pos, pos)
	}
	prg := "bcftools"
	var args []string
	if dataset == data.UKB {
		args = []string{
			"query",
			"-S", data.UKBBSamplesFile,
			"-f", "%POS~[%SAMPLE=%GT\t]\n",
			"-r", strings.Join(regions, ","),
			data.GetChromosomeFilepath(chr, dataset),
		}
	} else {
		log.Fatalf("getIndividualPositionSNPsBatched only supports UKB dataset")
	}

	posSet := make(map[string]struct{}, len(positions))
	for _, pos := range positions {
		posSet[pos] = struct{}{}
	}

	output, err := exec.Command(prg, args...).Output()
	if err != nil {
		log.Fatalf("Error executing bcftools command for chr %s: %v", chr, err)
	}

	var pos, gt, idv string
	var ok bool
	var snp uint8
	alleles := make(map[string]map[string]uint8)
	lines := strings.Split(string(output), "\n")
	fmt.Printf("  Chr %s: %d lines in output\n", chr, len(lines))
	for _, line := range lines {
		if line == "" {
			continue
		}
		fields := strings.SplitN(line, "~", 2)
		if len(fields) < 2 {
			log.Printf("Not enough fields in line: %s\n", line)
			continue
		}
		pos = fields[0]
		if _, ok = posSet[pos]; !ok {
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
				continue
			}
			snp, err = data.SnpToSum(gt)
			alleles[idv][pos] = snp
		}
	}
	fmt.Printf("  Chr %s: %d unique individuals\n", chr, len(alleles))
	return alleles
}

func getIndividualPositionStrings(chr string, positions []string, individuals []string, dataset string) map[string]map[string]string {
	prg := "bcftools"
	var args []string
	if dataset == data.UKB {
		args = []string{
			"query",
			"-S", data.UKBBSamplesFile,
			"-f", "%POS~[%SAMPLE=%GT\t]\n",
			data.GetChromosomeFilepath(chr, dataset),
			"-r", "",
		}
	} else {
		args = []string{
			"query",
			"-s", strings.Join(individuals, ","),
			"-f", "%POS~[%SAMPLE=%GT\t]\n",
			data.GetChromosomeFilepath(chr, dataset),
			"-r", "",
		}
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
			fields := strings.SplitN(line, "~", 2)
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

func measureDensityAccuracyEffect() {
	densityThreshold := 4.0
	variantThreshold := 30
	idsToNumVariants, _, err := loadValidatedPgsAndLoci()
	if err != nil {
		log.Println(err)
		return
	}

	individuals := getIndividualsSample()
	anc := os.Args[2]
	id, err := strconv.Atoi(os.Args[3])
	if err != nil || id < 0 || id >= len(individuals[anc]) {
		log.Printf("Invalid individual index: %s", os.Args[3])
	}
	idv := individuals[anc][id]
	fmt.Printf("=== Ancestry: %s, Idv: %d (%s) ===\n", anc, id, individuals[anc][id])

	type Result struct {
		PgsID            string
		NumVariants      int
		Individual       string
		Ancestry         string
		Precision        uint32
		Density          float64
		GuessAccuracy    float32
		BaselineAccuracy float32
		Latency          float64
		CPUTime          float64
		MemoryPeak       uint64
	}
	results := make([]Result, 0)

	allPgs := make([]string, 0, len(idsToNumVariants))
	for id := range idsToNumVariants {
		allPgs = append(allPgs, id)
	}
	// sort by number of variants
	sort.Slice(allPgs, func(i, j int) bool {
		return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
	})

	for _, pgsID := range allPgs {
		if idsToNumVariants[pgsID] > variantThreshold {
			break
		}
		fmt.Printf("Processing %s: %d variants\n", pgsID, idsToNumVariants[pgsID])
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			continue
		}
		p.LoadStats(data.GG)
		cohort := solver.NewCohort(p, data.GG)
		var precisionLimit uint32 = 9
		if len(p.Loci) >= DeterminismLimit {
			precisionLimit = 8
		}

		originalWeightPrecision := p.WeightPrecision
		originalWeights := make([]*apd.Decimal, len(p.Weights))
		for i, w := range p.Weights {
			originalWeights[i] = new(apd.Decimal).Set(w)
			fmt.Printf("%s\t", p.Weights[i].String())
		}
		fmt.Println()

		roundingCtx := apd.BaseContext.WithPrecision(p.Context.Precision)
		roundingCtx.Rounding = apd.RoundHalfUp
		// Vary density by gradually rounding weights and recalculating the target scores
		for prec := originalWeightPrecision; prec > 0; prec-- {
			// Calculate density based on the reduced precision
			maxw := p.FindMaxAbsoluteWeight() * math.Pow(10, float64(prec))
			density := float64(p.NumVariants) / Log3(maxw)
			if density > densityThreshold || density <= 0 {
				break
			}

			if prec < p.WeightPrecision {
				for i := range p.Weights {
					_, err = roundingCtx.Quantize(p.Weights[i], originalWeights[i], -int32(prec))
					if err != nil {
						log.Printf("Error rounding weight: %v\n", err)
					}
				}
				p.WeightPrecision = prec
			}
			fmt.Printf("%d: %.2f\t", prec, density)

			idvScore := solver.CalculateDecimalScore(p.Context, cohort[idv].Genotype, p.Weights, p.EffectAlleles)
			slv := solver.NewDP(idvScore, p, anc, precisionLimit, make(map[int]uint8))

			var usage syscall.Rusage
			// Memory measurement setup
			var maxMem uint64
			done := make(chan struct{})
			var wgMem sync.WaitGroup

			// Force GC to get a clean baseline
			runtime.GC()
			var m runtime.MemStats
			runtime.ReadMemStats(&m)
			maxMem = m.HeapAlloc

			wgMem.Add(1)
			go func() {
				defer wgMem.Done()
				ticker := time.NewTicker(time.Millisecond)
				defer ticker.Stop()
				var m runtime.MemStats
				for {
					select {
					case <-done:
						return
					case <-ticker.C:
						runtime.ReadMemStats(&m)
						if m.HeapAlloc > maxMem {
							maxMem = m.HeapAlloc
						}
					}
				}
			}()

			start := time.Now()
			syscall.Getrusage(syscall.RUSAGE_SELF, &usage)
			cpuStart := float64(usage.Utime.Sec) + float64(usage.Utime.Usec)/1e6 + float64(usage.Stime.Sec) + float64(usage.Stime.Usec)/1e6

			var solutions [][]uint8
			if p.NumVariants < DeterminismLimit {
				sm := slv.SolveDeterministic()
				solutions = solver.SortByLikelihood(sm, p.PopulationStats[anc], p.EffectAlleles)
			} else {
				sm := slv.SolveProbabilistic()
				solutions = solver.SortByLikelihood(sm, p.PopulationStats[anc], p.EffectAlleles)
			}
			clockElapsed := time.Since(start).Seconds()
			syscall.Getrusage(syscall.RUSAGE_SELF, &usage)
			cpuEnd := float64(usage.Utime.Sec) + float64(usage.Utime.Usec)/1e6 + float64(usage.Stime.Sec) + float64(usage.Stime.Usec)/1e6
			cpuElapsed := cpuEnd - cpuStart

			close(done)
			wgMem.Wait()

			acc := float32(0.0)
			if len(solutions) > 0 {
				acc = solver.Accuracy(solutions[0], cohort[idv].Genotype)
			}

			// Calculate major genotype accuracy
			majorAcc := float32(0.0)
			for i := 0; i < len(p.Loci); i++ {
				majorG := uint8(0)
				maxFreq := float32(-1.0)

				// Check frequency of 0, 1, 2
				for g := 0; g <= 2; g++ {
					if p.PopulationStats[anc].GF[i][g] > maxFreq {
						maxFreq = p.PopulationStats[anc].GF[i][g]
						majorG = uint8(g)
					}
				}

				// Compare with true genotype
				trueG := cohort[idv].Genotype[pgs.Ploidy*i] + cohort[idv].Genotype[pgs.Ploidy*i+1]
				if majorG == trueG {
					majorAcc++
				}
			}
			majorAcc /= float32(len(p.Loci))

			results = append(results, Result{
				PgsID:            pgsID,
				NumVariants:      idsToNumVariants[pgsID],
				Individual:       idv,
				Ancestry:         anc,
				Precision:        prec,
				Density:          density,
				GuessAccuracy:    acc,
				BaselineAccuracy: majorAcc,
				Latency:          clockElapsed,
				CPUTime:          cpuElapsed,
				MemoryPeak:       maxMem,
			})
		}
		fmt.Println()
	}

	filepath := path.Join("results/density/",
		fmt.Sprintf("%d_density_accuracy_%s_%d.json", int(densityThreshold), anc, id))
	resFile, err := os.Create(filepath)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Density accuracy measurement completed.")
}

func measureResources() {
	idsToNumVariants, _, err := loadValidatedPgsAndLoci()
	if err != nil {
		log.Println(err)
		return
	}

	individuals := getIndividualsSample()
	anc := os.Args[2]
	id, err := strconv.Atoi(os.Args[3])
	if err != nil || id < 0 || id >= len(individuals[anc]) {
		log.Printf("Invalid individual index: %s", os.Args[3])
	}
	fmt.Printf("=== Ancestry: %s, Idv: %d (%s) ===\n", anc, id, individuals[anc][id])

	type Result struct {
		PgsID         string
		NumVariants   int
		Individual    string
		Ancestry      string
		GuessAccuracy float32
		WallTime      float64
		CPUTime       float64
		MemoryPeak    uint64
	}
	results := make([]Result, 0)

	allPgs := make([]string, 0, len(idsToNumVariants))
	for id := range idsToNumVariants {
		allPgs = append(allPgs, id)
	}
	// sort by number of variants
	sort.Slice(allPgs, func(i, j int) bool {
		return idsToNumVariants[allPgs[i]] < idsToNumVariants[allPgs[j]]
	})

	for _, pgsID := range allPgs {
		fmt.Printf("Processing %s: %d variants\n", pgsID, idsToNumVariants[pgsID])
		p := pgs.NewPGS()
		err = p.LoadCatalogFile(path.Join(data.LocalInputFolder, pgsID+"_hmPOS_GRCh37.txt"))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			continue
		}
		p.LoadStats(data.GG)
		cohort := solver.NewCohort(p, data.GG)
		var precisionLimit uint32 = 9
		if len(p.Loci) >= DeterminismLimit {
			precisionLimit = 8
		}
		idv := individuals[anc][id]

		// Force GC to get a clean baseline
		runtime.GC()
		// Memory measurement setup
		var maxMem uint64
		done := make(chan struct{})
		var wgMem sync.WaitGroup
		var usage syscall.Rusage

		var m runtime.MemStats
		runtime.ReadMemStats(&m)
		maxMem = m.HeapAlloc
		start := time.Now()
		syscall.Getrusage(syscall.RUSAGE_SELF, &usage)
		cpuStart := float64(usage.Utime.Sec) + float64(usage.Utime.Usec)/1e6 + float64(usage.Stime.Sec) + float64(usage.Stime.Usec)/1e6

		wgMem.Add(1)
		go func() {
			defer wgMem.Done()
			ticker := time.NewTicker(time.Millisecond)
			defer ticker.Stop()
			var m runtime.MemStats
			for {
				select {
				case <-done:
					return
				case <-ticker.C:
					runtime.ReadMemStats(&m)
					if m.HeapAlloc > maxMem {
						maxMem = m.HeapAlloc
					}
				}
			}
		}()

		var solutions [][]uint8
		slv := solver.NewDP(cohort[idv].Score, p, anc, precisionLimit, make(map[int]uint8))
		if p.NumVariants < DeterminismLimit {
			sm := slv.SolveDeterministic()
			solutions = solver.SortByLikelihood(sm, p.PopulationStats[anc], p.EffectAlleles)
		} else {
			sm := slv.SolveProbabilistic()
			solutions = solver.SortByLikelihood(sm, p.PopulationStats[anc], p.EffectAlleles)
		}
		clockElapsed := time.Since(start).Seconds()
		syscall.Getrusage(syscall.RUSAGE_SELF, &usage)
		cpuEnd := float64(usage.Utime.Sec) + float64(usage.Utime.Usec)/1e6 + float64(usage.Stime.Sec) + float64(usage.Stime.Usec)/1e6
		cpuElapsed := cpuEnd - cpuStart

		close(done)
		wgMem.Wait()

		acc := float32(0.0)
		if len(solutions) > 0 {
			acc = solver.Accuracy(solutions[0], cohort[idv].Genotype)
		}

		results = append(results, Result{
			PgsID:         pgsID,
			NumVariants:   idsToNumVariants[pgsID],
			Individual:    idv,
			Ancestry:      anc,
			GuessAccuracy: acc,
			WallTime:      clockElapsed,
			CPUTime:       cpuElapsed,
			MemoryPeak:    maxMem,
		})
	}

	os.MkdirAll("results/resources/", 0755)
	filepath := path.Join("results/resources/",
		fmt.Sprintf("resources_%s_%d.json", anc, id))
	resFile, err := os.Create(filepath)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Println("Resources measurement completed.")
}

func genotypeConfusionMatrix() {
	fmt.Println("--- Computing genotype confusion matrix ---")

	guessed := loadGuessedGenotypes(perAncestryAf)

	positionsPerChr := make(map[string]map[string]struct{})
	for _, chrMap := range guessed {
		for chr, posMap := range chrMap {
			if _, ok := positionsPerChr[chr]; !ok {
				positionsPerChr[chr] = make(map[string]struct{})
			}
			for pos := range posMap {
				positionsPerChr[chr][pos] = struct{}{}
			}
		}
	}

	mainSamples := solver.All1000GenomesSamples()
	relatives := solver.AllRelativeSamples()
	trueSnps := make(map[string]map[string]map[string]uint8) // individual -> chr -> pos -> genotype
	for _, idv := range solver.All1000GenomesAndRelativeSamples() {
		trueSnps[idv] = make(map[string]map[string]uint8)
	}
	for chrID := 1; chrID <= 22; chrID++ {
		chr := strconv.Itoa(chrID)
		posMap, ok := positionsPerChr[chr]
		if !ok {
			continue
		}
		positions := make([]string, 0, len(posMap))
		for pos := range posMap {
			positions = append(positions, pos)
		}
		fmt.Printf("Chromosome %s: %d positions\n", chr, len(positions))
		retrieved := getIndividualPositionSNPs(chr, positions, mainSamples, data.GG)
		for idv := range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
		retrieved = getIndividualPositionSNPs(chr, positions, relatives, data.RL)
		for idv := range retrieved {
			trueSnps[idv][chr] = retrieved[idv]
		}
	}

	// counts[locus][trueGenotype][guessedGenotype] = count
	type ErrorLikelihood struct {
		Counts        [3][3]int     `json:"counts"`
		Probabilities [3][3]float64 `json:"probabilities"`
	}
	locusErrors := make(map[string]*ErrorLikelihood)

	for idv, chrMap := range guessed {
		if _, ok := trueSnps[idv]; !ok {
			continue
		}
		for chr, posMap := range chrMap {
			for pos, guessedGtp := range posMap {
				trueGtp, ok := trueSnps[idv][chr][pos]
				if !ok {
					continue
				}
				if guessedGtp > 2 || trueGtp > 2 {
					continue
				}
				locus := data.MergeLocus(chr, pos)
				if _, exists := locusErrors[locus]; !exists {
					locusErrors[locus] = &ErrorLikelihood{}
				}
				locusErrors[locus].Counts[trueGtp][guessedGtp]++
			}
		}
	}

	// Normalize each row to get probabilities
	for _, el := range locusErrors {
		for trueGtp := 0; trueGtp <= 2; trueGtp++ {
			rowSum := 0
			for guessedGtp := 0; guessedGtp <= 2; guessedGtp++ {
				rowSum += el.Counts[trueGtp][guessedGtp]
			}
			if rowSum > 0 {
				for guessedGtp := 0; guessedGtp <= 2; guessedGtp++ {
					el.Probabilities[trueGtp][guessedGtp] = float64(el.Counts[trueGtp][guessedGtp]) / float64(rowSum)
				}
			}
		}
	}

	resultFolder := "results"
	outPath := path.Join(resultFolder, "confusion_matrix.json")
	resFile, err := os.OpenFile(outPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	encoder.SetIndent("", "  ")
	if err = encoder.Encode(locusErrors); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	fmt.Printf("Saved confusion matrix for %d loci to %s\n", len(locusErrors), outPath)
}

func ukbLinking() {
	fmt.Println("--- UKB perturbed genotype linking ---")

	type ConfusionEntry struct {
		Counts        [3][3]int     `json:"counts"`
		Probabilities [3][3]float64 `json:"probabilities"`
	}
	cmFile, err := os.Open("results/confusion_matrix.json")
	if err != nil {
		log.Fatalf("Error opening confusion matrix: %v", err)
	}
	defer cmFile.Close()
	var confusionMatrix map[string]*ConfusionEntry
	if err = json.NewDecoder(cmFile).Decode(&confusionMatrix); err != nil {
		log.Fatalf("Error decoding confusion matrix: %v", err)
	}
	fmt.Printf("Loaded confusion matrix for %d loci\n", len(confusionMatrix))

	// Clamp probabilities: replace 0 with 0.01, 1 with 0.99, then normalize each row
	for _, entry := range confusionMatrix {
		for row := 0; row < 3; row++ {
			for col := 0; col < 3; col++ {
				if entry.Probabilities[row][col] == 0 {
					entry.Probabilities[row][col] = 0.01
				} else if entry.Probabilities[row][col] == 1 {
					entry.Probabilities[row][col] = 0.99
				}
			}
			// Normalize the row
			rowSum := 0.0
			for col := 0; col < 3; col++ {
				rowSum += entry.Probabilities[row][col]
			}
			if rowSum > 0 {
				for col := 0; col < 3; col++ {
					entry.Probabilities[row][col] /= rowSum
				}
			}
		}
	}

	// Group loci by chromosome
	positionsPerChr := make(map[string][]string)
	locusList := make([]string, 0, len(confusionMatrix))
	for locus := range confusionMatrix {
		locusList = append(locusList, locus)
		chr, pos := data.SplitLocus(locus)
		positionsPerChr[chr] = append(positionsPerChr[chr], pos)
	}
	// Sort positions within each chromosome
	for chr := range positionsPerChr {
		sort.Strings(positionsPerChr[chr])
	}

	// Retrieve true genotypes for all UK Biobank individuals
	individuals := solver.AllUKBiobankSamples()
	sort.Strings(individuals)
	fmt.Printf("Number of UKB individuals: %d\n", len(individuals))

	numWorkers := 12

	// trueGenotypes: individual -> locus -> genotype (0, 1, or 2)
	trueGenotypes := make(map[string]map[string]uint8)
	var mu sync.Mutex
	var wg sync.WaitGroup
	chrChan := make(chan int, numWorkers)
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chrID := range chrChan {
				chr := strconv.Itoa(chrID)
				positions, ok := positionsPerChr[chr]
				if !ok || len(positions) == 0 {
					continue
				}
				fmt.Printf("Chromosome %s: retrieving %d positions\n", chr, len(positions))
				retrieved := getIndividualPositionSNPsBatched(chr, positions, data.UKB)
				mu.Lock()
				for idv, posMap := range retrieved {
					if _, exists := trueGenotypes[idv]; !exists {
						trueGenotypes[idv] = make(map[string]uint8)
					}
					for pos, gtp := range posMap {
						locus := data.MergeLocus(chr, pos)
						trueGenotypes[idv][locus] = gtp
					}
				}
				mu.Unlock()
			}
		}()
	}
	for chrID := 1; chrID <= 22; chrID++ {
		chrChan <- chrID
	}
	close(chrChan)
	wg.Wait()
	fmt.Printf("Retrieved genotypes for %d individuals\n", len(trueGenotypes))

	// Apply genotype errors across individuals
	perturbedGenotypes := make(map[string]map[string]uint8, len(trueGenotypes))
	for idv, locusMap := range trueGenotypes {
		copied := make(map[string]uint8, len(locusMap))
		for locus, gtp := range locusMap {
			copied[locus] = gtp
		}
		perturbedGenotypes[idv] = copied
	}
	idvChan := make(chan string, numWorkers)
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			localRng := rand.New(rand.NewSource(rand.Int63()))
			for idv := range idvChan {
				perturbed := perturbedGenotypes[idv]
				for locus, trueGtp := range perturbed {
					entry, ok := confusionMatrix[locus]
					if !ok || trueGtp > 2 {
						continue
					}
					// Sample from the confusion matrix row
					r := localRng.Float64()
					cumulative := 0.0
					newGtp := uint8(2) // default to last column
					for col := 0; col < 3; col++ {
						cumulative += entry.Probabilities[trueGtp][col]
						if r < cumulative {
							newGtp = uint8(col)
							break
						}
					}
					if newGtp != trueGtp {
						perturbed[locus] = newGtp
					}
				}
			}
		}()
	}
	for idv := range perturbedGenotypes {
		idvChan <- idv
	}
	close(idvChan)
	wg.Wait()
	fmt.Printf("Applied genotype errors to %d individuals\n", len(perturbedGenotypes))

	// Load REF/ALT alleles from pre-computed file
	allelesFile := "results/ukb_alleles.json"
	af, err := os.Open(allelesFile)
	if err != nil {
		log.Fatalf("Error opening alleles file %s: %v. Run with -e=ukballeles first.", allelesFile, err)
	}
	var alleles map[string][]string
	if err = json.NewDecoder(af).Decode(&alleles); err != nil {
		log.Fatalf("Error decoding alleles file: %v", err)
	}
	af.Close()
	fmt.Printf("Loaded REF/ALT alleles for %d loci from %s\n", len(alleles), allelesFile)

	// Sort loci by chromosome and position for VCF output
	sort.Slice(locusList, func(i, j int) bool {
		chrI, posI := data.SplitLocus(locusList[i])
		chrJ, posJ := data.SplitLocus(locusList[j])
		chrIInt, _ := strconv.Atoi(chrI)
		chrJInt, _ := strconv.Atoi(chrJ)
		if chrIInt != chrJInt {
			return chrIInt < chrJInt
		}
		posIInt, _ := strconv.Atoi(posI)
		posJInt, _ := strconv.Atoi(posJ)
		return posIInt < posJInt
	})

	resultDir := "results/ukb_linking"
	if err = os.MkdirAll(resultDir, 0755); err != nil {
		log.Fatalf("Error creating result directory: %v", err)
	}

	vcfPath := path.Join(resultDir, "plink.vcf")
	vcfFile, err := os.OpenFile(vcfPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("Cannot create VCF file: %v", err)
	}
	defer vcfFile.Close()

	// Write VCF header
	writer := bufio.NewWriterSize(vcfFile, 64*1024*1024) // 64MB buffer
	fmt.Fprint(writer, "##fileformat=VCFv4.1\n")
	fmt.Fprint(writer, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")

	// Samples header: perturbed ($idv) then truth (idv)
	fmt.Fprint(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
	for _, idv := range individuals {
		fmt.Fprintf(writer, "\t$%s", idv)
	}
	for _, idv := range individuals {
		fmt.Fprintf(writer, "\t%s", idv)
	}
	fmt.Fprint(writer, "\n")
	writer.Flush()
	fmt.Printf("Saved the sample IDs into the VCF file\n")

	// Write genotype lines
	genotypeStr := func(gtp uint8) string {
		switch gtp {
		case 0:
			return "0|0"
		case 1:
			return "1/0"
		case 2:
			return "1|1"
		default:
			return ".|."
		}
	}

	for _, locus := range locusList {
		if alleles[locus][0] == "." || alleles[locus][1] == "." {
			continue
		}
		chr, pos := data.SplitLocus(locus)
		fmt.Fprintf(writer, "%s\t%s\t.\t%s\t%s\t100\tPASS\t.\tGT", chr, pos, alleles[locus][0], alleles[locus][1])

		// Perturbed genotypes
		for _, idv := range individuals {
			if gtp, ok := perturbedGenotypes[idv][locus]; ok {
				fmt.Fprintf(writer, "\t%s", genotypeStr(gtp))
			} else {
				fmt.Fprint(writer, "\t0|0")
			}
		}

		// Truth genotypes
		for _, idv := range individuals {
			if gtp, ok := trueGenotypes[idv][locus]; ok {
				fmt.Fprintf(writer, "\t%s", genotypeStr(gtp))
			} else {
				fmt.Fprint(writer, "\t0|0")
			}
		}

		fmt.Fprint(writer, "\n")
	}
	writer.Flush()
	vcfFile.Close()
	fmt.Printf("VCF written to %s\n", vcfPath)

	// Free large data structures before running plink
	trueGenotypes = nil
	perturbedGenotypes = nil
	alleles = nil
	runtime.GC()
	fmt.Println("Memory released")
	// Compress, index, and run KING
	compressedVcfPath := vcfPath + ".gz"
	compressVcfFile(vcfPath, compressedVcfPath)
	indexVcfFile(compressedVcfPath)

	runCommand("plink", []string{
		"--vcf", compressedVcfPath,
		"--make-king", "triangle",
		"--out", path.Join(resultDir, "plink"),
	})
	fmt.Printf("KING relatedness matrix saved to %s.king\n", path.Join(resultDir, "plink"))

	// Cleanup intermediate files
	removeFiles(vcfPath, compressedVcfPath, compressedVcfPath+".csi")
	fmt.Println("Done.")
}

func retrieveUkbAlleles() {
	fmt.Println("--- Retrieving UKB REF/ALT alleles ---")

	// Load confusion matrix to get the loci
	cmFile, err := os.Open("results/confusion_matrix.json")
	if err != nil {
		log.Fatalf("Error opening confusion matrix: %v", err)
	}
	var confusionMatrix map[string]json.RawMessage
	if err = json.NewDecoder(cmFile).Decode(&confusionMatrix); err != nil {
		log.Fatalf("Error decoding confusion matrix: %v", err)
	}
	cmFile.Close()

	// Group positions by chromosome
	positionsPerChr := make(map[string][]string)
	for locus := range confusionMatrix {
		chr, pos := data.SplitLocus(locus)
		positionsPerChr[chr] = append(positionsPerChr[chr], pos)
	}
	for chr := range positionsPerChr {
		sort.Strings(positionsPerChr[chr])
	}

	numWorkers := 4
	alleles := make(map[string][]string)
	var mu sync.Mutex
	var wg sync.WaitGroup
	chrChan := make(chan int, numWorkers)
	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for chrID := range chrChan {
				chr := strconv.Itoa(chrID)
				positions, ok := positionsPerChr[chr]
				if !ok || len(positions) == 0 {
					continue
				}
				fmt.Printf("Chromosome %s: retrieving alleles for %d positions\n", chr, len(positions))
				regions := make([]string, len(positions))
				for i, pos := range positions {
					regions[i] = fmt.Sprintf("%s:%s-%s", chr, pos, pos)
				}
				args := []string{
					"query",
					"-f", "%POS\t%REF\t%ALT\n",
					"-r", strings.Join(regions, ","),
					data.GetChromosomeFilepath(chr, data.UKB),
				}
				output, err := exec.Command("bcftools", args...).Output()
				if err != nil {
					log.Fatalf("Error executing bcftools for chr %s: %v", chr, err)
				}
				posSet := make(map[string]struct{}, len(positions))
				for _, pos := range positions {
					posSet[pos] = struct{}{}
				}
				lines := strings.Split(string(output), "\n")
				localAlleles := make(map[string][]string)
				for _, line := range lines {
					if line == "" {
						continue
					}
					fields := strings.Split(line, "\t")
					if len(fields) < 3 {
						continue
					}
					pos := fields[0]
					if _, ok := posSet[pos]; !ok {
						continue
					}
					locus := data.MergeLocus(chr, pos)
					if _, exists := localAlleles[locus]; !exists {
						localAlleles[locus] = []string{fields[1], fields[2]}
					}
				}
				for _, pos := range positions {
					locus := data.MergeLocus(chr, pos)
					if _, exists := localAlleles[locus]; !exists {
						localAlleles[locus] = []string{".", "."}
					}
				}
				mu.Lock()
				for locus, al := range localAlleles {
					alleles[locus] = al
				}
				mu.Unlock()
			}
		}()
	}
	for chrID := 1; chrID <= 22; chrID++ {
		chrChan <- chrID
	}
	close(chrChan)
	wg.Wait()
	fmt.Printf("Retrieved alleles for %d loci\n", len(alleles))

	outPath := "results/ukb_alleles.json"
	outFile, err := os.OpenFile(outPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("Error creating output file: %v", err)
	}
	defer outFile.Close()
	encoder := json.NewEncoder(outFile)
	encoder.SetIndent("", "  ")
	if err = encoder.Encode(alleles); err != nil {
		log.Fatalf("Error encoding alleles: %v", err)
	}
	fmt.Printf("Saved alleles to %s\n", outPath)
}

func ukbKingRelatedness() {
	fmt.Println("--- UKB KING relatedness ground truth ---")

	if len(os.Args) < 4 {
		log.Fatalf("Usage: %s -e=ukbking <bed_prefix> <snp_ids_file>\n"+
			"  bed_prefix:    path prefix for .bed/.bim/.fam files\n"+
			"  snp_ids_file:  file with one SNP ID per line to extract", os.Args[0])
	}
	bedPrefix := os.Args[2]
	snpFile := os.Args[3]

	resultDir := "results/ukb_king"
	if err := os.MkdirAll(resultDir, 0755); err != nil {
		log.Fatalf("Error creating result directory: %v", err)
	}
	outPrefix := path.Join(resultDir, "king")

	// Full square kinship matrix
	fmt.Println("Computing KING kinship matrix...")
	runCommand("plink", []string{
		"--bfile", bedPrefix,
		"--extract", snpFile,
		"--memory", "100000",
		"--make-king", "square",
		"--out", outPrefix,
	})
	fmt.Printf("KING matrix saved to %s.king\n", outPrefix)

	// Relationship table: twins, 1st and 2nd degree (kinship >= 0.0884)
	fmt.Println("Computing KING relationship table...")
	runCommand("plink", []string{
		"--bfile", bedPrefix,
		"--extract", snpFile,
		"--memory", "100000",
		"--make-king-table",
		"--king-table-filter", "0.0884",
		"--out", outPrefix,
	})

	// Post-process: classify relationships by degree
	kinFile := outPrefix + ".kin0"
	kinData, err := os.ReadFile(kinFile)
	if err != nil {
		log.Fatalf("Error reading relationship table %s: %v", kinFile, err)
	}
	lines := strings.Split(strings.TrimSpace(string(kinData)), "\n")
	if len(lines) < 2 {
		fmt.Println("No related pairs found.")
		fmt.Println("Done.")
		return
	}

	annotatedFile := outPrefix + "_relationships.tsv"
	out, err := os.Create(annotatedFile)
	if err != nil {
		log.Fatalf("Error creating annotated file: %v", err)
	}
	defer out.Close()

	// Write header with added Degree column
	fmt.Fprintf(out, "%s\tDegree\n", lines[0])

	twinCount, firstCount, secondCount := 0, 0, 0
	for _, line := range lines[1:] {
		fields := strings.Split(line, "\t")
		// KINSHIP is the last column in plink2 .kin0 output
		kinship, err := strconv.ParseFloat(fields[len(fields)-1], 64)
		if err != nil {
			log.Printf("Error parsing kinship value: %s", fields[len(fields)-1])
			continue
		}
		var degree string
		switch {
		case kinship > 0.354:
			degree = "Twin"
			twinCount++
		case kinship > 0.177:
			degree = "1st"
			firstCount++
		default:
			degree = "2nd"
			secondCount++
		}
		fmt.Fprintf(out, "%s\t%s\n", line, degree)
	}
	fmt.Printf("Relationships: %d twin/duplicate, %d 1st-degree, %d 2nd-degree\n",
		twinCount, firstCount, secondCount)
	fmt.Printf("Annotated relationship table saved to %s\n", annotatedFile)

	fmt.Println("Done.")
}

// ============================================================
// UKB Linking Statistics — streaming KING triangle matrix reader
// ============================================================

// BoxPlotStats holds pre-computed boxplot statistics for one category.
type BoxPlotStats struct {
	Label  string    `json:"label"`
	Method string    `json:"method"`
	N      int64     `json:"n"`
	Mean   float64   `json:"mean"`
	Median float64   `json:"median"`
	Q1     float64   `json:"q1"`
	Q3     float64   `json:"q3"`
	WhisLo float64   `json:"whislo"`
	WhisHi float64   `json:"whishi"`
	Fliers []float64 `json:"fliers"`
}

// streamingHistogram estimates quantiles from a stream of bounded values
// using a fixed-width histogram. KING kinship values lie in [-1, 1].
const (
	sqNumBins = 20000
	sqMin     = -1.0
	sqMax     = 1.0
)

type streamingHistogram struct {
	bins  []int64
	count int64
	sum   float64
}

func newStreamingHistogram() *streamingHistogram {
	return &streamingHistogram{bins: make([]int64, sqNumBins)}
}

func (sh *streamingHistogram) add(val float64) {
	sh.count++
	sh.sum += val
	bin := int((val - sqMin) / ((sqMax - sqMin) / float64(sqNumBins)))
	if bin < 0 {
		bin = 0
	}
	if bin >= sqNumBins {
		bin = sqNumBins - 1
	}
	sh.bins[bin]++
}

func (sh *streamingHistogram) merge(other *streamingHistogram) {
	for i := range sh.bins {
		sh.bins[i] += other.bins[i]
	}
	sh.count += other.count
	sh.sum += other.sum
}

func (sh *streamingHistogram) histPercentile(p float64) float64 {
	target := int64(math.Ceil(p / 100.0 * float64(sh.count)))
	if target <= 0 {
		target = 1
	}
	binWidth := (sqMax - sqMin) / float64(sqNumBins)
	cumulative := int64(0)
	for i, c := range sh.bins {
		cumulative += c
		if cumulative >= target {
			return sqMin + (float64(i)+0.5)*binWidth
		}
	}
	return sqMax
}

func (sh *streamingHistogram) computeBoxPlotStats(label, method string) *BoxPlotStats {
	if sh.count == 0 {
		return nil
	}
	binWidth := (sqMax - sqMin) / float64(sqNumBins)
	q1 := sh.histPercentile(25)
	med := sh.histPercentile(50)
	q3 := sh.histPercentile(75)
	iqr := q3 - q1
	whisLo := q1 - 1.5*iqr
	whisHi := q3 + 1.5*iqr

	// Find actual whisker endpoints and collect fliers from bins
	actualWhisLo, actualWhisHi := sqMin, sqMax
	whisLoFound, whisHiFound := false, false
	fliers := make([]float64, 0)

	for i := 0; i < sqNumBins; i++ {
		if sh.bins[i] == 0 {
			continue
		}
		binCenter := sqMin + (float64(i)+0.5)*binWidth
		if binCenter >= whisLo && !whisLoFound {
			actualWhisLo = binCenter
			whisLoFound = true
		} else if binCenter < whisLo {
			// Lower fliers
			for j := int64(0); j < sh.bins[i] && len(fliers) < 500; j++ {
				fliers = append(fliers, binCenter)
			}
		}
	}
	for i := sqNumBins - 1; i >= 0; i-- {
		if sh.bins[i] == 0 {
			continue
		}
		binCenter := sqMin + (float64(i)+0.5)*binWidth
		if binCenter <= whisHi && !whisHiFound {
			actualWhisHi = binCenter
			whisHiFound = true
		} else if binCenter > whisHi {
			// Upper fliers
			for j := int64(0); j < sh.bins[i] && len(fliers) < 500; j++ {
				fliers = append(fliers, binCenter)
			}
		}
	}

	return &BoxPlotStats{
		Label:  label,
		Method: method,
		N:      sh.count,
		Mean:   sh.sum / float64(sh.count),
		Median: med,
		Q1:     q1,
		Q3:     q3,
		WhisLo: actualWhisLo,
		WhisHi: actualWhisHi,
		Fliers: fliers,
	}
}

// statsAccumulator collects values for a single (category, method) key.
// Small categories store all values; large ones use a streaming histogram.
type statsAccumulator struct {
	label     string
	method    string
	values    []float64           // exact storage for small categories
	streaming *streamingHistogram // histogram-based streaming stats for large categories
}

func newExactAccumulator(label, method string) *statsAccumulator {
	return &statsAccumulator{label: label, method: method, values: make([]float64, 0, 1024)}
}

func newStreamingAccumulator(label, method string) *statsAccumulator {
	return &statsAccumulator{label: label, method: method, streaming: newStreamingHistogram()}
}

func (sa *statsAccumulator) add(val float64) {
	if sa.streaming != nil {
		sa.streaming.add(val)
	} else {
		sa.values = append(sa.values, val)
	}
}

func (sa *statsAccumulator) count() int64 {
	if sa.streaming != nil {
		return sa.streaming.count
	}
	return int64(len(sa.values))
}

func (sa *statsAccumulator) computeStats() *BoxPlotStats {
	if sa.streaming != nil {
		return sa.streaming.computeBoxPlotStats(sa.label, sa.method)
	}
	data := sa.values
	if len(data) == 0 {
		return nil
	}
	sort.Float64s(data)
	n := len(data)
	q1 := percentile(data, 25)
	med := percentile(data, 50)
	q3 := percentile(data, 75)
	iqr := q3 - q1
	whisLo := q1 - 1.5*iqr
	whisHi := q3 + 1.5*iqr
	actualWhisLo := data[0]
	for _, v := range data {
		if v >= whisLo {
			actualWhisLo = v
			break
		}
	}
	actualWhisHi := data[n-1]
	for i := n - 1; i >= 0; i-- {
		if data[i] <= whisHi {
			actualWhisHi = data[i]
			break
		}
	}
	fliers := make([]float64, 0)
	for _, v := range data {
		if v < actualWhisLo || v > actualWhisHi {
			fliers = append(fliers, v)
			if len(fliers) >= 500 {
				break
			}
		}
	}
	s := 0.0
	for _, v := range data {
		s += v
	}
	return &BoxPlotStats{
		Label:  sa.label,
		Method: sa.method,
		N:      int64(len(data)),
		Mean:   s / float64(len(data)),
		Median: med,
		Q1:     q1,
		Q3:     q3,
		WhisLo: actualWhisLo,
		WhisHi: actualWhisHi,
		Fliers: fliers,
	}
}

func percentile(sorted []float64, p float64) float64 {
	if len(sorted) == 0 {
		return 0
	}
	k := (p / 100) * float64(len(sorted)-1)
	f := math.Floor(k)
	c := math.Ceil(k)
	if f == c {
		return sorted[int(k)]
	}
	return sorted[int(f)]*(c-k) + sorted[int(c)]*(k-f)
}

func getBaseID(fullID string) string {
	clean := fullID
	if strings.HasPrefix(clean, "$") {
		clean = clean[1:]
	}
	parts := strings.SplitN(clean, "_", 2)
	if len(parts) == 2 && parts[0] == parts[1] {
		return parts[0]
	}
	return clean
}

func readKingIDs(filePath string) []string {
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatalf("Cannot open KING ID file %s: %v", filePath, err)
	}
	defer f.Close()
	var ids []string
	scanner := bufio.NewScanner(f)
	first := true
	for scanner.Scan() {
		if first {
			first = false // skip header
			continue
		}
		ids = append(ids, strings.TrimSpace(scanner.Text()))
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("Error reading KING ID file: %v", err)
	}
	return ids
}

type ukbRelationship struct {
	relBaseID string
	degree    string // "Twin", "1", "2"
}

func readUKBRelatedIndividualsGo(filePath string) map[string][]ukbRelationship {
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatalf("Cannot open relationship file %s: %v", filePath, err)
	}
	defer f.Close()
	reader := bufio.NewReader(f)
	headerLine, err := reader.ReadString('\n')
	if err != nil {
		log.Fatalf("Cannot read relationship header: %v", err)
	}
	header := strings.Split(strings.TrimSpace(headerLine), "\t")
	iid1Col, iid2Col, degreeCol := -1, -1, -1
	for i, h := range header {
		switch h {
		case "#IID1", "IID1":
			iid1Col = i
		case "IID2":
			iid2Col = i
		case "Degree":
			degreeCol = i
		}
	}
	if iid1Col < 0 || iid2Col < 0 || degreeCol < 0 {
		log.Fatalf("Missing columns in relationship file: IID1=%d IID2=%d Degree=%d", iid1Col, iid2Col, degreeCol)
	}
	related := make(map[string][]ukbRelationship)
	scanner := bufio.NewScanner(reader)
	scanner.Buffer(make([]byte, 1024*1024), 1024*1024)
	for scanner.Scan() {
		fields := strings.Split(scanner.Text(), "\t")
		if len(fields) <= degreeCol {
			continue
		}
		id1, id2, deg := fields[iid1Col], fields[iid2Col], fields[degreeCol]
		if deg != "Twin" && deg != "1st" && deg != "2nd" {
			continue
		}
		related[id1] = append(related[id1], ukbRelationship{relBaseID: id2, degree: deg})
		related[id2] = append(related[id2], ukbRelationship{relBaseID: id1, degree: deg})
	}
	return related
}

// workerResult holds per-category accumulated values from one worker.
type workerResult struct {
	// Small categories: key = "category|method" → values
	small map[string][]float64
	// Unrelated streaming histograms
	unrelKING     *streamingHistogram
	unrelTrueKING *streamingHistogram
}

func processTriangleRows(
	dataFile string,
	rowOffsets []struct {
		sampleIdx int
		offset    int64
	},
	numPerturbed int,
	relatedOf [][]struct {
		localIdx int
		degree   string
	},
	workerID int,
) *workerResult {
	res := &workerResult{
		small:         make(map[string][]float64),
		unrelKING:     newStreamingHistogram(),
		unrelTrueKING: newStreamingHistogram(),
	}

	f, err := os.Open(dataFile)
	if err != nil {
		log.Fatalf("Worker %d: cannot open data file: %v", workerID, err)
	}
	defer f.Close()

	for ri, ro := range rowOffsets {
		si := ro.sampleIdx
		k := si - numPerturbed // local truth index

		if _, err := f.Seek(ro.offset, 0); err != nil {
			log.Fatalf("Worker %d: seek error: %v", workerID, err)
		}

		// Read line
		reader := bufio.NewReaderSize(f, 32*1024*1024)
		line, err := reader.ReadBytes('\n')
		if err != nil && len(line) == 0 {
			log.Fatalf("Worker %d: read error at row %d: %v", workerID, si, err)
		}
		line = bytes.TrimRight(line, "\n\r")

		// Build related lookup for this truth sample
		perturbedRelated := make(map[int]string) // perturbed localIdx → degree
		truthRelated := make(map[int]string)     // truth localIdx → degree
		for _, r := range relatedOf[k] {
			perturbedRelated[r.localIdx] = r.degree
			truthRelated[r.localIdx] = r.degree
		}

		// Parse fields and process
		fieldIdx := 0
		start := 0
		for i := 0; i <= len(line); i++ {
			if i == len(line) || line[i] == '\t' {
				if fieldIdx < si { // si columns expected (0..si-1)
					val, perr := strconv.ParseFloat(string(line[start:i]), 64)
					if perr != nil {
						start = i + 1
						fieldIdx++
						continue
					}
					if fieldIdx < numPerturbed {
						// Perturbed column
						if fieldIdx == k {
							// Self pair: perturbed[k] vs truth[k]
							key := "Self|KING"
							res.small[key] = append(res.small[key], val)
						} else if deg, ok := perturbedRelated[fieldIdx]; ok {
							key := degreeToCategory(deg) + "|KING"
							res.small[key] = append(res.small[key], val)
						} else {
							// Unrelated
							res.unrelKING.add(val)
						}
					} else {
						// Truth column (fieldIdx >= numPerturbed, fieldIdx < si)
						truthLocalIdx := fieldIdx - numPerturbed
						if deg, ok := truthRelated[truthLocalIdx]; ok {
							key := degreeToCategory(deg) + "|TrueKING"
							res.small[key] = append(res.small[key], val)
						} else {
							// TrueKING unrelated
							res.unrelTrueKING.add(val)
						}
					}
				}
				start = i + 1
				fieldIdx++
			}
		}

		if ri > 0 && ri%5000 == 0 {
			fmt.Printf("  Worker %d: processed %d/%d rows\n", workerID, ri, len(rowOffsets))
		}
	}
	return res
}

func degreeToCategory(deg string) string {
	switch deg {
	case "Twin":
		return "Twin"
	case "1st":
		return "1st degree"
	case "2nd":
		return "2nd degree"
	}
	return deg
}

func ukbLinkingStats() {
	fmt.Println("--- UKB Linking Statistics (streaming) ---")

	resultDir := "results/ukb_linking"
	idFile := path.Join(resultDir, "plink.king.id")
	dataFile := path.Join(resultDir, "plink.king")
	relFile := "results/ukb_king/king_relationships.tsv"

	// 1. Read sample IDs
	ids := readKingIDs(idFile)
	n := len(ids)
	fmt.Printf("Total samples: %d\n", n)

	// Separate perturbed and truth IDs
	var perturbedIDs, truthIDs []string
	for _, sid := range ids {
		if strings.HasPrefix(sid, "$") {
			perturbedIDs = append(perturbedIDs, sid)
		} else {
			truthIDs = append(truthIDs, sid)
		}
	}
	numPerturbed := len(perturbedIDs)
	numTruth := len(truthIDs)
	fmt.Printf("Perturbed: %d, Truth: %d\n", numPerturbed, numTruth)

	// 2. Read ground-truth relationships (keyed by base IDs like "-100")
	related := readUKBRelatedIndividualsGo(relFile)
	fmt.Printf("Loaded relationships for %d individuals\n", len(related))

	// 3. Build baseID → truth local index mapping
	baseToLocalIdx := make(map[string]int, numTruth)
	for k, sid := range truthIDs {
		baseToLocalIdx[getBaseID(sid)] = k
	}

	// Pre-compute per-truth-sample related info
	type relInfo struct {
		localIdx int
		degree   string
	}
	relatedOf := make([][]relInfo, numTruth)
	for k := 0; k < numTruth; k++ {
		baseID := getBaseID(truthIDs[k])
		rels := related[baseID]
		relatedOf[k] = make([]relInfo, 0, len(rels))
		for _, r := range rels {
			rIdx, ok := baseToLocalIdx[r.relBaseID]
			if !ok {
				continue
			}
			relatedOf[k] = append(relatedOf[k], relInfo{localIdx: rIdx, degree: r.degree})
		}
	}

	// 4. Index byte offsets for truth-sample rows in the triangle file
	// Triangle format: file line f (0-indexed) → sample index f+1, with f+1 columns.
	// Truth samples have global indices [numPerturbed, n-1], so file lines [numPerturbed-1, n-2].
	fmt.Println("Indexing byte offsets for truth-sample rows...")
	type rowOffset struct {
		sampleIdx int
		offset    int64
	}
	truthRowOffsets := make([]rowOffset, 0, numTruth)
	func() {
		f, err := os.Open(dataFile)
		if err != nil {
			log.Fatalf("Cannot open KING data file: %v", err)
		}
		defer f.Close()
		reader := bufio.NewReaderSize(f, 32*1024*1024)
		fileLine := 0
		firstTruthLine := numPerturbed - 1 // file line for first truth sample (index numPerturbed)
		lastTruthLine := n - 2             // file line for last truth sample (index n-1)
		byteOffset := int64(0)
		for fileLine <= lastTruthLine {
			if fileLine >= firstTruthLine {
				sampleIdx := fileLine + 1
				truthRowOffsets = append(truthRowOffsets, rowOffset{sampleIdx: sampleIdx, offset: byteOffset})
			}
			line, err := reader.ReadBytes('\n')
			byteOffset += int64(len(line))
			if err != nil {
				break
			}
			fileLine++
			if fileLine%100000 == 0 {
				fmt.Printf("  Scanned %d/%d file lines, found %d truth rows\n", fileLine, n-1, len(truthRowOffsets))
			}
		}
	}()
	fmt.Printf("Indexed %d truth-sample rows\n", len(truthRowOffsets))

	// 5. Parallel processing
	numWorkers := runtime.NumCPU()
	if numWorkers > 12 {
		numWorkers = 12
	}
	fmt.Printf("Using %d workers\n", numWorkers)

	// Convert relatedOf to the type expected by processTriangleRows
	type relInfoExport struct {
		localIdx int
		degree   string
	}
	relatedOfExport := make([][]struct {
		localIdx int
		degree   string
	}, numTruth)
	for k := 0; k < numTruth; k++ {
		relatedOfExport[k] = make([]struct {
			localIdx int
			degree   string
		}, len(relatedOf[k]))
		for j, r := range relatedOf[k] {
			relatedOfExport[k][j] = struct {
				localIdx int
				degree   string
			}{r.localIdx, r.degree}
		}
	}

	chunkSize := (len(truthRowOffsets) + numWorkers - 1) / numWorkers
	var wg sync.WaitGroup
	results := make([]*workerResult, numWorkers)
	for w := 0; w < numWorkers; w++ {
		start := w * chunkSize
		end := start + chunkSize
		if start >= len(truthRowOffsets) {
			break
		}
		if end > len(truthRowOffsets) {
			end = len(truthRowOffsets)
		}
		wg.Add(1)
		go func(workerID, s, e int) {
			defer wg.Done()
			workerRows := make([]struct {
				sampleIdx int
				offset    int64
			}, e-s)
			for i, ro := range truthRowOffsets[s:e] {
				workerRows[i] = struct {
					sampleIdx int
					offset    int64
				}{ro.sampleIdx, ro.offset}
			}
			results[workerID] = processTriangleRows(
				dataFile, workerRows, numPerturbed, relatedOfExport, workerID,
			)
			fmt.Printf("  Worker %d finished (%d rows)\n", workerID, e-s)
		}(w, start, end)
	}
	wg.Wait()
	fmt.Println("All workers finished")

	// 6. Merge results
	categories := []string{"Self", "Twin", "1st degree", "2nd degree"}
	methods := []string{"KING", "TrueKING"}

	accumulators := make(map[string]*statsAccumulator)
	for _, cat := range categories {
		for _, m := range methods {
			key := cat + "|" + m
			accumulators[key] = newExactAccumulator(cat, m)
		}
	}
	accumulators["Unrelated|KING"] = newStreamingAccumulator("Unrelated", "KING")
	accumulators["Unrelated|TrueKING"] = newStreamingAccumulator("Unrelated", "TrueKING")
	// TrueSelf is always 0.5
	accumulators["Self|TrueKING"] = newExactAccumulator("Self", "TrueKING")

	for _, wr := range results {
		if wr == nil {
			continue
		}
		// Merge small category values
		for key, vals := range wr.small {
			acc, ok := accumulators[key]
			if !ok {
				parts := strings.SplitN(key, "|", 2)
				accumulators[key] = newExactAccumulator(parts[0], parts[1])
				acc = accumulators[key]
			}
			for _, v := range vals {
				acc.add(v)
			}
		}
		// Merge unrelated histograms
		accumulators["Unrelated|KING"].streaming.merge(wr.unrelKING)
		accumulators["Unrelated|TrueKING"].streaming.merge(wr.unrelTrueKING)
	}

	// Add TrueSelf values (diagonal = 0.5, not stored in triangle matrix)
	for i := 0; i < numTruth; i++ {
		accumulators["Self|TrueKING"].add(0.5)
	}

	// 7. Compute and output boxplot stats
	var allStats []*BoxPlotStats
	outputOrder := []string{
		"Self|KING", "Self|TrueKING",
		"Twin|KING", "Twin|TrueKING",
		"1st degree|KING", "1st degree|TrueKING",
		"2nd degree|KING", "2nd degree|TrueKING",
		"Unrelated|KING", "Unrelated|TrueKING",
	}
	for _, key := range outputOrder {
		acc, ok := accumulators[key]
		if !ok || acc.count() == 0 {
			continue
		}
		stats := acc.computeStats()
		if stats != nil {
			allStats = append(allStats, stats)
			fmt.Printf("  %s (%s): n=%d, median=%.4f, Q1=%.4f, Q3=%.4f, mean=%.4f\n",
				stats.Label, stats.Method, stats.N, stats.Median, stats.Q1, stats.Q3, stats.Mean)
		}
	}

	outPath := path.Join(resultDir, "kinship_stats.json")
	outFile, err := os.OpenFile(outPath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("Cannot create output file: %v", err)
	}
	defer outFile.Close()
	enc := json.NewEncoder(outFile)
	enc.SetIndent("", "  ")
	if err := enc.Encode(allStats); err != nil {
		log.Fatalf("Cannot encode stats: %v", err)
	}
	fmt.Printf("Stats written to %s\n", outPath)
	fmt.Println("Done.")
}
