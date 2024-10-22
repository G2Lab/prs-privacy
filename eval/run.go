package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"path"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"sync"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
	"github.com/nikirill/prs/tools"
	"gonum.org/v1/gonum/stat/distuv"
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
	case "impute":
		imputeWorkflow()
	case "accuracy":
		guessAccuracy(perAncestryAf)
	case "gaccuracy":
		guessAccuracy(globalAf)
	case "genfreq":
		calculateGenotypeFrequenciesOnlyGuessed()
	case "king":
		linkingWithKing()
	case "ibd":
		prepareIBD()
	case "uniqueness_gg":
		uniquenessExperiment(tools.GG)
	case "uniqueness_uk":
		uniquenessExperiment(tools.UKB)
	case "scores":
		calculateScoresAndStats()
	case "rounding":
		precisionToUniqueness()
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
	//	p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
	//	if p.WeightPrecision >= 17 {
	//		maxw := p.FindMaxAbsoluteWeight() * math.Pow(10, float64(p.WeightPrecision))
	//		density := float64(p.NumVariants) / tools.Log3(maxw)
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
	dataset := tools.UKB
	pgsID := "PGS000869"
	p := pgs.NewPGS()
	err := p.LoadCatalogFile(path.Join(params.LocalDataFolder, pgsID+"_hmPOS_GRCh37.txt"))
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
		density := float64(p.NumVariants) / tools.Log3(float64(solver.FindMaxAbsoluteBigInt(weights)))
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
	resultFolder := "results/defense"
	filepath := path.Join(resultFolder, fmt.Sprintf("%s.json", pgsID))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	if err = encoder.Encode(results); err != nil {
		log.Fatal("Cannot encode json", err)
	}
	filepath = path.Join(resultFolder, "scores.json")
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

	precisionEstimationLimit := 15
	numWorkers := 1
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

func estimateNumUniqueScores(p *pgs.PGS, M int, powerUpperBound int) (float64, int, uint64) {
	numPossibleSubsets := calculateNumPossibleSubsets(p.Weights)
	minScoreDecimal, maxScoreDecimal, secondMinDecimal, secondMaxDecimal := p.FindMinAndMaxScores()
	minVal, _ := minScoreDecimal.Float64()
	maxVal, _ := maxScoreDecimal.Float64()
	secondMin, _ := secondMinDecimal.Float64()
	secondMax, _ := secondMaxDecimal.Float64()
	predictedMeans, predictedStds := p.EstimateMeanAndStd()
	normalDist := distuv.Normal{Mu: predictedMeans["ALL"], Sigma: predictedStds["ALL"]}
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

	type result struct {
		unique float64
		sets   map[uint16]uint16
	}
	numWorkers := 21
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
		sm = slv.SolveDeterministic()
	} else {
		sm = slv.SolveProbabilistic()
	}
	return solver.SortByLikelihood(sm, p.PopulationStats[ppl], p.EffectAlleles)
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
