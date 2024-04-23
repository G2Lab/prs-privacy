package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
	"path"
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

const DeterminismLimit = 32

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
		pprof.WriteHeapProfile(f)
		defer f.Close()
	}

	//evaluateGA()
	//scoreDistribution()
	//likelihoodEffect()
	//scoreToLikelihoodDistribution()
	//scoreToLikelihood()
	//samples()
	//distribution()
	//evaluateReferences()
	//accuracyLikelihood()
	//buildDPTables()
	//likelihoodWeight()
	//sortingChoice()
	//sortingChoiceParallel()
	//accuracyParallel()
	//findAllSolutions()
	sequentialSolving()
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

func sequentialSolving() {
	//thresholds := []int{500, 1000, 1500, 2000, 2500}
	thresholds := []int{500, 1000, 1500, 2000}
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

	pgsToLoci := make(map[string]map[string]struct{})
	for locus, pgsIDs := range lociToPgs {
		for _, pgsID := range pgsIDs {
			if _, ok := pgsToLoci[pgsID]; !ok {
				pgsToLoci[pgsID] = make(map[string]struct{})
			}
			pgsToLoci[pgsID][locus] = struct{}{}
		}
	}

	//INDIVIDUAL := "HG00119"
	individuals := []string{"HG00119", "HG00524", "HG00581", "HG00656", "HG00731", "HG01936", "HG02025", "HG02026",
		"HG02067", "HG02353", "HG02371", "HG02250", "HG02373", "HG02386", "HG02375", "HG03713", "HG03673", "NA19238",
		"NA19239", "NA19027", "NA19334", "NA19331", "NA19664", "NA19678", "NA19661", "NA19713", "NA20321", "NA20334",
		"NA20289", "NA20792", "NA20868", "NA20895"}

	chunkNum, chunkSize := getChunkInfo(len(individuals))
	fmt.Println(individuals[chunkNum*chunkSize : (chunkNum+1)*chunkSize])

	type Result struct {
		Individual string
		Accuracies map[int]float32
	}
	results := make([]*Result, 0)
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
		sort.Slice(allPgs, func(i, j int) bool {
			return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
		})

		populations := tools.LoadAncestry()
		indPop := populations[individual]
		if strings.Contains(indPop, ",") {
			indPop = strings.Split(indPop, ",")[0]
		}
		var pgsID string
		var solmap map[string][]uint8
		var solutions [][]uint8
		var accuracyWithMissed, accuracyWithoutMissed float32
		guessedSnps := make(map[string]uint8)
		trueSnpsWithMissed := make(map[string]uint8)    // count all the snps
		trueSnpsWithoutMissed := make(map[string]uint8) // count only the snps that were "guessed", i.e., no solutions does not count
		accuraciesWithMissed := make(map[int]float32)
		accuraciesWithoutMissed := make(map[int]float32)
		thresholdPos := 0
		for {
			pgsID = allPgs[0]
			p := pgs.NewPGS()
			err = p.LoadCatalogFile(path.Join(params.DataFolder, pgsID+"_hmPOS_GRCh38.txt"))
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
			slv := solver.NewDP(cohort[individual].Score, p, indPop, recoveredSnps)
			if len(p.Loci)-len(recoveredSnps) < DeterminismLimit {
				solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihood)
			} else {
				solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
			}
			solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
			if len(solutions) == 0 && len(recoveredSnps) > 0 {
				fmt.Println("^^^ No solutions with extra loci, calculating from scratch ^^^")
				slv := solver.NewDP(cohort[individual].Score, p, indPop, make(map[int]uint8))
				if len(p.Loci) < DeterminismLimit {
					solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihood)
				} else {
					solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
				}
				solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
			}
			if len(solutions) > 0 {
				fmt.Printf("Top solution accuracy: %.3f\n", solver.Accuracy(solutions[0], cohort[individual].Genotype))
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
				if _, ok := trueSnpsWithMissed[locus]; ok {
					continue
				}
				trueSnpsWithMissed[locus] = cohort[individual].Genotype[pgs.Ploidy*i] + cohort[individual].Genotype[pgs.Ploidy*i+1]
				if len(solutions) > 0 {
					trueSnpsWithoutMissed[locus] = trueSnpsWithMissed[locus]
				}
			}
			if len(guessedSnps) >= thresholds[thresholdPos] {
				accuracyWithMissed, accuracyWithoutMissed = 0, 0
				for locus, snp := range trueSnpsWithMissed {
					if guessed, ok := guessedSnps[locus]; ok && snp == guessed {
						accuracyWithMissed += 1
					}
				}
				accuraciesWithMissed[thresholds[thresholdPos]] = accuracyWithMissed / float32(len(trueSnpsWithMissed))
				//
				for locus, snp := range trueSnpsWithoutMissed {
					if guessed, ok := guessedSnps[locus]; ok && snp == guessed {
						accuracyWithoutMissed += 1
					}
				}
				accuraciesWithoutMissed[thresholds[thresholdPos]] = accuracyWithoutMissed / float32(len(trueSnpsWithoutMissed))
				fmt.Printf("### Accuracy %d: w/ missed %.3f, w/o missed %.3f\n", thresholds[thresholdPos],
					accuracyWithMissed/float32(len(trueSnpsWithMissed)), accuracyWithoutMissed/float32(len(trueSnpsWithoutMissed)))
				thresholdPos++
				if thresholdPos == len(thresholds) {
					break
				}
			}
			fmt.Printf("Guessed %d\n", len(guessedSnps))
			if len(guessedSnps) > thresholds[len(thresholds)-1] {
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
			sort.Slice(allPgs, func(i, j int) bool {
				return len(pgsToLoci[allPgs[i]]) < len(pgsToLoci[allPgs[j]])
			})
		}
		results = append(results, &Result{Individual: individual, Accuracies: accuraciesWithMissed})
	}

	resultFolder := "results/sequential"
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
	//var accuracyWithMissed float32 = 0.0
	//for locus, snp := range trueSnpsWithMissed {
	//	if guessed, ok := guessedSnps[locus]; ok && snp == guessed {
	//		accuracyWithMissed += 1
	//	}
	//}
	//fmt.Printf("Accuracy: %.3f\n", accuracyWithMissed/float32(len(trueSnpsWithMissed)))
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
	catalogFile := "PGS001835_hmPOS_GRCh38.txt"
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
	individuals := solver.All1000GenomeSamples()

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
			solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihood)
		} else {
			solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
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
	//INDIVIDUAL := "NA12286"
	//
	//INDIVIDUAL := "HG01028"
	//INDIVIDUAL := "NA18531"
	//INDIVIDUAL := "NA20872"
	//INDIVIDUAL := "NA20507"
	//INDIVIDUAL := "NA07037"
	//INDIVIDUAL := "HG03015"
	//INDIVIDUAL := "HG01253"
	INDIVIDUAL := "HG00119"

	p := pgs.NewPGS()
	//catalogFile := "PGS000154_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000753_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000043_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000307_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000534_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000011_hmPOS_GRCh38.txt"
	//catalogFile := "PGS003436_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002264_hmPOS_GRCh38.txt"
	//catalogFile := "PGS003760_hmPOS_GRCh38.txt"
	//catalogFile := "PGS004246_hmPOS_GRCh38.txt"
	//catalogFile := "PGS004249_hmPOS_GRCh38.txt"
	catalogFile := "PGS001828_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
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
	indPop := populations[INDIVIDUAL]
	if strings.Contains(indPop, ",") {
		indPop = strings.Split(indPop, ",")[0]
	}
	//recovered := map[int]uint8{0: 0, 1: 1, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 1, 8: 2, 10: 1, 11: 0}
	recovered := make(map[int]uint8)
	slv := solver.NewDP(cohort[INDIVIDUAL].Score, p, indPop, recovered)

	majorReference := solver.AllReferenceAlleleSample(p.PopulationStats[indPop].AF)
	fmt.Printf("Accuracy with reference: %f\n",
		solver.Accuracy(majorReference, cohort[INDIVIDUAL].Genotype))

	fmt.Printf("Effect loci: ")
	for i := 0; i < len(cohort[INDIVIDUAL].Genotype); i += 2 {
		if cohort[INDIVIDUAL].Genotype[i] == p.EffectAlleles[i/2] && cohort[INDIVIDUAL].Genotype[i+1] == p.EffectAlleles[i/2] {
			fmt.Printf("%d ", i+1)
		}
		if cohort[INDIVIDUAL].Genotype[i]+cohort[INDIVIDUAL].Genotype[i+1] == 1 {
			fmt.Printf("%d ", i)
		}
	}
	fmt.Println()
	//solmap := slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
	solmap := slv.SolveFromScratchDeterministic(solver.UseLikelihood)
	//solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	//solutions := solver.SortByLikelihood(solmap, p.PopulationEAF[indPop])
	solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)

	fmt.Printf("\nTrue:\n%s, %.2f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
		solver.CalculateFullSequenceLikelihood(cohort[INDIVIDUAL].Genotype, p.PopulationStats[indPop].AF, p.EffectAlleles))
	//solver.CalculateSolutionSpectrumDistance(cohort[INDIVIDUAL].Genotype, p.PopulationStats[indPop], p.EffectAlleles))

	fmt.Printf("Guessed %d:\n", len(solutions))
	//target := solver.ScoreToTarget(cohort[INDIVIDUAL].Score, p)
	target := cohort[INDIVIDUAL].Score
	for _, solution := range solutions {
		diff := new(apd.Decimal)
		p.Context.Sub(diff, target, solver.CalculateDecimalScore(p.Context, solution, p.Weights, p.EffectAlleles))
		fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
			diff.String(), solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles))
		//	solver.CalculateSolutionSpectrumDistance(solution, p.PopulationStats[indPop], p.EffectAlleles))
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
//	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000753_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000154_hmPOS_GRCh38.txt"
//	catalogFile := "PGS003760_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000851_hmPOS_GRCh38.txt"
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
//	var count int
//	for _, ppl := range pgs.POPULATIONS {
//		score := minScore
//		count = 0
//		for i := 0; i < len(sortedSamples); i++ {
//			if sortedScores[i] < score || populations[sortedSamples[i]] != ppl || sortedScores[i] == 0 {
//				continue
//			}
//			individuals = append(individuals, sortedSamples[i])
//			score += step
//			count += 1
//			if count == samplesPerPopulation {
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
//			solver.Accuracy(majorReference, cohort[individual].Genotype))
//
//		slv := solver.NewDP(cohort[individual].Score, p, indPop)
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihood)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
//		}
//		solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
//		out.LikelihoodAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodAccuracies = append(out.LikelihoodAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseSpectrum)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseSpectrum)
//		out.SpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.SpectrumAccuracies = append(out.SpectrumAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihoodAndSpectrum)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihoodAndSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihoodAndSpectrum)
//		out.LikelihoodSpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodSpectrumAccuracies = append(out.LikelihoodSpectrumAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
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
//	//catalogFile := "PGS003760_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000753_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000154_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS003760_hmPOS_GRCh38.txt"
//	catalogFile := "PGS000851_hmPOS_GRCh38.txt"
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
//	var count int
//	for _, ppl := range pgs.POPULATIONS {
//		score := minScore
//		count = 0
//		for i := 0; i < len(sortedSamples); i++ {
//			if sortedScores[i] < score || populations[sortedSamples[i]] != ppl || sortedScores[i] == 0 {
//				continue
//			}
//			individuals = append(individuals, sortedSamples[i])
//			score += step
//			count += 1
//			if count == samplesPerPopulation {
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
//			solver.Accuracy(majorReference, cohort[individual].Genotype))
//
//		slv := solver.NewDP(cohort[individual].Score, p, indPop)
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihood)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihood)
//		}
//		solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihood)
//		out.LikelihoodAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodAccuracies = append(out.LikelihoodAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseSpectrum)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseSpectrum)
//		out.SpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.SpectrumAccuracies = append(out.SpectrumAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
//		}
//		//
//		if len(p.Loci) < DeterminismLimit {
//			solmap = slv.SolveFromScratchDeterministic(solver.UseLikelihoodAndSpectrum)
//		} else {
//			solmap = slv.SolveFromScratchProbabilistic(solver.UseLikelihoodAndSpectrum)
//		}
//		solutions = solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihoodAndSpectrum)
//		out.LikelihoodSpectrumAccuracies = make([]string, 0, len(solutions))
//		for _, solution := range solutions {
//			out.LikelihoodSpectrumAccuracies = append(out.LikelihoodSpectrumAccuracies, fmt.Sprintf("%.3f", solver.Accuracy(solution, cohort[individual].Genotype)))
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

func buildDPTables() {
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	populations := tools.LoadAncestry()
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	samples := cohort.SortByScore()

	fmt.Println(pgs.POPULATIONS)
	for _, ppl := range pgs.POPULATIONS {
		filepath := path.Join(params.DataFolder, fmt.Sprintf("%s-%s.dp", p.PgsID, ppl))
		if _, err = os.Stat(filepath); os.IsNotExist(err) {
			for i := len(samples) - 1; i > 0; i-- {
				indPop := populations[samples[i]]
				if indPop != ppl {
					continue
				}
				fmt.Printf("Creating DP table for %s using %s\n", ppl, samples[i])
				slv := solver.NewDP(cohort[samples[i]].Score, p, indPop, make(map[int]uint8))
				state := slv.BuildProbabilisticState(2, solver.UseLikelihood)
				solver.SaveProbabilisticState(state, filepath)
				state = nil
				break
			}
		}
	}
}

//func evaluateGA() {
//	//INDIVIDUAL := "NA18595"
//	INDIVIDUAL := "HG02182" // lowest score for PGS000040
//	//INDIVIDUAL := "HG02215" // highest score for PGS000040
//	//INDIVIDUAL := "HG02728" // middle 648
//	//INDIVIDUAL := "NA19780" // high 648
//	//INDIVIDUAL := "HG00551" // low 648
//	//
//	//INDIVIDUAL := "HG01028"
//	//INDIVIDUAL := "NA18531"
//	//INDIVIDUAL := "NA20872"
//
//	p := pgs.NewPGS()
//	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	catalogFile := "PGS000040_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
//	p.LoadStats()
//	cohort := solver.NewCohort(p)
//	populations := tools.LoadAncestry()
//	indPop := populations[INDIVIDUAL]
//
//	slv := solver.NewGenetic(cohort[INDIVIDUAL].Score, p, indPop)
//
//	solmap := slv.Solve()
//	//solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
//	//solutions := solver.SortByLikelihood(solmap, p.PopulationEAF[indPop])
//	solutions := solver.SortByLikelihoodAndFrequency(solmap, p.PopulationStats[indPop], p.EffectAlleles, solver.UseLikelihoodAndSpectrum)
//	fmt.Printf("\nTrue:\n%s, %.2f, %.2f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
//		solver.CalculateFullSequenceLikelihood(cohort[INDIVIDUAL].Genotype, p.PopulationStats[indPop].AF, p.EffectAlleles),
//		solver.CalculateSolutionSpectrumDistance(cohort[INDIVIDUAL].Genotype, p.PopulationStats[indPop], p.EffectAlleles))
//
//	fmt.Printf("Guessed %d:\n", len(solutions))
//	for _, solution := range solutions {
//		diff := new(apd.Decimal)
//		p.Context.Sub(diff, cohort[INDIVIDUAL].Score, solver.CalculateDecimalScore(p.Context, solution, p.Weights, p.EffectAlleles))
//		fmt.Printf("%s -- %.3f, %s, %.2f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
//			diff.String(), solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles),
//			solver.CalculateSolutionSpectrumDistance(solution, p.PopulationStats[indPop], p.EffectAlleles))
//	}
//}

//func evaluateReferences() {
//	resFolder := "results"
//	catalogFiles := []string{
//		"PGS002302_hmPOS_GRCh38.txt",
//		"PGS000639_hmPOS_GRCh38.txt",
//		"PGS000037_hmPOS_GRCh38.txt",
//		"PGS000073_hmPOS_GRCh38.txt",
//		"PGS001827_hmPOS_GRCh38.txt"}
//	// Create a result file
//	filepath := path.Join(resFolder, "references.json")
//	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
//	if err != nil {
//		log.Fatalf("Error opening result file: %v", err)
//	}
//	defer resFile.Close()
//	encoder := json.NewEncoder(resFile)
//	samples := All1000GenomeSamples()
//	references := make([]*Reference, 0, len(catalogFiles))
//	for _, catalogFile := range catalogFiles {
//		p := pgs.NewPGS()
//		err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//		if err != nil {
//			log.Printf("Error loading catalog file: %v\n", err)
//			return
//		}
//		p.LoadStats()
//		cohort := solver.NewCohort(p)
//		fmt.Printf("%s\n", p.PgsID)
//		major := p.AllReferenceAlleleSample()
//		reference := NewReference(p.PgsID, solver.CalculateDecimalScore(p.Context, major, p.Weights), p.CalculateFullSequenceLikelihood(major))
//		for _, sample := range samples {
//			reference.Accuracies = append(reference.Accuracies, fmt.Sprintf("%.3f", solver.Accuracy(major, cohort[sample].Genotype)))
//		}
//		references = append(references, reference)
//	}
//	if err = encoder.Encode(references); err != nil {
//		log.Fatal("Cannot encode json", err)
//	}
//}

func accuracyLikelihood() {
	resFolder := "results/accuracyLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
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
	samples := solver.All1000GenomeSamples()
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
		solmap := slv.SolveFromScratchDeterministic(solver.UseLikelihood)
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

func scoreToLikelihood() {
	resFolder := "results/scoreLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	populations := tools.LoadAncestry()
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	samples := solver.All1000GenomeSamples()
	// Create csv result file
	chunkNum, chunkSize := getChunkInfo(len(samples))
	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	writer := csv.NewWriter(resFile)
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d ", i)
		output := []string{samples[i], fmt.Sprintf("%.3f", cohort[samples[i]].Score)}
		indPop := populations[samples[i]]
		slv := solver.NewDP(cohort[samples[i]].Score, p, indPop, make(map[int]uint8))
		solmap := slv.SolveFromSavedProbabilistic(solver.UseLikelihood)
		solutions := solver.SortByAccuracy(solmap, cohort[samples[i]].Genotype)
		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[samples[i]].Genotype) != 1.0 {
			fmt.Printf("No solution for %s: %s, %.3f\n", samples[i], solver.ArrayToString(solutions[0]),
				solver.Accuracy(solutions[0], cohort[samples[i]].Genotype))
			continue
		}
		for _, solution := range solutions {
			output = append(output, fmt.Sprintf("%.3f",
				solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles)))
		}
		writer.Write(output)
		writer.Flush()
	}
}

//func scoreToLikelihoodDistribution() {
//	resFolder := "results/scoreLikelihood"
//	p := pgs.NewPGS()
//	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
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
//		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[sample].Genotype) != 1.0 {
//			fmt.Printf("No solution for %s: %s, %.3f\n", sample, solver.ArrayToString(solutions[0]),
//				solver.Accuracy(solutions[0], cohort[sample].Genotype))
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
//	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	fmt.Printf("%s\n", p.PgsID)
//	cohort := solver.NewCohort(p)
//	samples := All1000GenomeSamples()
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
//	majorScore := solver.CalculateDecimalScore(p.Context, majorReference, p.Weights)
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
//				firstAcc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
//				firstLikelihood = solver.CalculateFullSequenceLikelihood(solution)
//			}
//			acc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
//			if acc == 1.0 {
//				writer.Write([]string{samples[i], fmt.Sprintf("%.3f", solver.CalculateDecimalScore(p.Context, solution, p.Weights)),
//					fmt.Sprintf("%.3f", majorScore), fmt.Sprintf("%.3f", medianScore),
//					fmt.Sprintf("%.3f", firstAcc),
//					fmt.Sprintf("%.3f", solver.Accuracy(majorReference, cohort[samples[i]].Genotype)),
//					fmt.Sprintf("%d", j), fmt.Sprintf("%d", len(solutions)), fmt.Sprintf("%.3f", firstLikelihood),
//					fmt.Sprintf("%.3f", solver.CalculateFullSequenceLikelihood(solution)), fmt.Sprintf("%.3f", majorLikelihood)})
//				break solLoop
//			}
//		}
//		if acc != 1.0 {
//			diff := new(apd.Decimal)
//			p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalScore(p.Context, cohort[samples[i]].Genotype, p.Weights))
//			fmt.Printf("\nNo solution for %s\n", samples[i])
//			fmt.Printf("\nTrue:\n%s -- %s, %f\n", solver.ArrayToString(cohort[samples[i]].Genotype),
//				diff.String(),
//				solver.CalculateFullSequenceLikelihood(cohort[samples[i]].Genotype))
//			for _, sol := range solutions {
//				p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalScore(p.Context, sol, p.Weights))
//				fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(sol), solver.Accuracy(sol, cohort[samples[i]].Genotype),
//					diff.String(), solver.CalculateFullSequenceLikelihood(sol))
//			}
//		}
//	}
//}

//func samples() {
//	smp := selectSamples(10, 2504, "prs")
//	fmt.Println(smp)
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

func distribution() {
	p := pgs.NewPGS()
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	populations := tools.LoadAncestry()
	cohort := solver.NewCohort(p)
	sorted := cohort.SortByScore()
	fmt.Printf("Median score: %g\n", cohort[sorted[len(sorted)/2]].Score)
	step := 250
	//majorReference := solver.AllReferenceAlleleSample()
	//fmt.Printf("Full major %s: score %f, likelihood %f\n", solver.ArrayToString(majorReference), solver.CalculateDecimalScore(p.Context, majorReference, p.Weights), p.CalculateFullSequenceLikelihood(majorReference))
	for i := 0; i < len(sorted); i += step {
		indPop := populations[sorted[i]]
		slv := solver.NewDP(cohort[sorted[i]].Score, p, indPop, make(map[int]uint8))
		solmap := slv.SolveFromSavedProbabilistic(solver.UseLikelihood)
		solutions := solver.SortByLikelihood(solmap, p.PopulationStats[indPop].AF, p.EffectAlleles)
		//fmt.Printf("\n\n%s, %s\n", p.PgsID, sorted[i])
		//for _, solution := range solutions {
		//	fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[sorted[i]].Genotype),
		//		cohort[sorted[i]].Score-solver.CalculateDecimalScore(solution, p.Weights), p.CalculateFullSequenceLikelihood(solution))
		var acc float32 = 0.0
		var likelihood float32 = 0.0
		for j, solution := range solutions {
			if j == 0 {
				acc = solver.Accuracy(solution, cohort[sorted[i]].Genotype)
				likelihood = solver.CalculateFullSequenceLikelihood(solution, p.PopulationStats[indPop].AF, p.EffectAlleles)
			}
			if solver.Accuracy(solution, cohort[sorted[i]].Genotype) == 1.0 {
				fmt.Printf("%s, score %f -- likelihood %d/%d, first acc %.3f, first / true likelihood %f/%f\n",
					sorted[i], cohort[sorted[i]].Score, j, len(solmap), acc, likelihood,
					solver.CalculateFullSequenceLikelihood(cohort[sorted[i]].Genotype, p.PopulationStats[indPop].AF, p.EffectAlleles))
				break
			}
		}
	}
}

func likelihoodWeight() {
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000043_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000307_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000534_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000011_hmPOS_GRCh38.txt"
	//catalogFile := "PGS003436_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002264_hmPOS_GRCh38.txt"
	//catalogFile := "PGS003760_hmPOS_GRCh38.txt"
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
//	catalogFile := "PGS000073_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	p.LoadStats()
//	cohort := solver.NewCohort(p)
//	samples := All1000GenomeSamples()
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
//			score, _ := solver.CalculateDecimalScore(p.Context, cohort[sample].Genotype[indices[i]*pgs.Ploidy:indices[i+1]*pgs.Ploidy],
//				p.Weights[indices[i]:indices[i+1]]).Float64()
//			scores[i] = append(scores[i], score)
//		}
//		//fullscore, _ := solver.CalculateDecimalScore(p.Context, cohort[sample].Genotype, p.Weights).Float64()
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
