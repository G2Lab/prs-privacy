package main

import (
	"context"
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"os"
	"path"
	"strconv"

	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
)

func main() {
	//likelihoodEffect()
	//findAllSolutions()
	//scoreToLikelihoodDistribution()
	scoreToLikelihood()
	//samples()
	//distribution()
}

func scoreToLikelihood() {
	resFolder := "results/scoreLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	//fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	samples := allSamples()
	// Create csv result file
	chunkNum, chunkSize := getChunkInfo(len(samples))
	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	writer := csv.NewWriter(resFile)
	defer writer.Flush()
	ctx, cancel := context.WithCancel(context.Background())
	numThreads := getNumThreads()
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d ", i)
		output := []string{samples[i], fmt.Sprintf("%.3f", cohort[samples[i]].Score)}
		slv := solver.NewSolver(ctx, cohort[samples[i]].Score, p, numThreads)
		solmap := slv.DP(numThreads)
		solutions := solver.SortByAccuracy(solmap, cohort[samples[i]].Genotype)
		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[samples[i]].Genotype) != 1.0 {
			fmt.Printf("No solution for %s: %s, %.3f\n", samples[i], solver.ArrayToString(solutions[0]),
				solver.Accuracy(solutions[0], cohort[samples[i]].Genotype))
			continue
		}
		for _, solution := range solutions {
			output = append(output, fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)))
		}
		writer.Write(output)
	}
	cancel()
}

func scoreToLikelihoodDistribution() {
	resFolder := "results/scoreLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	//fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	sorted := cohort.SortByScore()
	numSamples := 5 * 2
	step := (cohort[sorted[len(sorted)-1]].Score - cohort[sorted[0]].Score) / (float64(numSamples) * 10 / 6)
	samples := make([]string, 0, numSamples)
	for i := 0; i < len(sorted); i++ {
		if len(samples) == numSamples {
			break
		}
		if cohort[sorted[i]].Score >= cohort[sorted[0]].Score+float64(len(samples))*step {
			samples = append(samples, sorted[i])
			score := cohort[sorted[i]].Score
			for {
				i++
				if i >= len(sorted) {
					break
				}
				if cohort[sorted[i]].Score != score {
					samples = append(samples, sorted[i])
					break
				}
			}
		}
	}
	for _, sample := range samples {
		fmt.Printf("%s(%.3f) ", sample, cohort[sample].Score)
	}
	// Create csv result file
	filepath := path.Join(resFolder, fmt.Sprintf("%s-distro.csv", p.PgsID))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	writer := csv.NewWriter(resFile)
	defer writer.Flush()
	ctx, cancel := context.WithCancel(context.Background())
	numThreads := getNumThreads()
	for _, sample := range samples {
		fmt.Printf("%s\n", sample)
		output := []string{sample, fmt.Sprintf("%.3f", cohort[sample].Score)}
		slv := solver.NewSolver(ctx, cohort[sample].Score, p, numThreads)
		solmap := slv.DP(numThreads)
		solutions := solver.SortByAccuracy(solmap, cohort[sample].Genotype)
		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[sample].Genotype) != 1.0 {
			fmt.Printf("No solution for %s: %s, %.3f\n", sample, solver.ArrayToString(solutions[0]),
				solver.Accuracy(solutions[0], cohort[sample].Genotype))
			continue
		}
		for _, solution := range solutions {
			output = append(output, fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)))
		}
		writer.Write(output)
	}
	cancel()
}

func likelihoodEffect() {
	resFolder := "results"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	samples := allSamples()
	chunkNum, chunkSize := getChunkInfo(len(samples))
	// Create csv result file
	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	writer := csv.NewWriter(resFile)
	defer writer.Flush()
	writer.Write([]string{"individual",
		"score", "major score", "median score",
		"first accuracy", "major accuracy",
		"true position", "total solutions",
		"first likelihood", "true likelihood", "major likelihood"})
	ctx, cancel := context.WithCancel(context.Background())
	numThreads := getNumThreads()
	majorReference := p.AllMajorSample()
	majorScore := solver.CalculateScore(majorReference, p.Weights)
	majorLikelihood := p.CalculateSequenceLikelihood(majorReference)
	medianScore := cohort[cohort.SortByScore()[len(cohort)/2]].Score
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d ", i)
		slv := solver.NewSolver(ctx, cohort[samples[i]].Score, p, numThreads)
		solmap := slv.DP(numThreads)
		solutions := solver.SortByLikelihood(solmap, p)
		firstAcc := 0.0
		firstLikelihood := 0.0
		acc := 0.0
	solLoop:
		for j, solution := range solutions {
			if j == 0 {
				firstAcc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
				firstLikelihood = p.CalculateSequenceLikelihood(solution)
			}
			acc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
			if acc == 1.0 {
				writer.Write([]string{samples[i], fmt.Sprintf("%.3f", solver.CalculateScore(solution, p.Weights)),
					fmt.Sprintf("%.3f", majorScore), fmt.Sprintf("%.3f", medianScore),
					fmt.Sprintf("%.3f", firstAcc),
					fmt.Sprintf("%.3f", solver.Accuracy(majorReference, cohort[samples[i]].Genotype)),
					fmt.Sprintf("%d", j), fmt.Sprintf("%d", len(solutions)), fmt.Sprintf("%.3f", firstLikelihood),
					fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)), fmt.Sprintf("%.3f", majorLikelihood)})
				break solLoop
			}
		}
		if acc != 1.0 {
			fmt.Printf("\nNo solution for %s\n", samples[i])
			fmt.Printf("\nTrue:\n%s -- %f, %f\n", solver.ArrayToString(cohort[samples[i]].Genotype),
				cohort[samples[i]].Score-solver.CalculateScore(cohort[samples[i]].Genotype, p.Weights),
				p.CalculateSequenceLikelihood(cohort[samples[i]].Genotype))
			for _, sol := range solutions {
				fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(sol), solver.Accuracy(sol, cohort[samples[i]].Genotype),
					cohort[samples[i]].Score-solver.CalculateScore(sol, p.Weights), p.CalculateSequenceLikelihood(sol))
			}
		}
	}
	cancel()
}

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
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	cohort := solver.NewCohort(p)
	sorted := cohort.SortByScore()
	fmt.Printf("Median score: %g\n", cohort[sorted[len(sorted)/2]].Score)
	ctx, cancel := context.WithCancel(context.Background())
	step := 250
	numThreads := getNumThreads()
	majorReference := p.AllMajorSample()
	fmt.Printf("Full major %s: score %f, likelihood %f\n", solver.ArrayToString(majorReference), solver.CalculateScore(majorReference, p.Weights), p.CalculateSequenceLikelihood(majorReference))
	for i := 0; i < len(sorted); i += step {
		slv := solver.NewSolver(ctx, cohort[sorted[i]].Score, p, numThreads)
		solmap := slv.DP(numThreads)
		solutions := solver.SortByLikelihood(solmap, p)
		//fmt.Printf("\n\n%s, %s\n", p.PgsID, sorted[i])
		//for _, solution := range solutions {
		//	fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[sorted[i]].Genotype),
		//		cohort[sorted[i]].Score-solver.CalculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
		acc := 0.0
		likelihood := 0.0
		for j, solution := range solutions {
			if j == 0 {
				acc = solver.Accuracy(solution, cohort[sorted[i]].Genotype)
				likelihood = p.CalculateSequenceLikelihood(solution)
			}
			if solver.Accuracy(solution, cohort[sorted[i]].Genotype) == 1.0 {
				fmt.Printf("%s, score %f -- likelihood %d/%d, first acc %.3f, first / true likelihood %f/%f\n",
					sorted[i], cohort[sorted[i]].Score, j, len(solmap), acc,
					likelihood, p.CalculateSequenceLikelihood(cohort[sorted[i]].Genotype))
				break
			}
		}
	}
	cancel()
}

func findAllSolutions() {
	//INDIVIDUAL := "NA18595"
	//INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	INDIVIDUAL := "HG00551" // low 648
	//
	//INDIVIDUAL := "HG00266"
	//INDIVIDUAL := "HG00112"

	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	catalogFile := "PGS000891_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	p.LoadStats()
	cohort := solver.NewCohort(p)

	numThreads := getNumThreads()
	ctx, cancel := context.WithCancel(context.Background())
	slv := solver.NewSolver(ctx, cohort[INDIVIDUAL].Score, p, numThreads)

	zeroReference := make([]uint8, len(p.Weights)*pgs.NumHaplotypes)
	majorReference := p.AllMajorSample()
	fmt.Printf("Accuracy with major: %f, with zero %f\n",
		solver.Accuracy(majorReference, cohort[INDIVIDUAL].Genotype), solver.Accuracy(zeroReference, cohort[INDIVIDUAL].Genotype))

	solmap := slv.DP(numThreads)
	//solmap := slv.solve(numThreads)
	//solmap := slv.recursive(numThreads)
	cancel()
	//solmap = findComplements(solmap, p, numThreads)
	solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f, %f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
		cohort[INDIVIDUAL].Score-solver.CalculateScore(cohort[INDIVIDUAL].Genotype, p.Weights),
		p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Major likelihood: %f\n", p.CalculateSequenceLikelihood(majorReference))
	sum := int64(0)
	multiplier := math.Pow(10, params.PrecisionsLimit)
	for i := 0; i < 2*len(p.Weights); i++ {
		switch cohort[INDIVIDUAL].Genotype[i] {
		case 0:
			continue
		case 1:
			sum += int64(p.Weights[i/2] * multiplier)
		}
	}
	fmt.Printf("True sum: %d\n", sum)
	fmt.Printf("Guessed %d:\n", len(solutions))
	//fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solutions[0]),
	//	solver.Accuracy(solutions[0], cohort[INDIVIDUAL].Genotype),
	//	cohort[INDIVIDUAL].Score-solver.CalculateScore(solutions[0], p.Weights), p.CalculateSequenceLikelihood(solutions[0]))
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-solver.CalculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	}
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
	for i := 0; i < len(seq); i += pgs.NumHaplotypes {
		if seq[i] == 1 && seq[i+1] == 1 {
			num++
		}
	}
	return num
}
