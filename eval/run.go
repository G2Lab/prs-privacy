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
	likelihoodEffect()
	//samples()
	//distribution()
	//findAllSolutions()
}

func likelihoodEffect() {
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
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	//samples := cohort.SortByScore()
	samples := allSamples()
	//
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
		chunkSize = len(samples)
	}
	// Create csv result file
	filename := fmt.Sprintf("%s-%d.csv", p.PgsID, chunkNum)
	resFile, err := os.OpenFile(filename, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	writer := csv.NewWriter(resFile)
	defer writer.Flush()
	writer.Write([]string{"individual", "first accuracy", "true position", "total solutions", "first likelihood",
		"true likelihood", "major likelihood"})
	ctx, cancel := context.WithCancel(context.Background())
	numThreads := getNumThreads()
	majorReference := p.AllMajorSample()
	majorLikelihood := p.CalculateSequenceLikelihood(majorReference)
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
				writer.Write([]string{samples[i], fmt.Sprintf("%.3f", firstAcc), fmt.Sprintf("%d", j),
					fmt.Sprintf("%d", len(solutions)), fmt.Sprintf("%.3f", firstLikelihood),
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
	//INDIVIDUAL := "HG00551" // low 648
	//
	INDIVIDUAL := "HG00266"
	//INDIVIDUAL := "NA18613"

	p := pgs.NewPGS()
	catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
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

	solmap := slv.DP(numThreads)
	//solmap := slv.solve(numThreads)
	//solmap := slv.recursive(numThreads)
	cancel()
	//solmap = findComplements(solmap, p, numThreads)
	solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f, %f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
		cohort[INDIVIDUAL].Score-solver.CalculateScore(cohort[INDIVIDUAL].Genotype, p.Weights),
		p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
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
	fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solutions[0]),
		solver.Accuracy(solutions[0], cohort[INDIVIDUAL].Genotype),
		cohort[INDIVIDUAL].Score-solver.CalculateScore(solutions[0], p.Weights), p.CalculateSequenceLikelihood(solutions[0]))
	//for _, solution := range solutions {
	//	fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
	//		cohort[INDIVIDUAL].Score-solver.CalculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	//	//fmt.Printf("%s -- %.3f, %.8f\n", arrayToStringDiploid(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
	//	//	cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights))
	//}
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
