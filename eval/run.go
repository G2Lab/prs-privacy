package main

import (
	"context"
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
	samples()
	//distribution()
	//findAllSolutions()
}

func samples() {
	smp := selectSamples(10, 2504, "prs")
	fmt.Println(smp)
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
	step := 50
	numThreads := getNumThreads()
	fullMajor := make([]uint8, 2*len(p.Weights))
	for i := 0; i < len(p.Weights); i++ {
		fmt.Printf("%.4f ", p.Maf[i][0])
	}
	fmt.Println()
	for i := 0; i < len(p.Weights); i++ {
		if p.Maf[i][0] > 0.5 {
			fullMajor[2*i] = 0
			fullMajor[2*i+1] = 0
		} else {
			fullMajor[2*i] = 1
			fullMajor[2*i+1] = 1
		}
	}
	fmt.Printf("Full major %s: score %f, likelihood %f\n", solver.ArrayToString(fullMajor), solver.CalculateScore(fullMajor, p.Weights), p.CalculateSequenceLikelihood(fullMajor))
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
		htrz := 0
		snps := 0
		for j, solution := range solutions {
			if j == 0 {
				acc = solver.Accuracy(solution, cohort[sorted[i]].Genotype)
				likelihood = p.CalculateSequenceLikelihood(solution)
				htrz = numHeterozygous(solution)
				snps = numSNPs(solution)
			}
			if solver.Accuracy(solution, cohort[sorted[i]].Genotype) == 1.0 {
				fmt.Printf("%s, score %f -- likelihood %d/%d, first acc %.3f, first / true likelihood %f/%f, "+
					"snps %d/%d, heterozygous %d/%d\n",
					sorted[i], cohort[sorted[i]].Score, j, len(solmap), acc,
					likelihood, p.CalculateSequenceLikelihood(cohort[sorted[i]].Genotype),
					snps, numSNPs(cohort[sorted[i]].Genotype), htrz, numHeterozygous(cohort[sorted[i]].Genotype))
				break
			}
		}
	}
	cancel()
}

func findAllSolutions() {
	//INDIVIDUAL := "NA18595"
	INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	//INDIVIDUAL := "HG00551" // low 648

	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	catalogFile := "PGS000639_hmPOS_GRCh38.txt"
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
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-solver.CalculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
		//fmt.Printf("%s -- %.3f, %.8f\n", arrayToStringDiploid(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
		//	cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights))
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
		log.Printf("Could not read numCpus env %s\n", params.NumCpusEnv)
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
