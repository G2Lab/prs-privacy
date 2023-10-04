package main

import (
	"context"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"

	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
)

func main() {
	findAllSolutions()
}

func findAllSolutions() {
	//INDIVIDUAL := "NA18595"
	INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	//INDIVIDUAL := "HG00551" // low 648

	p := pgs.NewPGS()
	catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
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

	numThreadsEnv := os.Getenv(params.NumCpusEnv)
	var numThreads int
	if numThreadsEnv != "" {
		numThreads, err = strconv.Atoi(numThreadsEnv)
		if err != nil {
			log.Printf("Error parsing numCpus %s: %v\n", params.NumCpusEnv, err)
			return
		}
	} else {
		log.Printf("Could not read numCpus env %s: %v\n", params.NumCpusEnv, err)
		numThreads = 1
	}

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
