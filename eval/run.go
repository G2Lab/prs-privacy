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
	"strconv"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/solver"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
var memprofile = flag.String("memprofile", "", "write cpu profile to file")

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

	//scoreDistribution()
	findAllSolutions()
	//likelihoodEffect()
	//scoreToLikelihoodDistribution()
	//scoreToLikelihood()
	//accuracyLikelihood()
	//samples()
	//distribution()
	//evaluateReferences()
	//evaluateGA()
}

type Result struct {
	Individual  string
	Score       string
	Accuracies  []string
	Likelihoods []string
}

func NewResult(ind string, score *apd.Decimal) *Result {
	return &Result{
		Individual:  ind,
		Score:       score.String(),
		Accuracies:  make([]string, 0),
		Likelihoods: make([]string, 0),
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
//	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
//	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
//	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
//	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
//	if err != nil {
//		log.Printf("Error loading catalog file: %v\n", err)
//		return
//	}
//	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
//	p.LoadStats()
//	cohort := solver.NewCohort(p)
//
//	slv := solver.NewGenetic(cohort[INDIVIDUAL].Score, cohort[INDIVIDUAL].Genotype, p)
//
//	solmap := slv.Solve(1)
//	solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
//	fmt.Printf("\nTrue:\n%s -- %f, %f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
//		cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(cohort[INDIVIDUAL].Genotype, p.Weights),
//		p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
//	fmt.Printf("Major:\n%s -- %f, %f\n", solver.ArrayToString(p.AllMajorSample()),
//		cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(p.AllMajorSample(), p.Weights),
//		p.CalculateSequenceLikelihood(p.AllMajorSample()))
//	fmt.Printf("Guessed %d:\n", len(solutions))
//	//fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solutions[0]),
//	//	solver.Accuracy(solutions[0], cohort[INDIVIDUAL].Genotype),
//	//	cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(solutions[0], p.Weights), p.CalculateSequenceLikelihood(solutions[0]))
//	for _, solution := range solutions {
//		fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
//			cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
//	}
//}

func evaluateReferences() {
	resFolder := "results"
	catalogFiles := []string{
		"PGS002302_hmPOS_GRCh38.txt",
		"PGS000639_hmPOS_GRCh38.txt",
		"PGS000037_hmPOS_GRCh38.txt",
		"PGS000073_hmPOS_GRCh38.txt",
		"PGS001827_hmPOS_GRCh38.txt"}
	// Create a result file
	filepath := path.Join(resFolder, "references.json")
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	samples := allSamples()
	references := make([]*Reference, 0, len(catalogFiles))
	for _, catalogFile := range catalogFiles {
		p := pgs.NewPGS()
		err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
		if err != nil {
			log.Printf("Error loading catalog file: %v\n", err)
			return
		}
		p.LoadStats()
		cohort := solver.NewCohort(p)
		fmt.Printf("%s\n", p.PgsID)
		major := p.AllMajorSample()
		reference := NewReference(p.PgsID, solver.CalculateDecimalScore(p.Context, major, p.Weights), p.CalculateSequenceLikelihood(major))
		for _, sample := range samples {
			reference.Accuracies = append(reference.Accuracies, fmt.Sprintf("%.3f", solver.Accuracy(major, cohort[sample].Genotype)))
		}
		references = append(references, reference)
	}
	if err = encoder.Encode(references); err != nil {
		log.Fatal("Cannot encode json", err)
	}
}

func accuracyLikelihood() {
	resFolder := "results/accuracyLikelihood"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	fmt.Printf("%s\n", p.PgsID)
	cohort := solver.NewCohort(p)
	//samples := cohort.SortByScore()[len(cohort)-100:]
	samples := allSamples()
	// Create csv result file
	chunkNum, chunkSize := getChunkInfo(len(samples))
	filepath := path.Join(resFolder, fmt.Sprintf("%s-%d.json", p.PgsID, chunkNum))
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)
	var solved bool
	var acc float64
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d\n", i)
		slv := solver.NewDP(cohort[samples[i]].Score, p)
		solmap := slv.Solve()
		solutions := solver.SortByLikelihood(solmap, p)
		result := NewResult(samples[i], cohort[samples[i]].Score)
		solved = false
		for _, solution := range solutions {
			acc = solver.Accuracy(solution, cohort[samples[i]].Genotype)
			if acc == 1.0 {
				solved = true
			}
			result.Accuracies = append(result.Accuracies, fmt.Sprintf("%.3f", acc))
			result.Likelihoods = append(result.Likelihoods, fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)))
		}
		if len(solutions) == 0 || !solved {
			fmt.Printf("No solution for %s: %v\n", samples[i], cohort[samples[i]].Genotype)
			continue
		}
		if err = encoder.Encode(result); err != nil {
			log.Fatal("Cannot encode json", err)
		}
	}
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
	fmt.Printf("%s\n", p.PgsID)
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
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d ", i)
		output := []string{samples[i], fmt.Sprintf("%.3f", cohort[samples[i]].Score)}
		slv := solver.NewDP(cohort[samples[i]].Score, p)
		solmap := slv.Solve()
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
//		solmap := slv.Solve(numThreads)
//		solutions := solver.SortByAccuracy(solmap, cohort[sample].Genotype)
//		if len(solutions) == 0 || solver.Accuracy(solutions[0], cohort[sample].Genotype) != 1.0 {
//			fmt.Printf("No solution for %s: %s, %.3f\n", sample, solver.ArrayToString(solutions[0]),
//				solver.Accuracy(solutions[0], cohort[sample].Genotype))
//			continue
//		}
//		for _, solution := range solutions {
//			output = append(output, fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)))
//		}
//		writer.Write(output)
//	}
//	cancel()
//}

func likelihoodEffect() {
	resFolder := "results"
	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
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
	majorReference := p.AllMajorSample()
	majorScore := solver.CalculateDecimalScore(p.Context, majorReference, p.Weights)
	majorLikelihood := p.CalculateSequenceLikelihood(majorReference)
	medianScore := cohort[cohort.SortByScore()[len(cohort)/2]].Score
	for i := chunkNum * chunkSize; i < (chunkNum+1)*chunkSize; i++ {
		if i >= len(samples) {
			break
		}
		fmt.Printf("%d ", i)
		slv := solver.NewDP(cohort[samples[i]].Score, p)
		solmap := slv.Solve()
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
				writer.Write([]string{samples[i], fmt.Sprintf("%.3f", solver.CalculateDecimalScore(p.Context, solution, p.Weights)),
					fmt.Sprintf("%.3f", majorScore), fmt.Sprintf("%.3f", medianScore),
					fmt.Sprintf("%.3f", firstAcc),
					fmt.Sprintf("%.3f", solver.Accuracy(majorReference, cohort[samples[i]].Genotype)),
					fmt.Sprintf("%d", j), fmt.Sprintf("%d", len(solutions)), fmt.Sprintf("%.3f", firstLikelihood),
					fmt.Sprintf("%.3f", p.CalculateSequenceLikelihood(solution)), fmt.Sprintf("%.3f", majorLikelihood)})
				break solLoop
			}
		}
		if acc != 1.0 {
			diff := new(apd.Decimal)
			p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalScore(p.Context, cohort[samples[i]].Genotype, p.Weights))
			fmt.Printf("\nNo solution for %s\n", samples[i])
			fmt.Printf("\nTrue:\n%s -- %s, %f\n", solver.ArrayToString(cohort[samples[i]].Genotype),
				diff.String(),
				p.CalculateSequenceLikelihood(cohort[samples[i]].Genotype))
			for _, sol := range solutions {
				p.Context.Sub(diff, cohort[samples[i]].Score, solver.CalculateDecimalScore(p.Context, sol, p.Weights))
				fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(sol), solver.Accuracy(sol, cohort[samples[i]].Genotype),
					diff.String(), p.CalculateSequenceLikelihood(sol))
			}
		}
	}
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
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	cohort := solver.NewCohort(p)
	sorted := cohort.SortByScore()
	fmt.Printf("Median score: %g\n", cohort[sorted[len(sorted)/2]].Score)
	step := 250
	majorReference := p.AllMajorSample()
	fmt.Printf("Full major %s: score %f, likelihood %f\n", solver.ArrayToString(majorReference), solver.CalculateDecimalScore(p.Context, majorReference, p.Weights), p.CalculateSequenceLikelihood(majorReference))
	for i := 0; i < len(sorted); i += step {
		slv := solver.NewDP(cohort[sorted[i]].Score, p)
		solmap := slv.Solve()
		solutions := solver.SortByLikelihood(solmap, p)
		//fmt.Printf("\n\n%s, %s\n", p.PgsID, sorted[i])
		//for _, solution := range solutions {
		//	fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[sorted[i]].Genotype),
		//		cohort[sorted[i]].Score-solver.CalculateDecimalScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
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
}

func findAllSolutions() {
	//INDIVIDUAL := "NA18595"
	//INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	INDIVIDUAL := "HG00551" // low 648
	//INDIVIDUAL := "NA12286"
	//
	//INDIVIDUAL := "HG01028"
	//INDIVIDUAL := "NA18531"
	//INDIVIDUAL := "NA20872"
	//INDIVIDUAL := "NA20507"
	//INDIVIDUAL := "HG00242"

	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000845_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000534_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	p.LoadStats()
	cohort := solver.NewCohort(p)

	slv := solver.NewDP(cohort[INDIVIDUAL].Score, p)

	//majorReference := p.AllMajorSample()
	//fmt.Printf("Accuracy with major: %f\n",
	//	solver.Accuracy(majorReference, cohort[INDIVIDUAL].Genotype))

	solmap := slv.Solve()
	solutions := solver.SortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	//solutions := solver.SortByLikelihood(solmap, p)
	//fmt.Printf("\nTrue:\n%s -- %f, %f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype),
	//	cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(cohort[INDIVIDUAL].Genotype, p.Weights),
	//	p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	//fmt.Printf("Major likelihood: %f\n", p.CalculateSequenceLikelihood(majorReference))

	fmt.Printf("\nTrue:\n%s, %.2f\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype), p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Guessed %d:\n", len(solutions))
	//fmt.Printf("%s -- %.3f, %.12f, %.2f\n", solver.ArrayToString(solutions[0]),
	//	solver.Accuracy(solutions[0], cohort[INDIVIDUAL].Genotype),
	//	cohort[INDIVIDUAL].Score-solver.CalculateDecimalScore(solutions[0], p.Weights), p.CalculateSequenceLikelihood(solutions[0]))
	for _, solution := range solutions {
		diff := new(apd.Decimal)
		p.Context.Sub(diff, cohort[INDIVIDUAL].Score, solver.CalculateDecimalScore(p.Context, solution, p.Weights))
		fmt.Printf("%s -- %.3f, %s, %.2f\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype),
			diff.String(), p.CalculateSequenceLikelihood(solution))
	}
	//fmt.Printf("\nTrue:\n%s\n", solver.ArrayToString(cohort[INDIVIDUAL].Genotype))
	//fmt.Printf("Guessed %d:\n", len(solutions))
	//for _, solution := range solutions {
	//	fmt.Printf("%s (acc %.3f)\n", solver.ArrayToString(solution), solver.Accuracy(solution, cohort[INDIVIDUAL].Genotype))
	//}
}

func scoreDistribution() {
	p := pgs.NewPGS()
	catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000639_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000891_hmPOS_GRCh38.txt"
	//catalogFile := "PGS001827_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000066_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(path.Join(params.DataFolder, catalogFile))
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	cohort := solver.NewCohort(p)
	samples := allSamples()
	indices := []int{0, len(p.Weights) / 2, len(p.Weights)}

	resFolder := "results"
	filepath := path.Join(resFolder, "distro.json")
	resFile, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		log.Fatalf("Error opening result file: %v", err)
	}
	defer resFile.Close()
	encoder := json.NewEncoder(resFile)

	multiplier := apd.New(1, int32(p.WeightPrecision))
	if p.WeightPrecision > params.PrecisionsLimit {
		multiplier.SetFinite(1, params.PrecisionsLimit)
	}
	weights := solver.DecimalsToInts(p.Context, p.Weights, multiplier)
	maxTotalPositive, maxTotalNegative := solver.GetMaxTotal(weights)
	mf, _ := multiplier.Float64()
	maxTotalPositiveF, maxTotalNegativeF := float64(maxTotalPositive)/mf, float64(maxTotalNegative)/mf
	if err = encoder.Encode([]float64{maxTotalPositiveF, maxTotalNegativeF}); err != nil {
		log.Fatal("Cannot encode json for major-minor", err)
	}

	betas := make(map[uint16]int64, len(p.Weights))
	for i := 0; i < len(weights); i++ {
		betas[uint16(i)] = weights[i]
	}
	sampledMax, sampledMin := solver.SampleMaxMinScores(0, len(p.Weights), 100*len(p.Weights), betas, p)
	sampledMaxF, sampledMinF := float64(sampledMax)/mf, float64(sampledMin)/mf

	if err = encoder.Encode([]float64{sampledMaxF, sampledMinF}); err != nil {
		log.Fatal("Cannot encode json for major-minor", err)
	}

	fmt.Println(sampledMaxF, sampledMinF)
	scores := make([][]float64, 2)
	for r := 0; r < len(scores); r++ {
		scores[r] = make([]float64, 0, len(samples))
	}
	for _, sample := range samples {
		for i := 0; i < len(indices)-1; i++ {
			score, _ := solver.CalculateDecimalScore(p.Context, cohort[sample].Genotype[indices[i]*pgs.NumHplt:indices[i+1]*pgs.NumHplt],
				p.Weights[indices[i]:indices[i+1]]).Float64()
			scores[i] = append(scores[i], score)
		}
		//fullscore, _ := solver.CalculateDecimalScore(p.Context, cohort[sample].Genotype, p.Weights).Float64()
		//fmt.Println(sample, fullscore, scores[0][len(scores[0])-1], scores[1][len(scores[1])-1])
	}
	for r := 0; r < len(scores); r++ {
		if err = encoder.Encode(scores[r]); err != nil {
			log.Fatal("Cannot encode scores json", err)
		}
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
	for i := 0; i < len(seq); i += pgs.NumHplt {
		if seq[i] == 1 && seq[i+1] == 1 {
			num++
		}
	}
	return num
}
