package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strings"
	"time"

	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

const (
	SCORE      = "score"
	ITERATIONS = 5000
	MARGIN     = 0.0000001
	numThreads = 16
)

func main() {
	//INDIVIDUAL := "NA20543"
	//INDIVIDUAL := "NA11881"
	INDIVIDUAL := "NA18595"
	//INDIVIDUAL := "HG03304"
	//INDIVIDUAL := "NA19082"
	//INDIVIDUAL := "HG03022"
	//INDIVIDUAL := "HG01767"
	//INDIVIDUAL := "HG01868"
	p := pgs.NewPGS()
	//err := p.LoadCatalogFile("PGS000073_hmPOS_GRCh38.txt")
	//err := p.LoadCatalogFile("PGS000037_hmPOS_GRCh38.txt")
	err := p.LoadCatalogFile("PGS000040_hmPOS_GRCh38.txt")
	if err != nil {
		log.Println("Error:", err)
		return
	}
	p.LoadPriors()
	cohort := NewCohort()
	cohort.CalculatePRS(p)

	err = cohort.SaveScores(p.PgsID)
	if err != nil {
		log.Println("Error saving scores:", err)
		return
	}
	target := make([]int, len(p.Variants))
	for i, locus := range p.Loci {
		target[i] = int(cohort[INDIVIDUAL][locus])
	}

	solutions := Solve(cohort[INDIVIDUAL][SCORE], p)
	//solutions = sortByAccuracy(solutions, target)
	fmt.Printf("True:\n%s -- %f\n", arrayTostring(target), p.CalculateSequenceLikelihood(target))
	fmt.Printf("Guessed:\n")
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.5f, %.2f\n", arrayTostring(solution), accuracy(solution, target),
			cohort[INDIVIDUAL][SCORE]-calculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	}
	//fmt.Printf("True score:%f", cohort[INDIVIDUAL][SCORE])
	//fmt.Printf("\nGuessed scores:%f\n", calculateScore(solution, p.Weights))
	//fmt.Printf("Accuracy: %.2f\n", accuracy(solution, target))
	//allPermutations(solution, p.Weights)
}

func Solve(targetScore float64, pgs *pgs.PGS) [][]int {
	var err error
	candidates := make([][]int, len(pgs.Variants)*4)
	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = pgs.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	// Evaluate candidates
	for k := 0; k < ITERATIONS; k++ {
		fitness := calculateFitness(candidates, pgs, targetScore, k)
		parents, children := crossover(candidates, fitness)
		childrenFitness := calculateFitness(children, pgs, targetScore, k)
		candidates = tournament(append(parents, children...), append(fitness, childrenFitness...))
		//	mutations happen independently of the score, they are based only on likelihood
		candidates = mutate(candidates, pgs)
	}
	candidates = sortByFitness(candidates, pgs, targetScore)
	return candidates
}

func crossover(parents [][]int, fitness []float64) ([][]int, [][]int) {
	rand.NewSource(time.Now().UnixNano())
	offspring := make([][]int, 0, len(parents))
	parents, fitness = tools.Shuffle(parents, fitness)
	cumulative := 0.0
	var splitPos, step float64
	splice := func(firstIndex, secondIndex int) {
		child := make([]int, len(parents[firstIndex]))
		cumulative = fitness[firstIndex] + fitness[secondIndex]
		splitPos = rand.Float64() * cumulative
		step = cumulative / float64(len(parents[firstIndex]))
		copy(child[:int(splitPos/step)], parents[firstIndex][:int(splitPos/step)])
		copy(child[int(splitPos/step):], parents[secondIndex][int(splitPos/step):])
		offspring = append(offspring, child)
	}
	// each parent has two offspring: one with the next parent, and one with the parent at + len(parents)/2
	for i := 0; i < len(parents); i += 2 {
		splice(i, i+1)
	}
	for i := 0; i < len(parents)/2; i++ {
		splice(i, i+len(parents)/2)
	}
	return parents, offspring
}

func tournament(population [][]int, fitness []float64) [][]int {
	survivors := make([][]int, 0, len(population)/2)
	population, fitness = tools.Shuffle(population, fitness)
	// Select the best half of the population
	for i := 0; i < len(population); i += 2 {
		if fitness[i] > fitness[i+1] {
			survivors = append(survivors, population[i])
		} else {
			survivors = append(survivors, population[i+1])
		}
	}
	return survivors
}

func mutate(population [][]int, pgs *pgs.PGS) [][]int {
	for j, individual := range population {
		population[j] = pgs.MutateGenome(individual)
	}
	return population
}

func calculateFitness(population [][]int, pgs *pgs.PGS, targetScore float64, iteration int) []float64 {
	var delta float64
	fitness := make([]float64, len(population))
	for i, individual := range population {
		delta = calculateScore(individual, pgs.Weights) - targetScore
		// smaller the delta, higher the fitness. the delta's significance increases with iterations
		fitness[i] = math.Exp(-math.Abs(delta) / float64(ITERATIONS-iteration))
	}
	return fitness
}

func calculateScore(snps []int, weights []float64) float64 {
	score := 0.0
	for i, snp := range snps {
		score += float64(snp) * weights[i]
	}
	return score
}

func sortByFitness(population [][]int, pgs *pgs.PGS, targetScore float64) [][]int {
	fitness := calculateFitness(population, pgs, targetScore, 0)
	for i := 0; i < len(population); i++ {
		for j := i; j < len(population); j++ {
			if fitness[i] < fitness[j] {
				population[i], population[j] = population[j], population[i]
				fitness[i], fitness[j] = fitness[j], fitness[i]
			}
		}
	}
	return population
}

func accuracy(solution []int, target []int) float64 {
	if len(solution) != len(target) {
		return 0.0
	}
	acc := 0.0
	for i := range solution {
		if solution[i] == target[i] {
			acc++
		}
	}
	return acc / float64(len(solution))
}

func arrayTostring(array []int) string {
	str := make([]string, len(array))
	for i, num := range array {
		str[i] = fmt.Sprint(num)
	}
	return strings.Join(str, "")
}

func sortByAccuracy(solutions [][]int, target []int) [][]int {
	accuracies := make([]float64, len(solutions))
	for i, solution := range solutions {
		accuracies[i] = accuracy(solution, target)
	}
	for i := 0; i < len(solutions)-1; i++ {
		for j := i + 1; j < len(solutions); j++ {
			if accuracies[i] < accuracies[j] {
				solutions[i], solutions[j] = solutions[j], solutions[i]
				accuracies[i], accuracies[j] = accuracies[j], accuracies[i]
			}
		}
	}
	return solutions
}
