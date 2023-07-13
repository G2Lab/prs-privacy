package main

import (
	"fmt"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"math/rand"
	"strings"
	"time"
)

const (
	ITERATIONS = 10000
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
	err := p.LoadCatalogFile("PGS000073_hmPOS_GRCh38.txt")
	//err := p.LoadCatalogFile("PGS000037_hmPOS_GRCh38.txt")
	//err := p.LoadCatalogFile("PGS000040_hmPOS_GRCh38.txt")
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	cohort := NewCohort(p)

	solutions := Solve(cohort[INDIVIDUAL].Score, p)
	solutions = sortByAccuracy(solutions, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f\n", arrayTostring(cohort[INDIVIDUAL].Genotype), p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Guessed %d:\n", len(solutions))
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.5f, %.2f\n", arrayTostring(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	}
}

func Solve(targetScore float64, p *pgs.PGS) [][]int {
	var err error
	rand.NewSource(time.Now().UnixNano())
	POOLSIZE := len(p.Variants) * 10
	solutions := make([][]int, 0)
	candidates := make([][]int, POOLSIZE, 2*POOLSIZE)
	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = p.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	fmt.Printf("Iteration: ")
	// Evaluate candidates
	for k := 0; k < ITERATIONS; k++ {
		if k%1000 == 0 {
			fmt.Printf("%d/%d ", k, ITERATIONS)
		}
		deltas, solved := calculateDeltas(candidates, p, targetScore)
		if len(solved) > 0 {
			placeholder := make([][]int, len(solved))
			solutions = append(solutions, placeholder...)
			copy(solutions[len(solutions)-len(solved):], solved)
		}
		//solutions = append(solutions, solved...)
		children := crossover(candidates, p, targetScore)
		//_, solved = calculateDeltas(children, p, targetScore)
		chDeltas, solved := calculateDeltas(children, p, targetScore)
		if len(solved) > 0 {
			placeholder := make([][]int, len(solved))
			solutions = append(solutions, placeholder...)
			copy(solutions[len(solutions)-len(solved):], solved)
		}
		////solutions = append(solutions, solved...)
		candidates = tournament(append(candidates, children...), append(deltas, chDeltas...), POOLSIZE)
		deltas, solved = calculateDeltas(candidates, p, targetScore)
		if len(solved) > 0 {
			placeholder := make([][]int, len(solved))
			solutions = append(solutions, placeholder...)
			copy(solutions[len(solutions)-len(solved):], solved)
		}
		//solutions = append(solutions, solved...)
		//	mutations happen independently of the score, they are based only on likelihood
		mutate(candidates, deltas, p)
	}
	//solutions = sortByFitness(solutions, p, targetScore)
	return solutions
}

func calculateDeltas(population [][]int, p *pgs.PGS, targetScore float64) ([]float64, [][]int) {
	var delta float64
	var err error
	deltas := make([]float64, len(population))
	matches := make([][]int, 0)
	for i := range population {
		delta = calculateScore(population[i], p.Weights) - targetScore
		for delta == 0 {
			match := make([]int, len(population[i]))
			copy(match, population[i])
			matches = append(matches, match)
			population[i], err = p.SampleFromPopulation()
			if err != nil {
				log.Printf("Error resampling in fitness calculation: %v\n", err)
				return nil, nil
			}
			delta = calculateScore(population[i], p.Weights) - targetScore
		}
		deltas[i] = delta
	}
	return deltas, matches
}

func calculateScore(snps []int, weights []float64) float64 {
	score := 0.0
	for i := 0; i < len(snps); i += pgs.NUM_HAPLOTYPES {
		for j := 0; j < pgs.NUM_HAPLOTYPES; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				score += weights[i/2]
			default:
				log.Printf("Invalid alelle value: %d", snps[i+j])
			}
		}
	}
	return score
}

func crossover(population [][]int, p *pgs.PGS, target float64) [][]int {
	parents := tools.Shuffle(population)
	splice := func(first, second int) [][]int {
		children := make([][]int, len(parents[first])/2)
		for k := 0; k < len(children); k++ {
			children[k] = make([]int, len(parents[first]))
			copy(children[k][:k*pgs.NUM_HAPLOTYPES], parents[first][:k*pgs.NUM_HAPLOTYPES])
			copy(children[k][k*pgs.NUM_HAPLOTYPES:], parents[second][k*pgs.NUM_HAPLOTYPES:])
		}
		deltas, matches := calculateDeltas(children, p, target)
		if len(matches) > 0 {
			return matches
		}
		fitness := deltasToFitness(deltas)
		pickedIndex := tools.SampleFromDistribution(fitness)
		return [][]int{children[pickedIndex]}
	}
	offspring := make([][]int, 0, len(parents))
	// each parent has two offspring: one with the next parent, and one with the parent at + len(parents)/2
	for i := 0; i < len(parents); i += 2 {
		offspring = append(offspring, splice(i, i+1)...)
	}
	//for i := 0; i < len(parents)/2; i++ {
	//	offspring = append(offspring, splice(i, i+len(parents)/2))
	//}
	return offspring
}

func tournament(population [][]int, deltas []float64, populationSize int) [][]int {
	fitness := deltasToFitness(deltas)
	scores := make([]float64, len(population))
	var opponentIdx int
	for i := range population {
		for k := 0; k < len(population)/10; k++ {
			opponentIdx = rand.Intn(len(population))
			for opponentIdx == i {
				opponentIdx = rand.Intn(len(population))
			}
			if fitness[i] >= fitness[opponentIdx] {
				scores[i]++
			}
		}
	}
	sortBy(population, scores)
	return population[:populationSize]
}

//func tournament(population [][][]int, deltas []float64) [][][]int {
//	fitness := deltasToFitness(deltas)
//	survivors := make([][][]int, 0, len(population)/2)
//	population, fitness = tools.ShuffleWithLabels(population, fitness)
//	// Select the best half of the population
//	for i := 0; i < len(population); i += 2 {
//		if fitness[i] > fitness[i+1] {
//			survivors = append(survivors, population[i])
//		} else {
//			survivors = append(survivors, population[i+1])
//		}
//	}
//	return survivors
//}

func mutate(population [][]int, deltas []float64, p *pgs.PGS) {
	for j, original := range population {
		population[j] = p.MutateGenome(original, deltas[j])
	}
}

func deltasToFitness(deltas []float64) []float64 {
	fitness := make([]float64, len(deltas))
	for i, delta := range deltas {
		// smaller the delta, higher the fitness
		fitness[i] = 1 / math.Abs(delta)
	}
	return fitness
}

func accuracy(solution []int, target []int) float64 {
	if len(solution) != len(target) {
		return 0.0
	}
	acc := 0.0
	for i := 0; i < len(solution); i += pgs.NUM_HAPLOTYPES {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.NUM_HAPLOTYPES / float64(len(solution))
}

func arrayTostring(array []int) string {
	str := make([]string, len(array))
	for i := 0; i < len(array); i += pgs.NUM_HAPLOTYPES {
		str[i] = fmt.Sprint(array[i] + array[i+1])
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

func sortBy(items [][]int, properties []float64) {
	for i := 0; i < len(items)-1; i++ {
		for j := i + 1; j < len(items); j++ {
			if properties[i] < properties[j] {
				items[i], items[j] = items[j], items[i]
				properties[i], properties[j] = properties[j], properties[i]
			}
		}
	}
}
