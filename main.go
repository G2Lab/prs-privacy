package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"strings"
	"sync"
	"time"

	"github.com/nikirill/prs/tools"
)

const (
	SCORE      = "score"
	ITERATIONS = 1000
	MARGIN     = 0.0000001
	numThreads = 16
)

var GENOTYPES = []int{0, 1, 2}

func main() {
	//INDIVIDUAL := "NA20543"
	//INDIVIDUAL := "NA11881"
	INDIVIDUAL := "NA18595"
	//INDIVIDUAL := "HG03304"
	//INDIVIDUAL := "NA19082"
	//INDIVIDUAL := "HG03022"
	//INDIVIDUAL := "HG01767"
	//INDIVIDUAL := "HG01868"
	pgs := NewPGS()
	//err := pgs.LoadCatalogFile("PGS000073_hmPOS_GRCh38.txt")
	//err := pgs.LoadCatalogFile("PGS000037_hmPOS_GRCh38.txt")
	err := pgs.LoadCatalogFile("PGS000040_hmPOS_GRCh38.txt")
	if err != nil {
		log.Println("Error:", err)
		return
	}
	pgs.LoadPriors()
	cohort := NewCohort()
	cohort.CalculatePRS(pgs)

	err = cohort.SaveScores(pgs.PgsID)
	if err != nil {
		log.Println("Error saving scores:", err)
		return
	}
	target := make([]int, len(pgs.Variants))
	for i, locus := range pgs.Loci {
		target[i] = int(cohort[INDIVIDUAL][locus])
	}

	solutions := Solve(cohort[INDIVIDUAL][SCORE], pgs)
	//for i := range pgs.Weights {
	//	fmt.Printf("%.5f ", pgs.Weights[i])
	//}
	sortedSolutions := sortByAccuracy(solutions, target)
	fmt.Printf("True:\n%s -- %f\n", arrayTostring(target), pgs.CalculateSequenceLikelihood(target))
	fmt.Printf("Guessed:\n")
	for _, solution := range sortedSolutions {
		fmt.Printf("%s -- %f, %f\n", arrayTostring(solution), accuracy(solution, target), pgs.CalculateSequenceLikelihood(solution))
	}
	//fmt.Printf("True score:%f", cohort[INDIVIDUAL][SCORE])
	//fmt.Printf("\nGuessed scores:%f\n", calculateScore(solution, pgs.Weights))
	//fmt.Printf("Accuracy: %.2f\n", accuracy(solution, target))
	//allPermutations(solution, pgs.Weights)
}

func Solve(targetScore float64, pgs *PGS) map[string][]int {
	var err error
	candidates := make([][]int, (len(pgs.Variants)/numThreads)*2*numThreads)
	candidatesPerThread := len(candidates) / numThreads
	deltas := make([]float64, len(candidates))
	weights := pgs.Weights

	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = pgs.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	var workerPool sync.WaitGroup
	crossOverDone := make([]chan struct{}, numThreads)
	collectChans := make([]chan map[string][]int, numThreads)
	for thread := 0; thread < numThreads; thread++ {
		collectChans[thread] = make(chan map[string][]int)
		go func(t int, wg sync.WaitGroup) {
			threadSolutions := make(map[string][]int, 0)
			// Evaluate candidates
			for k := 0; k < ITERATIONS; k++ {
				wg.Add(1)
				// Check for valid solutions
				for j, candidate := range candidates[t*candidatesPerThread : (t+1)*candidatesPerThread] {
					deltas[t*candidatesPerThread+j] = calculateScore(candidate, weights) - targetScore
					// Found a solution
					if math.Abs(deltas[t*candidatesPerThread+j]) < MARGIN {
						threadSolutions[arrayTostring(candidate)] = candidate
						candidates[t*candidatesPerThread+j], err = pgs.SampleFromPopulation()
						if err != nil {
							fmt.Println(err)
							wg.Done()
							return
						}
						deltas[t*candidatesPerThread+j] = calculateScore(candidate, weights) - targetScore
					}
				}
				wg.Done()
				<-crossOverDone[t]
			candidateLoop:
				for j, candidate := range candidates[t*candidatesPerThread : (t+1)*candidatesPerThread] {
					//the weights that cover the delta better, get higher probability of being selected
					// big delta -> bigger weights get higher probability.
					// probability = 1 / abs( delta - weight * snp_old + weight * snp_new )
					possibleMutations := make([]int, 0, len(candidate)*(len(GENOTYPES)-1))
					probs := make([]float64, 0)
					for i, snp := range candidate {
						tmp := deltas[t*candidatesPerThread+j]
						tmp -= weights[i] * float64(snp)
						for _, v := range GENOTYPES {
							if v == snp {
								continue
							}
							possibleMutations = append(possibleMutations, v)
							if math.Abs(tmp+weights[i]*float64(v)) < MARGIN {
								candidate[i] = v
								threadSolutions[arrayTostring(candidate)] = candidate
								candidates[t*candidatesPerThread+j], err = pgs.SampleFromPopulation()
								if err != nil {
									fmt.Println(err)
									return
								}
								continue candidateLoop
							}
							probs = append(probs, 1/math.Abs(tmp+weights[i]*float64(v)))
						}
					}
					candidates[t*candidatesPerThread+j] = pgs.MutateGenome(candidate, possibleMutations, probs)
				}
			}
			collectChans[t] <- threadSolutions
		}(thread, workerPool)
	}
	// Cross over among all candidates
	for k := 0; k < ITERATIONS; k++ {
		workerPool.Wait()
		parents, offspring := crossover(candidates, deltas)
		offspringDeltas := make([]float64, len(offspring))
		for i, child := range offspring {
			offspringDeltas[i] = calculateScore(child, weights) - targetScore
		}
		candidates, deltas = tournament(append(parents, offspring...), append(deltas, offspringDeltas...))
		for thread := 0; thread < numThreads; thread++ {
			crossOverDone[thread] <- struct{}{}
		}
	}
	// Collect solutions from all the threads
	solutions := make(map[string][]int)
	for thread := 0; thread < numThreads; thread++ {
		for k, v := range <-collectChans[thread] {
			solutions[k] = v
		}
	}
	return solutions
}

func calculateScore(snps []int, weights []float64) float64 {
	score := 0.0
	for i, snp := range snps {
		score += float64(snp) * weights[i]
	}
	return score
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

func sortByAccuracy(solutions map[string][]int, target []int) [][]int {
	sorted := make([][]int, 0, len(solutions))
	for _, solution := range solutions {
		sorted = append(sorted, solution)
	}
	for i := 0; i < len(sorted)-1; i++ {
		maxIndex := i
		maxAccuracy := accuracy(sorted[i], target)
		for j := i + 1; j < len(sorted); j++ {
			if accuracy(sorted[j], target) > maxAccuracy {
				maxIndex = j
				maxAccuracy = accuracy(sorted[j], target)
			}
		}
		if maxIndex != i {
			sorted[i], sorted[maxIndex] = sorted[maxIndex], sorted[i]
		}
	}
	return sorted
}

func crossover(parents [][]int, fitness []float64) ([][]int, [][]int) {
	rand.NewSource(time.Now().UnixNano())
	offspring := make([][]int, 0, len(parents))
	parents, fitness = tools.Shuffle(parents, fitness)
	cumulative := 0.0
	var splitPos, step float64
	splice := func(firstIndex, secondIndex int) {
		child := make([]int, len(parents[firstIndex]))
		cumulative = 1/fitness[firstIndex] + 1/fitness[secondIndex]
		splitPos = rand.Float64() * cumulative
		step = cumulative / float64(len(parents[firstIndex]))
		fmt.Println(cumulative, splitPos, step)
		copy(child[:int(splitPos/step)], parents[firstIndex][:int(splitPos/step)])
		copy(child[int(splitPos/step):], parents[secondIndex][int(splitPos/step):])
		offspring = append(offspring, child)
	}
	fmt.Println(fitness)
	// each parent has two offspring: one with the next parent, and one with the parent at + len(parents)/2
	for i := 0; i < len(parents); i += 2 {
		fmt.Println(i, i+1)
		splice(i, i+1)
	}
	for i := 0; i < len(parents)/2; i++ {
		splice(i, i+len(parents)/2)
	}
	return parents, offspring
}

func tournament(population [][]int, fitness []float64) ([][]int, []float64) {
	survivors := make([][]int, 0, len(population)/2)
	survFit := make([]float64, 0, len(fitness)/2)
	population, fitness = tools.Shuffle(population, fitness)
	// Select the best half of the population
	for i := 0; i < len(population); i += 2 {
		if fitness[i] > fitness[i+1] {
			survivors = append(survivors, population[i])
			survFit = append(survFit, fitness[i])
		} else {
			survivors = append(survivors, population[i+1])
			survFit = append(survFit, fitness[i+1])
		}
	}
	return survivors, survFit
}

func allPermutations(origin []int, weights []float64) [][]int {
	permutations := make([][]int, 1)
	permutations[0] = origin
	// Find all the positions by their weights
	duplicates := make(map[float64][]int)
	for i, weight := range weights {
		if _, ok := duplicates[weight]; !ok {
			duplicates[weight] = make([]int, 0)
		}
		duplicates[weight] = append(duplicates[weight], i)
	}
	// Remove all the unique weights, and leave only the duplicates
	for _, weight := range weights {
		if len(duplicates[weight]) < 2 {
			delete(duplicates, weight)
		}
	}
	//for i := range origin {
	//
	//}
	return permutations
}
