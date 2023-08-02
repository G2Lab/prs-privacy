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
	numThreads = 8
)

func main() {
	//INDIVIDUAL := "NA18595"
	INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040

	p := pgs.NewPGS()
	//err := p.LoadCatalogFile("PGS000073_hmPOS_GRCh38.txt")
	//err := p.LoadCatalogFile("PGS000037_hmPOS_GRCh38.txt")
	//err := p.LoadCatalogFile("PGS000040_hmPOS_GRCh38.txt")
	err := p.LoadCatalogFile("PGS000648_hmPOS_GRCh38.txt")
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	p.LoadStats()
	cohort := NewCohort(p)

	solmap := Solve(cohort[INDIVIDUAL].Score, p)
	solmap = findComplements(solmap, p)
	solutions := sortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f\n", arrayTostring(cohort[INDIVIDUAL].Genotype), p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Guessed %d:\n", len(solutions))
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.5f, %.2f\n", arrayTostring(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	}
}

func Solve(targetScore float64, p *pgs.PGS) map[string][]int {
	var err error
	rand.NewSource(time.Now().UnixNano())
	POOLSIZE := len(p.Variants) * 100
	solutions := make(map[string][]int)
	candidates := make([][]int, POOLSIZE, 2*POOLSIZE)
	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = p.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	//done := make([]chan struct{}, numThreads)
	//for i := 0; i < numThreads; i++ {
	//	done[i] = make(chan struct{})
	//}
	//Evaluate candidates
	for k := 0; k < ITERATIONS; k++ {
		if k%500 == 0 {
			fmt.Printf("Iteration %d/%d\n", k, ITERATIONS)
		}
		deltas, solved := calculateDeltas(candidates, p, targetScore)
		extend(solutions, solved)
		children := crossover(candidates, p, targetScore)
		chDeltas, solved := calculateDeltas(children, p, targetScore)
		extend(solutions, solved)
		candidates = tournament(append(candidates, children...), append(deltas, chDeltas...), POOLSIZE)
		deltas, solved = calculateDeltas(candidates, p, targetScore)
		extend(solutions, solved)
		candidates = mutate(candidates, deltas, p)
		//for thread := 0; thread < numThreads; thread++ {
		//	go func(t int) {
		//		mutate(candidates[t*len(candidates)/numThreads:(t+1)*len(candidates)/numThreads], deltas[t*len(candidates)/numThreads:(t+1)*len(candidates)/numThreads], p)
		//		done[t] <- struct{}{}
		//	}(thread)
		//}
		//for thread := 0; thread < numThreads; thread++ {
		//	<-done[thread]
		//}
	}

	//collectChans := make([]chan [][]int, numThreads)
	//for thread := 0; thread < numThreads; thread++ {
	//	collectChans[thread] = make(chan [][]int)
	//	go func(t int) {
	//		localSol := make([][]int, 0)
	//		cnd := make([][]int, len(candidates)/numThreads)
	//		copy(cnd, candidates[t*len(candidates)/numThreads:(t+1)*len(candidates)/numThreads])
	//		for k := 0; k < ITERATIONS; k++ {
	//			if k%1000 == 0 {
	//				fmt.Printf("Thread %d: %d/%d\n", t, k, ITERATIONS)
	//			}
	//			deltas, solved := calculateDeltas(cnd, p, targetScore)
	//			localSol = append(localSol, solved...)
	//			cnd = mutate(cnd, deltas, p)
	//		}
	//		collectChans[t] <- localSol
	//	}(thread)
	//}
	//for thread := 0; thread < numThreads; thread++ {
	//	collect := <-collectChans[thread]
	//	extend(solutions, collect)
	//}

	return solutions
}

func extend(base map[string][]int, extension [][]int) {
	for _, ext := range extension {
		if _, ok := base[arrayTostring(ext)]; !ok {
			base[arrayTostring(ext)] = make([]int, len(ext))
			copy(base[arrayTostring(ext)], ext)
		}
	}
}

func calculateDeltas(population [][]int, p *pgs.PGS, targetScore float64) ([]float64, [][]int) {
	var delta float64
	var err error
	deltas := make([]float64, len(population))
	matches := make([][]int, 0)
	for i := range population {
		delta = calculateScore(population[i], p.Weights) - targetScore
		//for delta == 0 {
		for math.Abs(delta) <= pgs.ErrorMargin {
			match := make([]int, len(population[i]))
			copy(match, population[i])
			matches = append(matches, match)
			fmt.Printf("Found match: %s %f\n", arrayTostring(match), calculateScore(match, p.Weights)-targetScore)
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
	for i := 0; i < len(snps); i += pgs.NumHaplotypes {
		for j := 0; j < pgs.NumHaplotypes; j++ {
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
			copy(children[k][:k*pgs.NumHaplotypes], parents[first][:k*pgs.NumHaplotypes])
			copy(children[k][k*pgs.NumHaplotypes:], parents[second][k*pgs.NumHaplotypes:])
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
	for i := 0; i < (len(parents)/2)*2; i += 2 {
		offspring = append(offspring, splice(i, i+1)...)
	}
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

func mutate(population [][]int, deltas []float64, p *pgs.PGS) [][]int {
	//// We shuffle the indices to avoid biasing the mutation towards the first SNPs
	//// which can happen when multiple SNPs have the same weight
	//shuffledIndices := make([]int, p.VariantCount)
	//for i := 0; i < len(shuffledIndices); i++ {
	//	shuffledIndices[i] = i
	//}
	//for i := len(shuffledIndices) - 1; i > 0; i-- {
	//	j := rand.Intn(i + 1)
	//	shuffledIndices[i], shuffledIndices[j] = shuffledIndices[j], shuffledIndices[i]
	//}
	mutated := make([][]int, len(population))
	copy(mutated, population)
	for j, original := range population {
		result := p.MutateGenome(original, deltas[j])
		mutated[j] = result[0]
		if len(result) > 1 {
			mutated = append(mutated, result[1:]...)
		}
	}
	return mutated
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
	for i := 0; i < len(solution); i += pgs.NumHaplotypes {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.NumHaplotypes / float64(len(solution))
}

func arrayTostring(array []int) string {
	str := make([]string, len(array))
	for i := 0; i < len(array); i += pgs.NumHaplotypes {
		str[i] = fmt.Sprint(array[i] + array[i+1])
	}
	return strings.Join(str, "")
}

func sortByAccuracy(solutions map[string][]int, target []int) [][]int {
	i := 0
	flattened := make([][]int, len(solutions))
	accuracies := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		accuracies[i] = accuracy(solution, target)
		i++
	}
	sortBy(flattened, accuracies)
	return flattened
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

func findComplements(solutions map[string][]int, p *pgs.PGS) map[string][]int {
	// Find which positions have the same weight, hence can be swapped
	weightGroups := make(map[float64][]int)
	for i, weight := range p.Weights {
		if _, ok := weightGroups[weight]; !ok {
			weightGroups[weight] = make([]int, 1)
			weightGroups[weight][0] = i
		} else {
			weightGroups[weight] = append(weightGroups[weight], i)
		}
	}
	// Get all duplicated weight positions in a list
	weightCopies := make([][]int, 0)
	for weight, positions := range weightGroups {
		if len(positions) > 1 {
			weightCopies = append(weightCopies, positions)
			fmt.Printf("Weight %f in %d positions: %v\n", weight, len(positions), positions)
		}
	}

	extended := make(map[string][]int)
	for k, solution := range solutions {
		extended[k] = make([]int, len(solution))
		copy(extended[k], solution)
	}

	var leftVal, rightVal int
	for _, solution := range solutions {
		for _, weights := range weightCopies {
			for _, leftPos := range weights {
				leftVal = solution[leftPos*pgs.NumHaplotypes] + solution[leftPos*pgs.NumHaplotypes+1]
				for _, rightPos := range weights {
					rightVal = solution[rightPos*pgs.NumHaplotypes] + solution[rightPos*pgs.NumHaplotypes+1]
					if leftVal != rightVal {
						swapped := make([]int, len(solution))
						copy(swapped, solution)
						swapped[leftPos*pgs.NumHaplotypes], solution[rightPos*pgs.NumHaplotypes] =
							solution[rightPos*pgs.NumHaplotypes], solution[leftPos*pgs.NumHaplotypes]
						swapped[leftPos*pgs.NumHaplotypes+1], solution[rightPos*pgs.NumHaplotypes+1] =
							solution[rightPos*pgs.NumHaplotypes+1], solution[leftPos*pgs.NumHaplotypes+1]
						extend(extended, [][]int{swapped})
					}
				}
			}
		}
	}
	return extended
}

func diploidToSum(diploid []int) []int {
	sum := make([]int, len(diploid)/pgs.NumHaplotypes)
	for i := 0; i < len(diploid); i += pgs.NumHaplotypes {
		sum[i/pgs.NumHaplotypes] = diploid[i] + diploid[i+1]
	}
	return sum
}

//func removeDuplicates(slices [][]int) [][]int {
//	// Create a map to store unique slices
//	uniqueSlices := make(map[string][]int)
//
//	// Iterate over the slices
//	for _, slice := range slices {
//		// Convert the slice to a string representation
//		key := fmt.Sprintf("%v", slice)
//
//		// Check if the slice already exists in the map
//		if _, ok := uniqueSlices[key]; !ok {
//			// If not, add it to the map
//			uniqueSlices[key] = slice
//		}
//	}
//
//	// Create a new slice to store the unique slices
//	uniqueSliceList := make([][]int, 0, len(uniqueSlices))
//
//	// Add the unique slices to the new slice
//	for _, slice := range uniqueSlices {
//		uniqueSliceList = append(uniqueSliceList, slice)
//	}
//
//	return uniqueSliceList
//}
