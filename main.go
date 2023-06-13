package main

import (
	"fmt"
	"log"
	"math"
	"strings"
)

const (
	SCORE      = "score"
	ITERATIONS = 1000
	MARGIN     = 0.0000001
	numThreads = 8
)

var GENOTYPES = [...]int{0, 1, 2}

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
	weights := pgs.Weights

	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = pgs.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	collectChans := make([]chan map[string][]int, numThreads)
	for thread := 0; thread < numThreads; thread++ {
		collectChans[thread] = make(chan map[string][]int)
		go func(t int) {
			solutions := make(map[string][]int, 0)
			// Evaluate candidates
			var delta float64
			//mutated := make([][]int, len(candidates))
			for k := 0; k < ITERATIONS; k++ {
			candidateLoop:
				for j, candidate := range candidates[t*len(candidates)/numThreads : (t+1)*len(candidates)/numThreads] {
					delta = calculateScore(candidate, weights) - targetScore
					if math.Abs(delta) < MARGIN {
						solutions[arrayTostring(candidate)] = candidate
						candidates[j], err = pgs.SampleFromPopulation()
						if err != nil {
							fmt.Println(err)
							return
						}
						continue
					}
					// the weights that cover the delta better, get higher probability of being selected
					// big delta -> bigger weights get higher probability.
					// probability = 1 / abs( delta - weight * snp_old + weight * snp_new )
					possibleMutations := make([]int, 0, len(candidate)*(len(GENOTYPES)-1))
					probs := make([]float64, 0)
					for i, snp := range candidate {
						tmp := delta
						tmp -= weights[i] * float64(snp)
						for _, v := range GENOTYPES {
							if v == snp {
								continue
							}
							possibleMutations = append(possibleMutations, v)
							if math.Abs(tmp+weights[i]*float64(v)) < MARGIN {
								candidate[i] = v
								solutions[arrayTostring(candidate)] = candidate
								candidates[j], err = pgs.SampleFromPopulation()
								if err != nil {
									fmt.Println(err)
									return
								}
								continue candidateLoop
							}
							probs = append(probs, 1/math.Abs(tmp+weights[i]*float64(v)))
						}
					}
					//mutated[j] = pgs.MutateGenome(candidate, possibleMutations, probs)
					candidates[j] = pgs.MutateGenome(candidate, possibleMutations, probs)
				}
			}
			collectChans[t] <- solutions
		}(thread)
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
