package main

import (
	"context"
	"fmt"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
)

const ITERATIONS = 5000
const numCpusEnv = "SLURM_CPUS_PER_TASK"

type Solver struct {
	target          float64
	p               *pgs.PGS
	populationChans []chan [][]int
	deltaChans      []chan []float64
}

func NewSolver(ctx context.Context, target float64, p *pgs.PGS, numThreads int) *Solver {
	s := &Solver{
		target: target,
		p:      p,
	}
	s.populationChans = make([]chan [][]int, numThreads)
	s.deltaChans = make([]chan []float64, numThreads)
	for i := 0; i < numThreads; i++ {
		s.populationChans[i] = make(chan [][]int)
		s.deltaChans[i] = make(chan []float64)
		go s.mutate(ctx, s.populationChans[i], s.deltaChans[i])
	}
	return s
}

func main() {
	//INDIVIDUAL := "NA18595"
	INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040

	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	p.LoadStats()
	cohort := NewCohort(p)

	numThreadsEnv := os.Getenv(numCpusEnv)
	var numThreads int
	if numThreadsEnv != "" {
		numThreads, err = strconv.Atoi(numThreadsEnv)
		if err != nil {
			log.Printf("Error parsing numCpus %s: %v\n", numCpusEnv, err)
			return
		}
	} else {
		log.Printf("Could not read numCpus env %s: %v\n", numCpusEnv, err)
		return
	}

	ctx, cancel := context.WithCancel(context.Background())
	solver := NewSolver(ctx, cohort[INDIVIDUAL].Score, p, numThreads)

	solmap := solver.solve(numThreads)
	cancel()
	solmap = findComplements(solmap, p, numThreads)
	solutions := sortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f, %f\n", arrayTostring(cohort[INDIVIDUAL].Genotype),
		cohort[INDIVIDUAL].Score-calculateScore(cohort[INDIVIDUAL].Genotype, p.Weights),
		p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Guessed %d:\n", len(solutions))
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.6f, %.2f\n", arrayTostring(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
	}
	cancel()
}

func (s *Solver) solve(numThreads int) map[string][]int {
	var err error
	rand.NewSource(time.Now().UnixNano())
	poolSize := len(s.p.Variants) * 100
	solutions := make(map[string][]int)
	candidates := make([][]int, poolSize, 2*poolSize)
	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = s.p.SampleFromPopulation()
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	//Evaluate candidates
	for k := 0; k < ITERATIONS; k++ {
		if k%500 == 0 {
			fmt.Printf("Iteration %d/%d\n", k, ITERATIONS)
		}
		deltas, solved := s.calculateDeltas(candidates)
		extend(solutions, solved)
		children := s.crossover(candidates)
		chDeltas, solved := s.calculateDeltas(children)
		extend(solutions, solved)
		candidates = tournament(append(candidates, children...), append(deltas, chDeltas...), poolSize)
		deltas, solved = s.calculateDeltas(candidates)
		extend(solutions, solved)
		for thread := 0; thread < numThreads; thread++ {
			s.populationChans[thread] <- candidates[thread*len(candidates)/numThreads : (thread+1)*len(candidates)/numThreads]
			s.deltaChans[thread] <- deltas[thread*len(candidates)/numThreads : (thread+1)*len(candidates)/numThreads]
		}
		mutated := make([][]int, 0, len(candidates))
		for thread := 0; thread < numThreads; thread++ {
			mutated = append(mutated, <-s.populationChans[thread]...)
		}
		candidates = make([][]int, len(mutated))
		copy(candidates, mutated)
	}

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

func (s *Solver) calculateDeltas(population [][]int) ([]float64, [][]int) {
	var delta float64
	var err error
	deltas := make([]float64, len(population))
	matches := make([][]int, 0)
	for i := range population {
		delta = calculateScore(population[i], s.p.Weights) - s.target
		//for delta == 0 {
		for math.Abs(delta) <= pgs.ErrorMargin {
			match := make([]int, len(population[i]))
			copy(match, population[i])
			matches = append(matches, match)
			//fmt.Printf("Found match: %s %f\n", arrayTostring(match), calculateScore(match, s.p.Weights)-s.target)
			population[i], err = s.p.SampleFromPopulation()
			if err != nil {
				log.Printf("Error resampling in fitness calculation: %v\n", err)
				return nil, nil
			}
			delta = calculateScore(population[i], s.p.Weights) - s.target
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

func (s *Solver) crossover(population [][]int) [][]int {
	parents := tools.Shuffle(population)
	splice := func(first, second int) [][]int {
		children := make([][]int, len(parents[first])/2)
		for k := 0; k < len(children); k++ {
			children[k] = make([]int, len(parents[first]))
			copy(children[k][:k*pgs.NumHaplotypes], parents[first][:k*pgs.NumHaplotypes])
			copy(children[k][k*pgs.NumHaplotypes:], parents[second][k*pgs.NumHaplotypes:])
		}
		deltas, matches := s.calculateDeltas(children)
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

func (s *Solver) mutate(ctx context.Context, popChan chan [][]int, deltaChan chan []float64) {
	timeout := 1 * time.Second
	t := time.NewTimer(timeout)
	for {
		var population [][]int
		select {
		case <-ctx.Done():
			return
		case <-t.C:
			t.Reset(timeout)
			continue
		case population = <-popChan:
			t.Reset(timeout)
		}
		deltas := <-deltaChan
		mutated := make([][]int, 0, len(population))
		for j, original := range population {
			result := s.p.MutateGenome(original, deltas[j])
			mutated = append(mutated, result...)
		}
		popChan <- mutated
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
	for i := 0; i < len(solution); i += pgs.NumHaplotypes {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.NumHaplotypes / float64(len(solution))
}

func findComplements(solutions map[string][]int, p *pgs.PGS, numThreads int) map[string][]int {
	// Find which positions have the same weight, hence can be swapped
	weightGroups := make(map[float64][]int)
	for i, weight := range p.Weights {
		if _, ok := weightGroups[weight]; !ok {
			weightGroups[weight] = make([]int, 1, 2)
			weightGroups[weight][0] = i
		} else {
			weightGroups[weight] = append(weightGroups[weight], i)
		}
	}
	// Get all identical weight positions in a list
	weightCopies := make([][]int, 0)
	for weight, positions := range weightGroups {
		if len(positions) > 1 {
			weightCopies = append(weightCopies, positions)
			fmt.Printf("Weight %f in %d positions: %v\n", weight, len(positions), positions)
		}
	}

	mutexed := NewMutexMap(solutions)
	for _, solution := range solutions {
		//fmt.Printf("Exploring %s\n", arrayTostring(solution))
		explore(solution, weightCopies, mutexed)
	}
	//solSlice := make([][]int, 0, len(solutions))
	//for _, solution := range solutions {
	//	solSlice = append(solSlice, solution)
	//}
	//var wg sync.WaitGroup
	//for thread := 0; thread < numThreads; thread++ {
	//	go func() {
	//		wg.Add(1)
	//		defer wg.Done()
	//		for _, solution := range solSlice[thread*len(solutions)/numThreads : (thread+1)*len(solutions)/numThreads] {
	//			//fmt.Printf("Exploring %s\n", arrayTostring(solution))
	//			explore(solution, weightCopies, mutexed)
	//		}
	//	}()
	//}
	//wg.Wait()
	return mutexed.RetrieveMapOnly()
}

func explore(source []int, positions [][]int, saver *MutexMap) {
	var leftPos, rightPos, leftVal, rightVal int
	//fmt.Printf("Exploring %s\n", arrayTostring(source))
	//fmt.Printf("Num positions: %d\n", len(positions[0]))
	if len(positions) > 1 {
		// Sending further one unchanged version
		explore(source, positions[1:], saver)
	}
	// pairwise search
	for i := 0; i < len(positions[0])-1; i++ {
		leftPos = positions[0][i]
		leftVal = source[leftPos*pgs.NumHaplotypes] + source[leftPos*pgs.NumHaplotypes+1]
		for j := i + 1; j < len(positions[0]); j++ {
			rightPos = positions[0][j]
			rightVal = source[rightPos*pgs.NumHaplotypes] + source[rightPos*pgs.NumHaplotypes+1]
			switch {
			case leftVal != rightVal:
				// first, simple swapping
				swapped := make([]int, len(source))
				copy(swapped, source)
				swapped[leftPos*pgs.NumHaplotypes], swapped[rightPos*pgs.NumHaplotypes] =
					swapped[rightPos*pgs.NumHaplotypes], swapped[leftPos*pgs.NumHaplotypes]
				swapped[leftPos*pgs.NumHaplotypes+1], swapped[rightPos*pgs.NumHaplotypes+1] =
					swapped[rightPos*pgs.NumHaplotypes+1], swapped[leftPos*pgs.NumHaplotypes+1]
				saver.Put(arrayTostring(swapped), swapped)
				if len(positions) > 1 {
					explore(swapped, positions[1:], saver)
				}
				// second, splitting 2 + 0 into 1 + 1
				if (leftVal == 2 && rightVal == 0) || (leftVal == 0 && rightVal == 2) {
					split := make([]int, len(source))
					copy(split, source)
					split[leftPos*pgs.NumHaplotypes], split[leftPos*pgs.NumHaplotypes+1] = 1, 0
					split[rightPos*pgs.NumHaplotypes], split[rightPos*pgs.NumHaplotypes+1] = 0, 1
					saver.Put(arrayTostring(split), split)
					if len(positions) > 1 {
						explore(split, positions[1:], saver)
					}
				}
			//	converting 1 + 1 into 2 + 0 and 0 + 2
			case leftVal == 1 && rightVal == 1:
				leftPushed, rightPushed := make([]int, len(source)), make([]int, len(source))
				copy(leftPushed, source)
				copy(rightPushed, source)
				leftPushed[leftPos*pgs.NumHaplotypes], leftPushed[leftPos*pgs.NumHaplotypes+1] = 1, 1
				leftPushed[rightPos*pgs.NumHaplotypes], leftPushed[rightPos*pgs.NumHaplotypes+1] = 0, 0
				saver.Put(arrayTostring(leftPushed), leftPushed)
				rightPushed[leftPos*pgs.NumHaplotypes], rightPushed[leftPos*pgs.NumHaplotypes+1] = 0, 0
				rightPushed[rightPos*pgs.NumHaplotypes], rightPushed[rightPos*pgs.NumHaplotypes+1] = 1, 1
				saver.Put(arrayTostring(rightPushed), rightPushed)
				if len(positions) > 1 {
					explore(leftPushed, positions[1:], saver)
					explore(rightPushed, positions[1:], saver)
				}
			default:
				continue
			}
		}
	}
	if len(positions[0]) < 3 {
		return
	}
	triplets := getTriplets(positions[0])
	for _, triplet := range triplets {
		values := []int{source[triplet[0]*pgs.NumHaplotypes] + source[triplet[0]*pgs.NumHaplotypes+1],
			source[triplet[1]*pgs.NumHaplotypes] + source[triplet[1]*pgs.NumHaplotypes+1],
			source[triplet[2]*pgs.NumHaplotypes] + source[triplet[2]*pgs.NumHaplotypes+1]}
		sortInts(triplet, values)
		if values[0] == 0 && values[1] == 0 && values[2] == 2 {
			distributed := make([]int, len(source))
			copy(distributed, source)
			distributed[triplet[0]*pgs.NumHaplotypes], distributed[triplet[0]*pgs.NumHaplotypes+1] = 0, 1
			distributed[triplet[1]*pgs.NumHaplotypes], distributed[triplet[1]*pgs.NumHaplotypes+1] = 0, 1
			distributed[triplet[2]*pgs.NumHaplotypes], distributed[triplet[2]*pgs.NumHaplotypes+1] = 0, 0
			saver.Put(arrayTostring(distributed), distributed)
			if len(positions) > 1 {
				explore(distributed, positions[1:], saver)
			}
		}
	}
}

func sortInts(positions []int, values []int) {
	for i := 0; i < len(values)-1; i++ {
		for j := i + 1; j < len(values); j++ {
			if values[i] > values[j] {
				values[i], values[j] = values[j], values[i]
				positions[i], positions[j] = positions[j], positions[i]
			}
		}
	}
}

func getTriplets(nums []int) [][]int {
	triplets := make([][]int, 0)

	for i := 0; i < len(nums)-2; i++ {
		for j := i + 1; j < len(nums)-1; j++ {
			for k := j + 1; k < len(nums); k++ {
				triplet := []int{nums[i], nums[j], nums[k]}
				triplets = append(triplets, triplet)
			}
		}
	}

	return triplets
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

func diploidToSum(diploid []int) []int {
	sum := make([]int, len(diploid)/pgs.NumHaplotypes)
	for i := 0; i < len(diploid); i += pgs.NumHaplotypes {
		sum[i/pgs.NumHaplotypes] = diploid[i] + diploid[i+1]
	}
	return sum
}
