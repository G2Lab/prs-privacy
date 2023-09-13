package main

import (
	"context"
	"fmt"
	"log"
	"math"
	"math/rand"
	"strings"
	"time"

	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

const ITERATIONS = 5000
const numCpusEnv = "SLURM_CPUS_PER_TASK"

type Solver struct {
	target          float64
	p               *pgs.PGS
	populationChans []chan [][]uint8
	deltaChans      []chan []float64
}

func NewSolver(ctx context.Context, target float64, p *pgs.PGS, numThreads int) *Solver {
	s := &Solver{
		target: target,
		p:      p,
	}
	//s.populationChans = make([]chan [][]int, numThreads)
	//s.deltaChans = make([]chan []float64, numThreads)
	//for i := 0; i < numThreads; i++ {
	//	s.populationChans[i] = make(chan [][]int)
	//	s.deltaChans[i] = make(chan []float64)
	//	go s.mutate(ctx, s.populationChans[i], s.deltaChans[i])
	//}
	return s
}

func main() {
	//INDIVIDUAL := "NA18595"
	INDIVIDUAL := "HG02182" // lowest score for PGS000040
	//INDIVIDUAL := "HG02215" // highest score for PGS000040
	//INDIVIDUAL := "HG02728" // middle 648
	//INDIVIDUAL := "NA19780" // high 648
	//INDIVIDUAL := "HG00551" // low 648

	p := pgs.NewPGS()
	//catalogFile := "PGS000073_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000037_hmPOS_GRCh38.txt"
	catalogFile := "PGS000040_hmPOS_GRCh38.txt"
	//catalogFile := "PGS000648_hmPOS_GRCh38.txt"
	//catalogFile := "PGS002302_hmPOS_GRCh38.txt"
	err := p.LoadCatalogFile(catalogFile)
	if err != nil {
		log.Printf("Error loading catalog file: %v\n", err)
		return
	}
	fmt.Printf("%s, %s\n", p.PgsID, INDIVIDUAL)
	p.LoadStats()
	cohort := NewCohort(p)

	//numThreadsEnv := os.Getenv(numCpusEnv)
	//var numThreads int
	//if numThreadsEnv != "" {
	//	numThreads, err = strconv.Atoi(numThreadsEnv)
	//	if err != nil {
	//		log.Printf("Error parsing numCpus %s: %v\n", numCpusEnv, err)
	//		return
	//	}
	//} else {
	//	log.Printf("Could not read numCpus env %s: %v\n", numCpusEnv, err)
	//	return
	//}

	numThreads := 1

	ctx, cancel := context.WithCancel(context.Background())
	solver := NewSolver(ctx, cohort[INDIVIDUAL].Score, p, numThreads)

	solmap := solver.dp(numThreads)
	//solmap := solver.solve(numThreads)
	//solmap := solver.recursive(numThreads)
	cancel()
	//solmap = findComplements(solmap, p, numThreads)
	solutions := sortByAccuracy(solmap, cohort[INDIVIDUAL].Genotype)
	fmt.Printf("\nTrue:\n%s -- %f, %f\n", arrayToString(cohort[INDIVIDUAL].Genotype),
		cohort[INDIVIDUAL].Score-calculateScore(cohort[INDIVIDUAL].Genotype, p.Weights),
		p.CalculateSequenceLikelihood(cohort[INDIVIDUAL].Genotype))
	fmt.Printf("Guessed %d:\n", len(solutions))
	for _, solution := range solutions {
		fmt.Printf("%s -- %.3f, %.8f, %.2f\n", arrayToString(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
			cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights), p.CalculateSequenceLikelihood(solution))
		//fmt.Printf("%s -- %.3f, %.8f\n", arrayToStringDiploid(solution), accuracy(solution, cohort[INDIVIDUAL].Genotype),
		//	cohort[INDIVIDUAL].Score-calculateScore(solution, p.Weights))
	}
}

func (s *Solver) dp(numThreads int) map[string][]uint8 {
	multiplier := math.Pow(10, float64(s.p.WeightPrecision))

	errorMargin := int64(math.Pow(10, float64(s.p.WeightPrecision/4)))
	target := int64(s.target * multiplier)
	weights := make([]int64, len(s.p.Weights))
	for i, w := range s.p.Weights {
		weights[i] = int64(w * multiplier)
	}
	minWeight := findMin(weights)
	negativeTotal := sumUpNegatives(weights)

	fmt.Printf("Target: %d±%d\n", target, errorMargin)
	fmt.Printf("Weights: %v\n", weights)

	// Fill the dp table using dynamic programming
	table := make(map[int64][]uint16)
	// add the zero weight
	table[0] = make([]uint16, 0)
	var lowerBound int64
	if negativeTotal == 0 {
		lowerBound = target - errorMargin - minWeight
	} else {
		lowerBound = target + errorMargin - negativeTotal
	}
	targetLeft := target - errorMargin
	targetRight := target + errorMargin
	var prevSum, nextSum int64
	for i := 0; i < len(weights); i++ {
		fmt.Printf("Position %d/%d\n", i+1, len(weights))
		existingSums := make([]int64, 0, len(table))
		for w := range table {
			existingSums = append(existingSums, w)
		}
		for _, prevSum = range existingSums {
			for k := 1; k <= pgs.NumHaplotypes; k++ {
				nextSum = prevSum + int64(k)*weights[i]
				if nextSum <= lowerBound || (nextSum >= targetLeft && nextSum <= targetRight) {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
					}
					table[nextSum] = append(table[nextSum], uint16(pgs.NumHaplotypes*i+k-1))
				}
			}
		}
	}

	solutions := make(map[string][]uint8)
	addToSolutions := func(path []uint16) {
		sol := make([]uint8, len(s.p.Weights)*pgs.NumHaplotypes)
		for _, pos := range path {
			sol[pos] = 1
			if pos%pgs.NumHaplotypes == 1 {
				sol[pos-1] = 1
			}
		}
		solutions[arrayToString(sol)] = sol
	}
	fmt.Printf("Number of weights in the table %d\n", len(table))
	//fmt.Printf("%d: %v\n", target, table[target])

	//backtracking
	//var prevLen int
	var weight int64
	var backtrack func(int64, []uint16)
	backtrack = func(sum int64, path []uint16) {
		pointers, ok := table[sum]
		if int64(math.Abs(float64(sum))) <= errorMargin {
			addToSolutions(path)
			//if len(solutions)%10 == 0 && len(solutions) != prevLen {
			//	fmt.Printf("Found %d solutions\n", len(solutions))
			//}
			//prevLen = len(solutions)
			return
		}
		if !ok {
			pointers = make([]uint16, 0)
			// Due to the loss of precision, sometimes the sum in the table is slightly different
			for w := sum - errorMargin; w <= sum+errorMargin; w++ {
				if pos, exists := table[w]; exists {
					for _, p := range pos {
						pointers = append(pointers, p)
					}
				}
			}
		}
		for _, p := range pointers {
			// Make sure that we do not have paths with two values for the same snp
			if isNotPermitted(p, path) {
				continue
			}
			newPath := make([]uint16, len(path)+1)
			copy(newPath, path)
			newPath[len(path)] = p
			weight = weights[p/2] * (1 + int64(p%2))
			backtrack(sum-weight, newPath)
		}
	}
	backtrack(target, make([]uint16, 0))

	//var printer func(float64)
	//printer = func(w float64) {
	//	fmt.Printf("%g: ", w)
	//	for _, pos := range table[w] {
	//		fmt.Printf("%d ", pos)
	//	}
	//	fmt.Println()
	//	for _, pos := range table[w] {
	//		fmt.Printf("%d: ", pos)
	//		weight = s.p.Weights[pos/2] * (1 + float64(pos%2))
	//		printer(w - weight)
	//	}
	//}

	return solutions
}

func isNotPermitted(v uint16, array []uint16) bool {
	for _, a := range array {
		if a == v || (v%pgs.NumHaplotypes == 0 && a == v+1) || (v%pgs.NumHaplotypes == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func findMin(values []int64) int64 {
	minV := values[0]
	for _, v := range values {
		if v < minV {
			minV = v
		}
	}
	return minV
}

func sumUpNegatives(values []int64) int64 {
	sum := int64(0)
	for _, v := range values {
		if v < 0 {
			sum += v
		}
	}
	return sum
}

//func (s *Solver) recursive(numThreads int) map[string][]int {
//	var wg sync.WaitGroup
//	solutions := NewMutexMap(make(map[string][]int))
//	initial := make([]int, 0, len(s.p.Weights))
//	s.branching(initial,0.0, s.target, 0, solutions, numThreads, &wg)
//	wg.Wait()
//	return solutions.RetrieveMapOnly()
//}

//func (s *Solver) branching(current []int, sum, target float64, pos int, solutions *MutexMap, threadsLeft int, wg *sync.WaitGroup) {
//	if math.Abs(sum - target) < pgs.ErrorMargin {
//		if len(current) < len(s.p.Weights) {
//			// pad with zeros
//			current = append(current, make([]int, len(s.p.Weights)-len(current))...)
//		}
//		//fmt.Printf("Found solution: %s\n", arrayToStringDiploid(current))
//		solutions.Put(arrayToStringDiploid(current), current)
//		return
//	}
//	if sum > target || pos >= len(s.p.Weights) {
//		return
//	}
//
//	current = append(current, 0)
//	if threadsLeft > 0 {
//		wg.Add(1)
//		threadsLeft -= 1
//		go func() {
//			s.branching(current, sum, target, pos+1, solutions, threadsLeft, wg)
//			wg.Done()
//		}()
//	} else {
//		s.branching(current, sum, target, pos+1, solutions, threadsLeft, wg)
//	}
//	for snp := 1; snp <= 2; snp++ {
//		branched := make([]int, len(current))
//		copy(branched, current)
//		branched[len(branched)-1] = snp
//		if threadsLeft > 0 {
//			wg.Add(1)
//			threadsLeft -= 1
//			go func() {
//				s.branching(branched, sum+float64(snp)*s.p.Weights[pos], target, pos+1, solutions, threadsLeft, wg)
//				wg.Done()
//			}()
//		} else {
//			s.branching(branched, sum+float64(snp)*s.p.Weights[pos], target, pos+1, solutions, threadsLeft, wg)
//		}
//	}
//}

func (s *Solver) solve(numThreads int) map[string][]uint8 {
	var err error
	rand.NewSource(time.Now().UnixNano())
	poolSize := len(s.p.Variants) * 100
	solutions := make(map[string][]uint8)
	candidates := make([][]uint8, poolSize, 2*poolSize)
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
		mutated := make([][]uint8, 0, len(candidates))
		for thread := 0; thread < numThreads; thread++ {
			mutated = append(mutated, <-s.populationChans[thread]...)
		}
		candidates = make([][]uint8, len(mutated))
		copy(candidates, mutated)
	}

	return solutions
}

func extend(base map[string][]uint8, extension [][]uint8) {
	for _, ext := range extension {
		if _, ok := base[arrayToString(ext)]; !ok {
			base[arrayToString(ext)] = make([]uint8, len(ext))
			copy(base[arrayToString(ext)], ext)
		}
	}
}

func (s *Solver) calculateDeltas(population [][]uint8) ([]float64, [][]uint8) {
	var delta float64
	var err error
	deltas := make([]float64, len(population))
	matches := make([][]uint8, 0)
	for i := range population {
		delta = calculateScore(population[i], s.p.Weights) - s.target
		//for delta == 0 {
		for math.Abs(delta) <= pgs.ErrorMargin {
			match := make([]uint8, len(population[i]))
			copy(match, population[i])
			matches = append(matches, match)
			//fmt.Printf("Found match: %s %f\n", arrayToString(match), calculateScore(match, s.p.Weights)-s.target)
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

//func calculateScore(snps []int, weights []float64) float64 {
//	score := 0.0
//	for i := 0; i < len(snps); i++ {
//			switch snps[i] {
//			case 0:
//				continue
//			case 1:
//				score += weights[i]
//			case 2:
//				score += weights[i] + weights[i]
//			default:
//				log.Printf("Invalid alelle value: %d", snps[i])
//			}
//	}
//	return score
//}

func calculateScore(snps []uint8, weights []float64) float64 {
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

func (s *Solver) crossover(population [][]uint8) [][]uint8 {
	parents := tools.Shuffle(population)
	splice := func(first, second int) [][]uint8 {
		children := make([][]uint8, len(parents[first])/2)
		for k := 0; k < len(children); k++ {
			children[k] = make([]uint8, len(parents[first]))
			copy(children[k][:k*pgs.NumHaplotypes], parents[first][:k*pgs.NumHaplotypes])
			copy(children[k][k*pgs.NumHaplotypes:], parents[second][k*pgs.NumHaplotypes:])
		}
		deltas, matches := s.calculateDeltas(children)
		if len(matches) > 0 {
			return matches
		}
		fitness := deltasToFitness(deltas)
		pickedIndex := tools.SampleFromDistribution(fitness)
		return [][]uint8{children[pickedIndex]}
	}
	offspring := make([][]uint8, 0, len(parents))
	for i := 0; i < (len(parents)/2)*2; i += 2 {
		offspring = append(offspring, splice(i, i+1)...)
	}
	return offspring
}

func tournament(population [][]uint8, deltas []float64, populationSize int) [][]uint8 {
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

func (s *Solver) mutate(ctx context.Context, popChan chan [][]uint8, deltaChan chan []float64) {
	timeout := 1 * time.Second
	t := time.NewTimer(timeout)
	for {
		var population [][]uint8
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
		mutated := make([][]uint8, 0, len(population))
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

//func accuracy(solution []int, target []int) float64 {
//	if len(solution) != len(target) {
//		return 0.0
//	}
//	acc := 0.0
//	for i := 0; i < len(solution); i++ {
//		if solution[i] == target[i] {
//			acc++
//		}
//	}
//	return acc / float64(len(solution))
//}

func accuracy(solution []uint8, target []uint8) float64 {
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

func findComplements(solutions map[string][]uint8, p *pgs.PGS, numThreads int) map[string][]uint8 {
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
		//fmt.Printf("Exploring %s\n", arrayToString(solution))
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
	//			//fmt.Printf("Exploring %s\n", arrayToString(solution))
	//			explore(solution, weightCopies, mutexed)
	//		}
	//	}()
	//}
	//wg.Wait()
	return mutexed.RetrieveMapOnly()
}

func explore(source []uint8, positions [][]int, saver *MutexMap) {
	var leftPos, rightPos int
	var leftVal, rightVal uint8
	//fmt.Printf("Exploring %s\n", arrayToString(source))
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
				swapped := make([]uint8, len(source))
				copy(swapped, source)
				swapped[leftPos*pgs.NumHaplotypes], swapped[rightPos*pgs.NumHaplotypes] =
					swapped[rightPos*pgs.NumHaplotypes], swapped[leftPos*pgs.NumHaplotypes]
				swapped[leftPos*pgs.NumHaplotypes+1], swapped[rightPos*pgs.NumHaplotypes+1] =
					swapped[rightPos*pgs.NumHaplotypes+1], swapped[leftPos*pgs.NumHaplotypes+1]
				saver.Put(arrayToString(swapped), swapped)
				if len(positions) > 1 {
					explore(swapped, positions[1:], saver)
				}
				// second, splitting 2 + 0 into 1 + 1
				if (leftVal == 2 && rightVal == 0) || (leftVal == 0 && rightVal == 2) {
					split := make([]uint8, len(source))
					copy(split, source)
					split[leftPos*pgs.NumHaplotypes], split[leftPos*pgs.NumHaplotypes+1] = 1, 0
					split[rightPos*pgs.NumHaplotypes], split[rightPos*pgs.NumHaplotypes+1] = 0, 1
					saver.Put(arrayToString(split), split)
					if len(positions) > 1 {
						explore(split, positions[1:], saver)
					}
				}
			//	converting 1 + 1 into 2 + 0 and 0 + 2
			case leftVal == 1 && rightVal == 1:
				leftPushed, rightPushed := make([]uint8, len(source)), make([]uint8, len(source))
				copy(leftPushed, source)
				copy(rightPushed, source)
				leftPushed[leftPos*pgs.NumHaplotypes], leftPushed[leftPos*pgs.NumHaplotypes+1] = 1, 1
				leftPushed[rightPos*pgs.NumHaplotypes], leftPushed[rightPos*pgs.NumHaplotypes+1] = 0, 0
				saver.Put(arrayToString(leftPushed), leftPushed)
				rightPushed[leftPos*pgs.NumHaplotypes], rightPushed[leftPos*pgs.NumHaplotypes+1] = 0, 0
				rightPushed[rightPos*pgs.NumHaplotypes], rightPushed[rightPos*pgs.NumHaplotypes+1] = 1, 1
				saver.Put(arrayToString(rightPushed), rightPushed)
				if len(positions) > 1 {
					explore(leftPushed, positions[1:], saver)
					explore(rightPushed, positions[1:], saver)
				}
			default:
				continue
			}
		}
	}
	//if len(positions[0]) < 3 {
	//	return
	//}
	//triplets := getTriplets(positions[0])
	//for _, triplet := range triplets {
	//	values := []int{source[triplet[0]*pgs.NumHaplotypes] + source[triplet[0]*pgs.NumHaplotypes+1],
	//		source[triplet[1]*pgs.NumHaplotypes] + source[triplet[1]*pgs.NumHaplotypes+1],
	//		source[triplet[2]*pgs.NumHaplotypes] + source[triplet[2]*pgs.NumHaplotypes+1]}
	//	sortInts(triplet, values)
	//	if values[0] == 0 && values[1] == 0 && values[2] == 2 {
	//		distributed := make([]uint8, len(source))
	//		copy(distributed, source)
	//		distributed[triplet[0]*pgs.NumHaplotypes], distributed[triplet[0]*pgs.NumHaplotypes+1] = 0, 1
	//		distributed[triplet[1]*pgs.NumHaplotypes], distributed[triplet[1]*pgs.NumHaplotypes+1] = 0, 1
	//		distributed[triplet[2]*pgs.NumHaplotypes], distributed[triplet[2]*pgs.NumHaplotypes+1] = 0, 0
	//		saver.Put(arrayToString(distributed), distributed)
	//		if len(positions) > 1 {
	//			explore(distributed, positions[1:], saver)
	//		}
	//	}
	//}
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

func arrayToString(array []uint8) string {
	str := make([]string, len(array)/2)
	for i := 0; i < len(array); i += pgs.NumHaplotypes {
		str[i/2] = fmt.Sprint(array[i] + array[i+1])
	}
	return strings.Join(str, "")
}

func arrayToStringDiploid(array []int) string {
	str := make([]string, len(array))
	for i := 0; i < len(array); i++ {
		str[i] = fmt.Sprint(array[i])
	}
	return strings.Join(str, "")
}

func sortByAccuracy(solutions map[string][]uint8, target []uint8) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	accuracies := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		accuracies[i] = accuracy(solution, target)
		i++
	}
	sortBy(flattened, accuracies)
	return flattened
}

func sortBy(items [][]uint8, properties []float64) {
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
