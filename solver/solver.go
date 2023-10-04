package solver

import (
	"context"
	"fmt"
	"github.com/nikirill/prs/params"
	"log"
	"math"
	"math/rand"
	"strings"
	"sync"
	"time"

	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

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

func (s *Solver) DP(numThreads int) map[string][]uint8 {
	var roundedMode = false
	var multiplier float64
	var roundingError int64 = 0
	var errorMargin int64 = 0
	var preciseTarget int64
	var preciseWeights []int64
	var preciseMultiplier float64
	if s.p.WeightPrecision > params.PrecisionsLimit {
		multiplier = math.Pow(10, params.PrecisionsLimit)
		preciseMultiplier = math.Pow(10, float64(s.p.WeightPrecision))
		preciseTarget = int64(s.target * preciseMultiplier)
		preciseWeights = make([]int64, len(s.p.Weights))
		for i, w := range s.p.Weights {
			preciseWeights[i] = int64(w * preciseMultiplier)
		}
		roundedMode = true
		roundingError = int64(s.p.VariantCount) * 4 / 5
		if s.p.WeightPrecision > params.FloatPrecision {
			errorMargin = 1e2
		}
	} else {
		multiplier = math.Pow(10, float64(s.p.WeightPrecision))
	}

	target := int64(s.target * multiplier)
	weights := make([]int64, len(s.p.Weights))
	for i, w := range s.p.Weights {
		weights[i] = int64(w * multiplier)
	}
	targetLeft := target - roundingError
	targetRight := target + roundingError
	maxTotalPositive, maxTotalNegative := getMaxTotal(weights)
	upperBound := targetRight - maxTotalNegative
	lowerBound := targetLeft - maxTotalPositive
	//remainingMinValues := findNextPositiveMins(weights)

	fmt.Printf("Target: %d±%d\n", target, roundingError)
	fmt.Printf("Weights: %v\n", weights)
	//fmt.Printf("Upper bound: %d\n", maxTotalPositive)
	//fmt.Printf("Lower bound: %d\n", maxTotalNegative)
	//fmt.Printf("Next mins: %v\n", remainingMinValues)

	// Fill the DP table using dynamic programming
	table := make(map[int64][]uint16)
	// add the zero weight
	table[0] = make([]uint16, 0)
	var prevSum, nextSum int64
	for i := 0; i < len(weights); i++ {
		fmt.Printf("Position %d/%d\n", i+1, len(weights))
		if weights[i] > 0 {
			lowerBound += pgs.NumHaplotypes * weights[i]
		} else {
			upperBound += pgs.NumHaplotypes * weights[i]
		}
		existingSums := make([]int64, 0, len(table))
		for w := range table {
			existingSums = append(existingSums, w)
		}
		for _, prevSum = range existingSums {
			for k := 1; k <= pgs.NumHaplotypes; k++ {
				nextSum = prevSum + int64(k)*weights[i]
				//if (nextSum <= upperBound-remainingMinValues[i] && nextSum >= lowerBound) || (nextSum >= targetLeft && nextSum <= targetRight) {
				if (nextSum >= lowerBound && nextSum <= upperBound) || (nextSum >= targetLeft && nextSum <= targetRight) {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
					}
					table[nextSum] = append(table[nextSum], uint16(pgs.NumHaplotypes*i+k-1))
				}
			}
		}
	}

	var solMutex, threadMutex sync.Mutex
	//sf, err := os.OpenFile(fmt.Sprintf("solutions_%s.txt", s.p.PgsID), os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	//if err != nil {
	//	log.Printf("Error opening solutions file: %v\n", err)
	//	return nil
	//}
	//defer sf.Close()
	solutions := make(map[string][]uint8)
	addToSolutions := func(path []uint16) {
		sol := make([]uint8, len(s.p.Weights)*pgs.NumHaplotypes)
		for _, pos := range path {
			sol[pos] = 1
			if pos%pgs.NumHaplotypes == 1 {
				sol[pos-1] = 1
			}
		}
		// If we are in the roundedMode mode, check that the precise weights add up to the target
		if roundedMode && int64(math.Abs(float64(getScore(sol, preciseWeights)-preciseTarget))) > errorMargin {
			return
		}
		solMutex.Lock()
		//_, err = sf.WriteString(fmt.Sprintf("%s\n", ArrayToString(sol)))
		//if err != nil {
		//	log.Printf("Error writing %s to solutions file: %v\n", ArrayToString(sol), err)
		//}
		solutions[ArrayToString(sol)] = sol
		solMutex.Unlock()
	}
	fmt.Printf("Number of weights in the table %d\n", len(table))

	//backtracking
	var wg sync.WaitGroup
	var weight int64
	var backtrack func(int64, []uint16, bool, int)
	backtrack = func(sum int64, path []uint16, firstLevel bool, threadsAvailable int) {
		if int64(math.Abs(float64(sum))) <= roundingError {
			addToSolutions(path)
			return
		}
		stage := make(map[int64][]uint16)
		pointers, ok := table[sum]
		if ok {
			stage[sum] = pointers
		}
		// Due to the loss of precision, sometimes the sum in the table is slightly different
		if !ok || (firstLevel && roundedMode) {
			for w := sum - roundingError; w < sum+roundingError; w++ {
				if w == sum {
					continue
				}
				if pos, exists := table[w]; exists {
					stage[w] = pos
				}
			}
		}
		for w, pts := range stage {
			for j, p := range pts {
				// Make sure that we do not have paths with two values for the same snp
				if locusAlreadyInSlice(p, path) || (len(path) > 0 && p > path[len(path)-1]) {
					continue
				}
				if firstLevel {
					fmt.Printf("%d: %d/%d ", w, j+1, len(pts))
				}
				newPath := make([]uint16, len(path)+1)
				copy(newPath, path)
				newPath[len(path)] = p
				weight = weights[p/2] * (1 + int64(p%2))
				if firstLevel && threadsAvailable > 1 {
					threadMutex.Lock()
					threadsAvailable -= 1
					threadMutex.Unlock()
					wg.Add(1)
					go func(t int64) {
						defer wg.Done()
						backtrack(t, newPath, false, threadsAvailable)
						threadMutex.Lock()
						threadsAvailable += 1
						threadMutex.Unlock()
					}(w - weight)
				} else {
					backtrack(w-weight, newPath, false, threadsAvailable)
				}
			}
		}
	}
	backtrack(target, make([]uint16, 0), true, numThreads)
	wg.Wait()

	return solutions
}

func locusAlreadyInSlice(v uint16, array []uint16) bool {
	for _, a := range array {
		if a == v || (v%pgs.NumHaplotypes == 0 && a == v+1) || (v%pgs.NumHaplotypes == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func findNextPositiveMins(values []int64) []int64 {
	mins := make([]int64, len(values))
	j := len(values) - 1
	for {
		if values[j] > 0 || j == 1 {
			mins[j] = values[j]
			j--
			break
		}
		j--
		mins[j] = 0
	}
	for i := j; i >= 0; i-- {
		if values[i] > 0 && values[i] < mins[i+1] {
			mins[i] = values[i]
		} else {
			mins[i] = mins[i+1]
		}
	}
	return mins
}

func getMaxTotal(values []int64) (int64, int64) {
	positive, negative := int64(0), int64(0)
	for _, v := range values {
		if v > 0 {
			positive += pgs.NumHaplotypes * v
		} else {
			negative += pgs.NumHaplotypes * v
		}
	}
	return positive, negative
}

func getScore(snps []uint8, weights []int64) int64 {
	var score int64 = 0
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

func findPositiveMin(values []int64) int64 {
	minV := values[0]
	j := 0
	for {
		if values[j] >= 0 {
			minV = values[j]
			break
		}
		j++
	}
	for _, v := range values {
		if v > 0 && v < minV {
			minV = v
		}
	}
	return minV
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
	for k := 0; k < params.ITERATIONS; k++ {
		if k%500 == 0 {
			fmt.Printf("Iteration %d/%d\n", k, params.ITERATIONS)
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
		if _, ok := base[ArrayToString(ext)]; !ok {
			base[ArrayToString(ext)] = make([]uint8, len(ext))
			copy(base[ArrayToString(ext)], ext)
		}
	}
}

func (s *Solver) calculateDeltas(population [][]uint8) ([]float64, [][]uint8) {
	var delta float64
	var err error
	deltas := make([]float64, len(population))
	matches := make([][]uint8, 0)
	for i := range population {
		delta = CalculateScore(population[i], s.p.Weights) - s.target
		//for delta == 0 {
		for math.Abs(delta) <= pgs.ErrorMargin {
			match := make([]uint8, len(population[i]))
			copy(match, population[i])
			matches = append(matches, match)
			//fmt.Printf("Found match: %s %f\n", ArrayToString(match), CalculateScore(match, s.p.Weights)-s.target)
			population[i], err = s.p.SampleFromPopulation()
			if err != nil {
				log.Printf("Error resampling in fitness calculation: %v\n", err)
				return nil, nil
			}
			delta = CalculateScore(population[i], s.p.Weights) - s.target
		}
		deltas[i] = delta
	}
	return deltas, matches
}

//func CalculateScore(snps []int, weights []float64) float64 {
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

func CalculateScore(snps []uint8, weights []float64) float64 {
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

//func Accuracy(solution []int, target []int) float64 {
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

func Accuracy(solution []uint8, target []uint8) float64 {
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

	mutexed := tools.NewMutexMap(solutions)
	for _, solution := range solutions {
		//fmt.Printf("Exploring %s\n", ArrayToString(solution))
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
	//			//fmt.Printf("Exploring %s\n", ArrayToString(solution))
	//			explore(solution, weightCopies, mutexed)
	//		}
	//	}()
	//}
	//wg.Wait()
	return mutexed.RetrieveMapOnly()
}

func explore(source []uint8, positions [][]int, saver *tools.MutexMap) {
	var leftPos, rightPos int
	var leftVal, rightVal uint8
	//fmt.Printf("Exploring %s\n", ArrayToString(source))
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
				saver.Put(ArrayToString(swapped), swapped)
				if len(positions) > 1 {
					explore(swapped, positions[1:], saver)
				}
				// second, splitting 2 + 0 into 1 + 1
				if (leftVal == 2 && rightVal == 0) || (leftVal == 0 && rightVal == 2) {
					split := make([]uint8, len(source))
					copy(split, source)
					split[leftPos*pgs.NumHaplotypes], split[leftPos*pgs.NumHaplotypes+1] = 1, 0
					split[rightPos*pgs.NumHaplotypes], split[rightPos*pgs.NumHaplotypes+1] = 0, 1
					saver.Put(ArrayToString(split), split)
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
				saver.Put(ArrayToString(leftPushed), leftPushed)
				rightPushed[leftPos*pgs.NumHaplotypes], rightPushed[leftPos*pgs.NumHaplotypes+1] = 0, 0
				rightPushed[rightPos*pgs.NumHaplotypes], rightPushed[rightPos*pgs.NumHaplotypes+1] = 1, 1
				saver.Put(ArrayToString(rightPushed), rightPushed)
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
	//		saver.Put(ArrayToString(distributed), distributed)
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

func ArrayToString(array []uint8) string {
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

func SortByAccuracy(solutions map[string][]uint8, target []uint8) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	accuracies := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		accuracies[i] = Accuracy(solution, target)
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
