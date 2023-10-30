package solver

import (
	"context"
	"fmt"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"log"
	"math"
	"sync"
)

type DP struct {
	target float64
	p      *pgs.PGS
}

func NewDP(ctx context.Context, target float64, p *pgs.PGS, numThreads int) *DP {
	s := &DP{
		target: target,
		p:      p,
	}
	return s
}

func (s *DP) Solve(numThreads int) map[string][]uint8 {
	var roundedMode = false
	var multiplier float64
	var roundingError int64 = 0
	var errorMargin int64 = 1
	var preciseTarget int64
	var preciseWeights []int64
	var preciseMultiplier float64
	if s.p.WeightPrecision > params.PrecisionsLimit {
		multiplier = math.Pow(10, params.PrecisionsLimit)
		preciseMultiplier = math.Pow(10, float64(s.p.WeightPrecision))
		preciseTarget = int64(s.target * preciseMultiplier)
		//fmt.Printf("Precise target: %d\n", preciseTarget)
		preciseWeights = make([]int64, len(s.p.Weights))
		for i, w := range s.p.Weights {
			preciseWeights[i] = int64(w * preciseMultiplier)
		}
		roundedMode = true
		roundingError = int64(s.p.VariantCount)
		if s.p.WeightPrecision > params.FloatPrecision {
			errorMargin = 1e3
		}
	} else {
		multiplier = math.Pow(10, float64(s.p.WeightPrecision))
	}

	target := int64(s.target * multiplier)
	weights := make([]int64, len(s.p.Weights))
	for i, w := range s.p.Weights {
		weights[i] = int64(w * multiplier)
	}
	maxTotalPositive, maxTotalNegative := getMaxTotal(weights)
	targetLeft := target - roundingError
	var targetRight int64
	var remainingMinValues []int64
	if maxTotalNegative < 0 {
		targetRight = target + roundingError
		remainingMinValues = make([]int64, len(weights))
	} else {
		targetRight = target
		remainingMinValues = findNextPositiveMins(weights)
	}
	upperBound := targetRight - maxTotalNegative
	lowerBound := targetLeft - maxTotalPositive

	//fmt.Printf("Target: %d±%d\n", target, roundingError)
	//fmt.Printf("Weights: %v\n", weights)
	//fmt.Printf("Upper bound: %d\n", maxTotalPositive)
	//fmt.Printf("Lower bound: %d\n", maxTotalNegative)
	//fmt.Printf("Next mins: %v\n", remainingMinValues)

	// Fill the Solve table using dynamic programming
	table := make(map[int64][]uint16)
	// add the zero weight
	table[0] = make([]uint16, 0)
	var prevSum, nextSum int64
	for i := 0; i < len(weights); i++ {
		//fmt.Printf("Position %d/%d\n", i+1, len(weights))
		if weights[i] > 0 {
			lowerBound += pgs.NumHaplotypes * weights[i]
		} else {
			upperBound += pgs.NumHaplotypes * weights[i]
		}
		//fmt.Printf("%d) upper: %d, lower: %d, min: %d\n", i, upperBound, lowerBound, remainingMinValues[i])
		existingSums := make([]int64, 0, len(table))
		for w := range table {
			existingSums = append(existingSums, w)
		}
		for _, prevSum = range existingSums {
			for k := 1; k <= pgs.NumHaplotypes; k++ {
				nextSum = prevSum + int64(k)*weights[i]
				if (nextSum <= upperBound-remainingMinValues[i] && nextSum >= lowerBound) || (nextSum >= targetLeft && nextSum <= targetRight) {
					//if (nextSum >= lowerBound && nextSum <= upperBound) || (nextSum >= targetLeft && nextSum <= targetRight) {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
					}
					table[nextSum] = append(table[nextSum], uint16(pgs.NumHaplotypes*i+k-1))
				}
			}
		}
	}
	//fmt.Printf("Number of weights in the table %d\n", len(table))

	results := make(chan []uint8)
	validateSolution := func(path []uint16) {
		sol := make([]uint8, len(s.p.Weights)*pgs.NumHaplotypes)
		for _, pos := range path {
			sol[pos] = 1
			if pos%pgs.NumHaplotypes == 1 {
				sol[pos-1] = 1
			}
		}
		//// If we are in the roundedMode mode, check that the precise weights add up to the target
		if roundedMode && int64(math.Abs(float64(getScore(sol, preciseWeights)-preciseTarget))) > errorMargin {
			return
		}
		results <- sol
	}

	//backtracking
	var weight int64
	var backtrack func(int64, []uint16)
	backtrack = func(sum int64, path []uint16) {
		var ok bool
		if int64(math.Abs(float64(sum))) <= roundingError {
			validateSolution(path)
			return
		}
		tail := make(map[int64][]uint16)
		tail[sum], ok = table[sum]
		if !ok {
			// Due to the loss of precision with either rounding or floating point, the sum in the table might be
			// slightly different.
			for w := sum - roundingError; w < sum+roundingError; w++ {
				if w == sum {
					continue
				}
				if pos, exists := table[w]; exists {
					tail[w] = pos
				}
			}
		}
		for w, pts := range tail {
			for _, p := range pts {
				// Make sure that we do not have paths with two values for the same snp
				if locusAlreadyInSlice(p, path) || (len(path) > 0 && p > path[len(path)-1]) {
					continue
				}
				newPath := make([]uint16, len(path)+1)
				copy(newPath, path)
				newPath[len(path)] = p
				weight = weights[p/2] * (1 + int64(p%2))
				backtrack(w-weight, newPath)
			}
		}
	}

	var wg sync.WaitGroup
	type state struct {
		sum       int64
		positions []uint16
	}
	taskChan := make(chan state)
	traverser := func(tasks <-chan state) {
		for task := range tasks {
			backtrack(task.sum, task.positions)
			wg.Done()
		}
	}

	for i := 0; i < numThreads; i++ {
		go traverser(taskChan)
	}

	start := make(map[int64][]uint16)
	count := 0
	pointers, exists := table[target]
	if exists {
		start[target] = pointers
		count += len(pointers)
	}
	if !exists || roundedMode {
		for w := targetLeft; w < targetRight; w++ {
			//for w := target - roundingError; w < target; w++ {
			if w == target {
				continue
			}
			if pos, exists := table[w]; exists {
				start[w] = pos
				count += len(pos)
			}
		}
	}
	fmt.Println(start)
	go func() {
		if target != 0 {
			wg.Add(count)
			for w, pts := range start {
				for _, p := range pts {
					taskChan <- state{sum: w - weights[p/2]*(1+int64(p%2)), positions: []uint16{p}}
				}
			}
		} else {
			wg.Add(1)
			taskChan <- state{sum: target, positions: []uint16{}}
		}
		wg.Wait()
		close(results)
		close(taskChan)
	}()

	solutions := make(map[string][]uint8)
	for res := range results {
		//sf, err := os.OpenFile(fmt.Sprintf("solutions_%s.txt", s.p.PgsID), os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
		//if err != nil {
		//	log.Printf("Error opening solutions file: %v\n", err)
		//	return nil
		//}
		//defer sf.Close()
		//_, err = sf.WriteString(fmt.Sprintf("%s\n", ArrayToString(sol)))
		//if err != nil {
		//	log.Printf("Error writing %s to solutions file: %v\n", ArrayToString(sol), err)
		//}
		solutions[ArrayToString(res)] = res
	}

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
