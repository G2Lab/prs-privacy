package solver

import (
	"context"
	"fmt"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"log"
	"math"
	"sort"
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

type Beta struct {
	Pos    uint16
	Weight int64
}

func newBeta(pos int, weight int64) *Beta {
	return &Beta{
		Pos:    uint16(pos),
		Weight: weight,
	}
}

type Rounder struct {
	RoundedMode  bool
	Weights      []int64
	Target       int64
	ErrorMargin  int64
	RounderError int64
}

func (s *DP) Solve(numThreads int) map[string][]uint8 {
	var multiplier float64
	var roundingErrorLeft, roundingErrorRight int64 = 0, 0
	rounder := new(Rounder)
	if s.p.WeightPrecision > params.PrecisionsLimit {
		multiplier = math.Pow(10, params.PrecisionsLimit)
		roundingErrorLeft = int64(s.p.VariantCount)
		preciseMultiplier := math.Pow(10, float64(s.p.WeightPrecision))
		rounder.RoundedMode = true
		rounder.Target = int64(s.target * preciseMultiplier)
		rounder.Weights = make([]int64, len(s.p.Weights))
		rounder.RounderError = roundingErrorLeft
		for i, w := range s.p.Weights {
			rounder.Weights[i] = int64(w * preciseMultiplier)
		}
		if s.p.WeightPrecision > params.FloatPrecision {
			rounder.ErrorMargin = 1e3
		}
	} else {
		multiplier = math.Pow(10, float64(s.p.WeightPrecision))
	}

	target := int64(s.target * multiplier)
	betasLeft := betasFromWeights(s.p.Weights[:len(s.p.Weights)/2], 0, multiplier)
	betasRight := betasFromWeights(s.p.Weights[len(s.p.Weights)/2:], pgs.NumHaplotypes*(len(s.p.Weights)/2), multiplier)
	maxTotalPositiveLeft, maxTotalNegativeLeft := getMaxTotal(betasLeft)
	maxTotalPositiveRight, maxTotalNegativeRight := getMaxTotal(betasRight)
	if rounder.RoundedMode && (maxTotalNegativeLeft < 0 || maxTotalNegativeRight < 0) {
		roundingErrorRight = roundingErrorLeft
	}
	//fmt.Println(roundingErrorRight)

	tableLeft := calculateSubsetSumsTable(betasLeft, target-maxTotalNegativeRight, target-maxTotalPositiveRight-maxTotalPositiveLeft)
	tableRight := calculateSubsetSumsTable(betasRight, target-maxTotalNegativeLeft, target-maxTotalPositiveLeft-maxTotalPositiveRight)
	fmt.Printf("%d, %d\n", len(tableLeft), len(tableRight))

	// Combine and backtrack
	targets := []int64{target}
	if rounder.RoundedMode {
		for w := target + roundingErrorRight; w > target-roundingErrorLeft; w-- {
			if w == target {
				continue
			}
			targets = append(targets, w)
		}
	}
	fmt.Printf("Target: %d\n", target)
	subsets := make([][]uint16, 0)
	for rightSum := range tableRight {
		for _, t := range targets {
			if _, ok := tableLeft[t-rightSum]; ok {
				//fmt.Printf("\n%d, %d, %d\n", t, rightSum, t-rightSum)
				//for _, beta := range tableRight[rightSum] {
				//	fmt.Printf("[%d, %d] ", beta.Pos, beta.Weight)
				//}
				//fmt.Println()
				//for _, beta := range tableLeft[t-rightSum] {
				//	fmt.Printf("[%d, %d] ", beta.Pos, beta.Weight)
				//}
				start := make([]uint16, 0)
				rightHalfSolutions := backtrack(start, rightSum, tableRight, false, rounder)
				for _, rightHalfSolution := range rightHalfSolutions {
					subsets = append(subsets, backtrack(rightHalfSolution, t-rightSum, tableLeft, true, rounder)...)
				}
			}
		}
	}

	solutions := make(map[string][]uint8)
	for _, subset := range subsets {
		sol := lociToGenotype(subset, pgs.NumHaplotypes*len(s.p.Weights))
		solutions[ArrayToString(sol)] = sol
	}
	return solutions

	////backtracking
	//var weight int64
	//var backtrack func(int64, []uint16)
	//backtrack = func(sum int64, path []uint16) {
	//	var ok bool
	//	if int64(math.Abs(float64(sum))) <= roundingError {
	//		validateSolution(path)
	//		return
	//	}
	//	tail := make(map[int64][]uint16)
	//	tail[sum], ok = table[sum]
	//	if !ok {
	//		// Due to the loss of precision with either rounding or floating point, the sum in the table might be
	//		// slightly different.
	//		for w := sum - roundingError; w < sum+roundingError; w++ {
	//			if w == sum {
	//				continue
	//			}
	//			if pos, exists := table[w]; exists {
	//				tail[w] = pos
	//			}
	//		}
	//	}
	//	for w, pts := range tail {
	//		for _, p := range pts {
	//			// Make sure that we do not have paths with two values for the same snp
	//			if locusAlreadyExists(p, path) || (len(path) > 0 && p > path[len(path)-1]) {
	//				continue
	//			}
	//			newPath := make([]uint16, len(path)+1)
	//			copy(newPath, path)
	//			newPath[len(path)] = p
	//			weight = weights[p/2] * (1 + int64(p%2))
	//			backtrack(w-weight, newPath)
	//		}
	//	}
	//}

	//var wg sync.WaitGroup
	//type state struct {
	//	sum       int64
	//	positions []uint16
	//}
	//taskChan := make(chan state)
	//traverser := func(tasks <-chan state) {
	//	for task := range tasks {
	//		backtrack(task.sum, task.positions)
	//		wg.Done()
	//	}
	//}
	//
	//for i := 0; i < numThreads; i++ {
	//	go traverser(taskChan)
	//}
	//
	//start := make(map[int64][]uint16)
	//count := 0
	//pointers, exists := table[target]
	//if exists {
	//	start[target] = pointers
	//	count += len(pointers)
	//}
	//if !exists || roundedMode {
	//	for w := targetLeft; w < targetRight; w++ {
	//		//for w := target - roundingError; w < target; w++ {
	//		if w == target {
	//			continue
	//		}
	//		if pos, exists := table[w]; exists {
	//			start[w] = pos
	//			count += len(pos)
	//		}
	//	}
	//}
	//go func() {
	//	if target != 0 {
	//		wg.Add(count)
	//		for w, pts := range start {
	//			for _, p := range pts {
	//				taskChan <- state{sum: w - weights[p/2]*(1+int64(p%2)), positions: []uint16{p}}
	//			}
	//		}
	//	} else {
	//		wg.Add(1)
	//		taskChan <- state{sum: target, positions: []uint16{}}
	//	}
	//	wg.Wait()
	//	close(results)
	//	close(taskChan)
	//}()
	//
	//solutions := make(map[string][]uint8)
	//for res := range results {
	//	//sf, err := os.OpenFile(fmt.Sprintf("solutions_%s.txt", s.p.PgsID), os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	//	//if err != nil {
	//	//	log.Printf("Error opening solutions file: %v\n", err)
	//	//	return nil
	//	//}
	//	//defer sf.Close()
	//	//_, err = sf.WriteString(fmt.Sprintf("%s\n", ArrayToString(sol)))
	//	//if err != nil {
	//	//	log.Printf("Error writing %s to solutions file: %v\n", ArrayToString(sol), err)
	//	//}
	//	solutions[ArrayToString(res)] = res
	//}
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumsTable(betas []*Beta, upperBound, lowerBound int64) map[int64][]*Beta {
	// Fill out the table using dynamic programming
	table := make(map[int64][]*Beta)
	// add the zero weight
	table[0] = make([]*Beta, 0)
	var prevSum, nextSum int64
	for i := 0; i < len(betas); i += pgs.NumHaplotypes {
		//fmt.Printf("Position %d/%d\n", i+1, len(betas))
		if betas[i].Weight > 0 {
			lowerBound += betas[i].Weight
		}
		existingSums := make([]int64, 0, len(table))
		for s := range table {
			existingSums = append(existingSums, s)
		}
		for k := 0; k < pgs.NumHaplotypes; k++ {
			for _, prevSum = range existingSums {
				nextSum = prevSum + betas[i+k].Weight
				if nextSum <= upperBound && nextSum >= lowerBound {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]*Beta, 0)
					}
					table[nextSum] = append(table[nextSum], betas[i+k])
				}
			}
		}
	}
	return table
}

func backtrack(path []uint16, sum int64, table map[int64][]*Beta, validation bool, rounder *Rounder) [][]uint16 {
	if int64(math.Abs(float64(sum))) <= rounder.RounderError {
		if validation && rounder.RoundedMode {
			preciseSum := lociToScore(path, rounder.Weights)
			if int64(math.Abs(float64(preciseSum-rounder.Target))) > rounder.ErrorMargin {
				return nil
			}
		}
		return [][]uint16{path}
	}
	output := make([][]uint16, 0)
	for _, beta := range table[sum] {
		if locusAlreadyExists(beta.Pos, path) {
			continue
		}
		newState := make([]uint16, len(path)+1)
		copy(newState, path)
		newState[len(path)] = beta.Pos
		if res := backtrack(newState, sum-beta.Weight, table, validation, rounder); res != nil {
			output = append(output, res...)
		}
	}
	return output
}

func locusAlreadyExists(v uint16, array []uint16) bool {
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

func getMaxTotal(betas []*Beta) (int64, int64) {
	lower, upper := int64(0), int64(0)
	// Consider the double-allele weights to find the maximum possible
	for i := 1; i < len(betas); i += 2 {
		if betas[i].Weight > 0 {
			lower += betas[i].Weight
		} else {
			upper += betas[i].Weight
		}
	}
	return lower, upper
}

func betasFromWeights(weights []float64, startIdx int, multiplier float64) []*Beta {
	indices := sortedIndices(weights)
	betas := make([]*Beta, len(weights)*pgs.NumHaplotypes)
	for i, pos := range indices {
		for j := 0; j < pgs.NumHaplotypes; j++ {
			betas[pgs.NumHaplotypes*i+j] = newBeta(startIdx+pgs.NumHaplotypes*pos+j, int64(j+1)*int64(weights[pos]*multiplier))
		}
	}
	return betas
}

func betasToScore(betas []*Beta, weights []int64) int64 {
	var score int64 = 0
	for _, beta := range betas {
		switch beta.Pos % pgs.NumHaplotypes {
		case 0:
			score += weights[beta.Pos/2]
		case 1:
			score += weights[beta.Pos/2] * pgs.NumHaplotypes
		}
	}
	return score
}

func lociToScore(loci []uint16, weights []int64) int64 {
	var score int64 = 0
	for _, locus := range loci {
		switch locus % pgs.NumHaplotypes {
		case 0:
			score += weights[locus/2]
		case 1:
			score += weights[locus/2] * pgs.NumHaplotypes
		}
	}
	return score
}

func lociToGenotype(loci []uint16, total int) []uint8 {
	sol := make([]uint8, total)
	for _, pos := range loci {
		sol[pos] = 1
		if pos%pgs.NumHaplotypes == 1 {
			sol[pos-1] = 1
		}
	}
	return sol
}

func sortedIndices(values []float64) []int {
	// Create a slice of indices.
	indices := make([]int, len(values))
	for i := range indices {
		indices[i] = i
	}

	// Sort the indices based on the values.
	sort.Slice(indices, func(i, j int) bool {
		return values[indices[i]] < values[indices[j]]
	})

	return indices
}

//func getMaxTotal(values []int64) (int64, int64) {
//	positive, negative := int64(0), int64(0)
//	for _, v := range values {
//		if v > 0 {
//			positive += pgs.NumHaplotypes * v
//		} else {
//			negative += pgs.NumHaplotypes * v
//		}
//	}
//	return positive, negative
//}

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
