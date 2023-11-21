package solver

import (
	"context"
	"fmt"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"log"
	"math"
	"math/big"
	"sort"
)

type DP struct {
	target *big.Rat
	p      *pgs.PGS
}

func NewDP(ctx context.Context, target *big.Rat, p *pgs.PGS, numThreads int) *DP {
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
	RounderError int64
}

func (s *DP) Solve(numThreads int) map[string][]uint8 {
	var roundingErrorLeft, roundingErrorRight int64 = 0, 0
	rounder := new(Rounder)
	//var multiplier float64
	multiplier := new(big.Rat).SetFloat64(math.Pow(10, float64(s.p.WeightPrecision)))
	if s.p.WeightPrecision > params.PrecisionsLimit {
		multiplier.SetFloat64(math.Pow(10, params.PrecisionsLimit))
		roundingErrorLeft = int64(s.p.VariantCount)
		rounder.RoundedMode = true
		rounder.RounderError = roundingErrorLeft
	}

	tf, _ := new(big.Rat).Mul(s.target, multiplier).Float64()
	target := int64(tf)
	betasLeft := betasFromWeights(s.p.Weights[:len(s.p.Weights)/2], 0, multiplier)
	betasRight := betasFromWeights(s.p.Weights[len(s.p.Weights)/2:], pgs.NumHaplotypes*(len(s.p.Weights)/2), multiplier)
	maxTotalPositiveLeft, maxTotalNegativeLeft := getMaxTotal(betasLeft)
	maxTotalPositiveRight, maxTotalNegativeRight := getMaxTotal(betasRight)
	if rounder.RoundedMode && (maxTotalNegativeLeft < 0 || maxTotalNegativeRight < 0) {
		roundingErrorRight = roundingErrorLeft
	}
	fmt.Printf("Precision: %d\n", s.p.WeightPrecision)
	//fmt.Printf("Multiplier: %f, \n", multiplier.Float64())
	fmt.Printf("Target: %d\n", target)

	tableLeft := calculateSubsetSumsTable(betasLeft, target-maxTotalNegativeRight-maxTotalNegativeLeft,
		target-maxTotalPositiveRight-maxTotalPositiveLeft)
	tableRight := calculateSubsetSumsTable(betasRight, target-maxTotalNegativeLeft-maxTotalNegativeRight,
		target-maxTotalPositiveLeft-maxTotalPositiveRight)
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
	//lowestAbsoluteLikelihood := math.MaxFloat64
	//var lkl float64
	subsets := make([][]uint16, 0)
	for rightSum := range tableRight {
		for _, t := range targets {
			if _, ok := tableLeft[t-rightSum]; ok {
				right, left := make([]uint16, 0), make([]uint16, 0)
				rightHalfSolutions := backtrack(right, rightSum, tableRight, rounder)
				leftHalfSolutions := backtrack(left, t-rightSum, tableLeft, rounder)
				for _, rightHalfSolution := range rightHalfSolutions {
					for _, leftHalfSolution := range leftHalfSolutions {
						joint := append(leftHalfSolution, rightHalfSolution...)
						if rounder.RoundedMode {
							preciseSum := lociToScore(joint, s.p.Weights)
							if preciseSum.Cmp(s.target) != 0 {
								continue
							}
						}
						subsets = append(subsets, joint)
						//lkl = calculateNegativeLikelihood(joint, len(s.p.Weights)*pgs.NumHaplotypes, s.p)
						//if lkl <= lowestAbsoluteLikelihood {
						//	subsets = append(subsets, joint)
						//	lowestAbsoluteLikelihood = lkl
						//}
					}
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
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumsTable(betas []*Beta, upperBound, lowerBound int64) map[int64][]*Beta {
	// Fill out the table using dynamic programming
	table := make(map[int64][]*Beta)
	// add the zero weight
	table[0] = make([]*Beta, 0)
	var prevSum, nextSum int64
	for i := 0; i < len(betas); i += pgs.NumHaplotypes {
		fmt.Printf("Position %d/%d\n", i/2+1, len(betas)/2)
		if betas[i].Weight > 0 {
			lowerBound += betas[i].Weight
		} else {
			upperBound += betas[i].Weight
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

func backtrack(path []uint16, sum int64, table map[int64][]*Beta, rounder *Rounder) [][]uint16 {
	if int64(math.Abs(float64(sum))) <= rounder.RounderError {
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
		if res := backtrack(newState, sum-beta.Weight, table, rounder); res != nil {
			output = append(output, res...)
		}
	}
	return output
}

func selectWithLowestNegativeLikelihood(candidates [][]uint16, totalLen int, p *pgs.PGS) [][]uint16 {
	likelihoods := make([]float64, len(candidates))
	for i, candidate := range candidates {
		likelihoods[i] = calculateNegativeLikelihood(candidate, totalLen, p)
	}
	sortBy(candidates, likelihoods)
	return candidates[len(candidates)*2/3:]
}

func calculateNegativeLikelihood(sequence []uint16, totalLen int, p *pgs.PGS) float64 {
	var likelihood float64 = 0
	indexed := make(map[uint16]struct{})
	for _, pos := range sequence {
		indexed[pos] = struct{}{}
	}
	for j := 0; j < totalLen; j += pgs.NumHaplotypes {
		if _, ok := indexed[uint16(j)]; ok {
			likelihood += mafToLikelihood(p.Maf[j/2][0]) + mafToLikelihood(p.Maf[j/2][1])
			continue
		}
		if _, ok := indexed[uint16(j+1)]; ok {
			likelihood += pgs.NumHaplotypes * mafToLikelihood(p.Maf[j/2][1])
			continue
		}
		likelihood += pgs.NumHaplotypes * mafToLikelihood(p.Maf[j/2][0])
	}
	return likelihood
}

func mafToLikelihood(maf float64) float64 {
	return -math.Log(maf)
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

func betasFromWeights(weights []*big.Rat, startIdx int, multiplier *big.Rat) []*Beta {
	indices := sortedIndices(weights)
	betas := make([]*Beta, len(weights)*pgs.NumHaplotypes)
	for i, pos := range indices {
		for j := 0; j < pgs.NumHaplotypes; j++ {
			tmp, _ := new(big.Rat).Mul(weights[pos], multiplier).Float64()
			betas[pgs.NumHaplotypes*i+j] = newBeta(startIdx+pgs.NumHaplotypes*pos+j, int64(j+1)*int64(tmp))
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

func lociToScore(loci []uint16, weights []*big.Rat) *big.Rat {
	score := new(big.Rat).SetInt64(0)
	for _, locus := range loci {
		score.Add(score, weights[locus/2])
		if locus%pgs.NumHaplotypes == 1 {
			score.Add(score, weights[locus/2])
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

func sortedIndices(values []*big.Rat) []int {
	// Create a slice of indices.
	indices := make([]int, len(values))
	for i := range indices {
		indices[i] = i
	}

	// Sort the indices based on the values.
	sort.Slice(indices, func(i, j int) bool {
		if values[indices[i]].Cmp(values[indices[j]]) < 0 {
			return true
		} else {
			return false
		}
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
