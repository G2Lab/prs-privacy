package solver

import (
	"container/heap"
	"context"
	"fmt"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
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

type Rounder struct {
	RoundedMode  bool
	RounderError int64
}

func (s *DP) Solve(numThreads int) map[string][]uint8 {
	var roundingErrorLeft, roundingErrorRight int64 = 0, 0
	rounder := new(Rounder)
	multiplier := new(big.Rat).SetFloat64(math.Pow(10, float64(s.p.WeightPrecision)))
	if s.p.WeightPrecision > params.PrecisionsLimit {
		multiplier.SetFloat64(math.Pow(10, params.PrecisionsLimit))
		//roundingErrorLeft = int64(s.p.VariantCount) * 5 / 4
		roundingErrorLeft = int64(s.p.VariantCount)
		rounder.RoundedMode = true
		rounder.RounderError = roundingErrorLeft
	}

	tf, _ := new(big.Rat).Mul(s.target, multiplier).Float64()
	target := int64(tf)
	allBetas := betasFromWeights(s.p.Weights, multiplier)
	maxTotalPositive, maxTotalNegative := getMaxTotal(allBetas)
	if rounder.RoundedMode && maxTotalNegative < 0 {
		roundingErrorRight = roundingErrorLeft
	}
	//fmt.Printf("Precision: %d\n", s.p.WeightPrecision)
	//fmt.Printf("Target: %d\n", target)

	numSegments := 4
	splitIdxHalf := len(s.p.Weights) / 2
	betas := make([]map[uint16]int64, numSegments)
	tables := make([]map[int64][]uint16, numSegments)
	betas[0] = makeBetaMap(allBetas, 0, splitIdxHalf/2)
	betas[1] = makeBetaMap(allBetas, splitIdxHalf/2, splitIdxHalf)
	betas[2] = makeBetaMap(allBetas, splitIdxHalf, splitIdxHalf+(len(allBetas)-splitIdxHalf)/2)
	betas[3] = makeBetaMap(allBetas, splitIdxHalf+(len(allBetas)-splitIdxHalf)/2, len(allBetas))
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumsTable(betas[i], target-maxTotalNegative+roundingErrorLeft, target-maxTotalPositive-roundingErrorRight)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
	}

	var modulo int64 = 0
	if math.Abs(float64(target)) >= math.Pow(10, 7) {
		modulo = tools.FindNextBiggerPrime(int64(math.Sqrt(math.Abs(float64(target)))))
	} else {
		modulo = tools.FindNextSmallerPrime(int64(math.Abs(float64(target))))
	}
	fmt.Printf("Modulo: %d\n", modulo)
	moduloMaps := make([]map[int32][]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		moduloMaps[i] = buildModuloMap(modulo, tables[i])
	}

	// If we have rounded the target, search for other targets nearby
	targets := []int64{target}
	if rounder.RoundedMode {
		for w := target + roundingErrorRight; w > target-roundingErrorLeft; w-- {
			if w == target {
				continue
			}
			targets = append(targets, w)
		}
	}

	// Combine and backtrack
	umodulo := int32(modulo)
	subsets := make([][]uint16, 0)
	var midValue, modTarget int32
	for midValue = 0; midValue < umodulo; midValue++ {
		if midValue%50 == 0 {
			fmt.Printf("MidValue: %d\n", midValue)
		}
		combL := allSumCombinations(midValue, umodulo, moduloMaps[:numSegments/2], tables[:numSegments/2],
			betas[:numSegments/2], rounder)
		// No pair adds up to this midValue
		if len(combL) == 0 {
			continue
		}
		for _, t := range targets {
			modTarget = tools.SubMod(int32(tools.Mod(t, modulo)), midValue, umodulo)
			combR := allSumCombinations(modTarget, umodulo, moduloMaps[numSegments/2:], tables[numSegments/2:],
				betas[numSegments/2:], rounder)
			for sumL, solsL := range combL {
				for sumR, solsR := range combR {
					if sumR+sumL != t {
						continue
					}
					//fmt.Printf("Found a solution!\n")
					//fmt.Println(sumL, sumR, t)
					//fmt.Println(solsL)
					//fmt.Println(solsR)
					//fmt.Println("//////////")
					for j := range solsL {
						for i := range solsR {
							joint := make([]uint16, len(solsL[j])+len(solsR[i]))
							copy(joint, solsL[j])
							copy(joint[len(solsL[j]):], solsR[i])
							//joint := append(solsL[j], solsR[i]...)
							if rounder.RoundedMode {
								preciseSum := lociToScore(joint, s.p.Weights)
								if preciseSum.Cmp(s.target) != 0 {
									continue
								}
							}
							//fmt.Println(combL)
							//fmt.Println(combR)
							//fmt.Println("-----", ArrayToString(lociToGenotype(joint, pgs.NumHaplotypes*len(s.p.Weights))),
							//	new(big.Rat).Sub(s.target, lociToScore(joint, s.p.Weights)).FloatString(12), sumR+sumL)
							subsets = append(subsets, joint)
						}
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

	//lowestAbsoluteLikelihood := math.MaxFloat64
	//var lkl float64
	//lkl = calculateNegativeLikelihood(joint, len(s.p.Weights)*pgs.NumHaplotypes, s.p)
	//if lkl <= lowestAbsoluteLikelihood {
	//	subsets = append(subsets, joint)
	//	lowestAbsoluteLikelihood = lkl
	//}
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumsTable(betas map[uint16]int64, upperBound, lowerBound int64) map[int64][]uint16 {
	// Fill out the table using dynamic programming
	table := make(map[int64][]uint16)
	// add the zero weight
	table[0] = make([]uint16, 0)
	indices := make([]uint16, 0, len(betas))
	for i := range betas {
		indices = append(indices, i)
	}
	sort.Slice(indices, func(i, j int) bool {
		return indices[i] < indices[j]
	})
	i := 0
	var prevSum, nextSum int64
	var k uint16
	for _, pos := range indices {
		i++
		fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[pos] > 0 {
			lowerBound += pgs.NumHaplotypes * betas[pos]
		} else {
			upperBound += pgs.NumHaplotypes * betas[pos]
		}
		existingSums := make([]int64, 0, len(table))
		for s := range table {
			existingSums = append(existingSums, s)
		}
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.NumHaplotypes; k++ {
				nextSum = prevSum + int64(k)*betas[pos]
				if nextSum >= lowerBound && nextSum <= upperBound {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
					}
					table[nextSum] = append(table[nextSum], pgs.NumHaplotypes*pos+k-1)
				}
			}
		}
	}
	return table
}

func backtrack(path []uint16, sum int64, table map[int64][]uint16, weights map[uint16]int64, rounder *Rounder) [][]uint16 {
	if int64(math.Abs(float64(sum))) <= rounder.RounderError {
		return [][]uint16{path}
	}
	output := make([][]uint16, 0)
	//newStates := make([][]uint16, 0)
	//newSums := make([]int64, 0)
	var weight int64
	for _, ptr := range table[sum] {
		if locusAlreadyExists(ptr, path) || (len(path) > 0 && ptr > path[len(path)-1]) {
			continue
		}
		//newStates = append(newStates, append(path, ptr))
		//newSums = append(newSums, sum-weights[ptr/2]*(1+int64(ptr%2)))
		//newStates = append(newStates, make([]uint16, len(path)+1))
		//copy(newStates[len(newStates)-1], path)
		//newStates[len(newStates)-1][len(path)] = ptr
		//
		//
		newState := make([]uint16, len(path)+1)
		copy(newState, path)
		newState[len(path)] = ptr
		weight = weights[ptr/2] * (1 + int64(ptr%2))
		if res := backtrack(newState, sum-weight, table, weights, rounder); res != nil {
			output = append(output, res...)
		}
	}
	//path = nil
	//for i := range newStates {
	//	if res := backtrack(newStates[i], newSums[i], table, weights, rounder); res != nil {
	//		output = append(output, res...)
	//	}
	//}
	return output
}

func buildModuloMap(modulo int64, table map[int64][]uint16) map[int32][]int64 {
	moduloMap := make(map[int32][]int64)
	var reduced int32
	for sum := range table {
		reduced = int32(tools.Mod(sum, modulo))
		if _, ok := moduloMap[reduced]; !ok {
			moduloMap[reduced] = make([]int64, 0)
		}
		moduloMap[reduced] = append(moduloMap[reduced], sum)
	}
	return moduloMap
}

func allSumCombinations(modSum, umodulo int32, modTables []map[int32][]int64, fullTables []map[int64][]uint16,
	betas []map[uint16]int64, rounder *Rounder) map[int64][][]uint16 {
	combinations := make(map[int64][][]uint16)
	var valueEntry, valueExit int32
	var tmpSum int64
	for valueEntry = range modTables[0] {
		valueExit = tools.SubMod(modSum, valueEntry, umodulo)
		if _, ok := modTables[1][valueExit]; !ok {
			continue
		}
		leftSolComb := backtrackedSolutions(modTables[0][valueEntry], fullTables[0], betas[0], rounder)
		rightSolComb := backtrackedSolutions(modTables[1][valueExit], fullTables[1], betas[1], rounder)
		for nonModSumL, lSols := range leftSolComb {
			for nonModSumR, rSols := range rightSolComb {
				tmpSum = nonModSumL + nonModSumR
				for i := range lSols {
					for j := range rSols {
						if _, ok := combinations[tmpSum]; !ok {
							combinations[tmpSum] = make([][]uint16, 0)
						}
						combinations[tmpSum] = append(combinations[tmpSum], append(lSols[i], rSols[j]...))
						//joint := make([]uint16, len(lSols[i])+len(rSols[j]))
						//copy(joint, lSols[i])
						//copy(joint[len(lSols[i]):], rSols[j])
						//combinations[tmpSum] = append(combinations[tmpSum], joint)
					}
				}
			}
		}
	}
	return combinations
}

func backtrackedSolutions(sums []int64, table map[int64][]uint16, betas map[uint16]int64, rounder *Rounder) map[int64][][]uint16 {
	solutions := make(map[int64][][]uint16)
	for _, sum := range sums {
		input := make([]uint16, 0)
		solutions[sum] = backtrack(input, sum, table, betas, rounder)
	}
	return solutions
}

func makeBetaMap(betas []int64, start, end int) map[uint16]int64 {
	bmap := make(map[uint16]int64)
	for i := start; i < end; i++ {
		bmap[uint16(i)] = betas[i]
	}
	return bmap
}

func sortByLikelihood(candidates [][]uint16, totalLen int, p *pgs.PGS) [][]uint16 {
	likelihoods := make([]float64, len(candidates))
	for i, candidate := range candidates {
		likelihoods[i] = calculateNegativeLikelihood(candidate, totalLen, p)
	}
	sortBy(candidates, likelihoods)
	return candidates
}

type candidateWithLikelihood struct {
	candidate  []uint16
	likelihood float64
}

type likelihoodHeap []candidateWithLikelihood

func (lh likelihoodHeap) Len() int { return len(lh) }
func (lh likelihoodHeap) Less(i, j int) bool {
	return lh[i].likelihood > lh[j].likelihood
}
func (lh likelihoodHeap) Swap(i, j int) {
	lh[i].candidate, lh[j].candidate = lh[j].candidate, lh[i].candidate
	lh[i].likelihood, lh[j].likelihood = lh[j].likelihood, lh[i].likelihood
}

func (lh *likelihoodHeap) Push(x interface{}) {
	*lh = append(*lh, x.(candidateWithLikelihood))
}

func (lh *likelihoodHeap) Pop() interface{} {
	old := *lh
	n := len(old)
	x := old[n-1]
	*lh = old[0 : n-1]
	return x
}

func selectTopLikelihoodCandidates(candidates [][]uint16, totalLen int, p *pgs.PGS, top int) [][]uint16 {
	if len(candidates) <= top {
		return candidates
	}
	h := &likelihoodHeap{}
	// Push the first N slices onto the heap
	for i := 0; i < top; i++ {
		heap.Push(h, candidateWithLikelihood{candidates[i], calculateNegativeLikelihood(candidates[i], totalLen, p)})
	}
	var lkl float64
	// Update the heap with remaining slices
	for i := top; i < len(candidates); i++ {
		lkl = calculateNegativeLikelihood(candidates[i], totalLen, p)
		if lkl > (*h)[0].likelihood {
			// Replace the minimum element with the current slice
			heap.Pop(h)
			heap.Push(h, candidateWithLikelihood{candidates[i], lkl})
		}
	}

	// Extract slices from the heap
	result := make([][]uint16, top)
	for i := top - 1; i >= 0; i-- {
		result[i] = heap.Pop(h).(candidateWithLikelihood).candidate
	}

	return result
}

func calculateNegativeLikelihood(sequence []uint16, totalLen int, p *pgs.PGS) float64 {
	var likelihood float64 = 0
	indexed := make(map[uint16]struct{})
	for _, pos := range sequence {
		indexed[pos] = struct{}{}
	}
	for j := 0; j < totalLen; j++ {
		if _, ok := indexed[uint16(j)]; ok {
			likelihood += mafToLikelihood(p.Maf[j/2][1])
		} else {
			likelihood += mafToLikelihood(p.Maf[j/2][0])
		}
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

func betasFromWeights(weights []*big.Rat, multiplier *big.Rat) []int64 {
	betas := make([]int64, len(weights))
	for i := range betas {
		tmp, _ := new(big.Rat).Mul(weights[i], multiplier).Float64()
		betas[i] = int64(tmp)
	}
	return betas
}

func lociToScore(loci []uint16, weights []*big.Rat) *big.Rat {
	score := new(big.Rat).SetInt64(0)
	for _, locus := range loci {
		score.Add(score, weights[locus/2])
		if locus%2 == 1 {
			score.Add(score, weights[locus/2])
		}
	}
	return score
}

func lociToGenotype(loci []uint16, total int) []uint8 {
	sol := make([]uint8, total)
	for _, locus := range loci {
		sol[locus] = 1
		if locus%2 == 1 {
			sol[locus-1] = 1
		}
		//if sol[pgs.NumHaplotypes*locus] == 0 {
		//	sol[pgs.NumHaplotypes*locus] = 1
		//} else {
		//	sol[pgs.NumHaplotypes*locus+1] = 1
		//}
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
