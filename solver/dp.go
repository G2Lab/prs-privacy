package solver

import (
	"container/heap"
	"fmt"
	"log"
	"math"
	"math/big"
	"sort"

	"github.com/ericlagergren/decimal"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

type DP struct {
	target *decimal.Big
	p      *pgs.PGS
}

func NewDP(target *decimal.Big, p *pgs.PGS, numThreads int) *DP {
	s := &DP{
		target: target,
		p:      p,
	}
	return s
}

func (s *DP) Solve(numThreads int) map[string][]uint8 {
	scale := s.target.Scale()
	multiplier := decimal.WithContext(s.p.Context).SetMantScale(1, -scale)
	fmt.Println("Scale:", scale)
	fmt.Println("Multiplier:", multiplier.String())
	scaledWeights := scaleWeights(s.p.Context, s.p.Weights, multiplier)
	scaledTarget := decimal.WithContext(s.p.Context)
	scaledTarget.Copy(s.target)
	scaledTarget.Mul(scaledTarget, multiplier)
	fmt.Println("Scaled target:", scaledTarget.String())
	fmt.Println("Scaled weights:", scaledWeights)

	numSegments := 4
	splitIdxHalf := len(scaledWeights) / 2
	betas := make([]map[uint16]*decimal.Big, numSegments)
	betas[0] = makeBetaMap(scaledWeights, 0, splitIdxHalf/2)
	betas[1] = makeBetaMap(scaledWeights, splitIdxHalf/2, splitIdxHalf)
	betas[2] = makeBetaMap(scaledWeights, splitIdxHalf, splitIdxHalf+(len(scaledWeights)-splitIdxHalf)/2)
	betas[3] = makeBetaMap(scaledWeights, splitIdxHalf+(len(scaledWeights)-splitIdxHalf)/2, len(scaledWeights))

	tables := make([]map[string][]uint16, numSegments)
	maxTotalPositive, maxTotalNegative := getMaxTotal(s.p.Context, scaledWeights)
	upper, lower := decimal.WithContext(s.p.Context), decimal.WithContext(s.p.Context)
	s.p.Context.Sub(upper, scaledTarget, maxTotalNegative)
	s.p.Context.Sub(lower, scaledTarget, maxTotalPositive)
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumsTable(s.p.Context, betas[i], upper, lower)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		//fmt.Println(tables[i])
		//fmt.Println(betas[i])
	}

	var umodulo int32
	if len(s.p.Weights)/3 < params.MaxNumModuloBits {
		umodulo = int32(tools.FindNextBiggerPrime(uint64(math.Pow(2, float64(len(s.p.Weights)/3)))))
	} else {
		umodulo = params.DefaultPrimeModulo
	}
	modulo := decimal.WithContext(s.p.Context)
	modulo.SetUint64(uint64(umodulo))
	fmt.Printf("Modulo: %s\n", modulo.String())
	moduloMaps := make([]map[int32][]string, numSegments)
	for i := 0; i < numSegments; i++ {
		moduloMaps[i] = buildModuloMap(s.p.Context, modulo, tables[i])
		//fmt.Println(moduloMaps[i])
	}

	// Combine and backtrack
	subsets := make([][]uint16, 0)
	tmp, ok := tools.BigMod(s.p.Context, scaledTarget, modulo).Int64()
	if !ok {
		log.Fatalf("Failed to convert %s to int64", scaledTarget.String())
	}
	modTarget := int32(tmp)
	fmt.Printf("Mod target: %d\n", modTarget)
	sl, sr, slr := decimal.WithContext(s.p.Context), decimal.WithContext(s.p.Context), decimal.WithContext(s.p.Context)
	var midValue, targetDiff int32
	for midValue = 0; midValue < umodulo; midValue++ {
		if midValue%1000000 == 0 {
			fmt.Printf("MidValue: %d\n", midValue)
		}
		combL := allSumCombinations(s.p.Context, midValue, umodulo, moduloMaps[:numSegments/2], tables[:numSegments/2],
			betas[:numSegments/2])
		// No pair adds up to this midValue
		if len(combL) == 0 {
			continue
		}
		targetDiff = tools.SubMod(modTarget, midValue, umodulo)
		combR := allSumCombinations(s.p.Context, targetDiff, umodulo, moduloMaps[numSegments/2:], tables[numSegments/2:],
			betas[numSegments/2:])
		// No pair adds up to targetDiff
		if len(combR) == 0 {
			continue
		}
		for sumL, solsL := range combL {
			sl.SetUint64(0)
			s.p.Context.SetString(sl, sumL)
			for sumR, solsR := range combR {
				sr.SetUint64(0)
				s.p.Context.SetString(sr, sumR)
				s.p.Context.Add(slr, sl, sr)
				if slr.Cmp(scaledTarget) != 0 {
					continue
				}
				for j := range solsL {
					for i := range solsR {
						//joint := make([]uint16, len(solsL[j])+len(solsR[i]))
						//copy(joint, solsL[j])
						//copy(joint[len(solsL[j]):], solsR[i])
						joint := append(solsL[j], solsR[i]...)
						//fmt.Println(combL)
						//fmt.Println(combR)
						//fmt.Println("-----", ArrayToString(lociToGenotype(joint, pgs.NumHaplotypes*len(s.p.Weights))),
						subsets = append(subsets, joint)
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
func calculateSubsetSumsTable(ctx decimal.Context, betas map[uint16]*decimal.Big, upperBoundOG, lowerBoundOG *decimal.Big) map[string][]uint16 {
	// Fill out the table using dynamic programming
	table := make(map[string][]uint16)
	// add the zero weight
	zero := decimal.WithContext(ctx)
	table[zero.String()] = make([]uint16, 0)
	indices := make([]uint16, 0, len(betas))
	for i := range betas {
		indices = append(indices, i)
	}
	sort.Slice(indices, func(i, j int) bool {
		return indices[i] < indices[j]
	})
	lowerBound, upperBound := decimal.WithContext(ctx), decimal.WithContext(ctx)
	lowerBound.Copy(lowerBoundOG)
	upperBound.Copy(upperBoundOG)
	prevSum, nextSum, weight := decimal.WithContext(ctx), decimal.WithContext(ctx), decimal.WithContext(ctx)
	existingSums := make([]*decimal.Big, 1)
	existingSums[0] = zero
	i := 0
	var k uint16
	var nextss string
	for _, pos := range indices {
		i++
		fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[pos].Sign() > 0 {
			ctx.Add(lowerBound, lowerBound, betas[pos])
			ctx.Add(lowerBound, lowerBound, betas[pos])
		} else {
			ctx.Add(upperBound, upperBound, betas[pos])
			ctx.Add(upperBound, upperBound, betas[pos])
		}
		newSums := make([]*decimal.Big, 0)
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.NumHaplotypes; k++ {
				ctx.Mul(weight, betas[pos], decimal.New(int64(k), 0))
				ctx.Add(nextSum, prevSum, weight)
				if nextSum.Cmp(lowerBound) >= 0 && nextSum.Cmp(upperBound) <= 0 {
					nextss = nextSum.String()
					if _, ok := table[nextss]; !ok {
						table[nextss] = make([]uint16, 0)
						newSums = append(newSums, decimal.WithContext(ctx).Copy(nextSum))
					}
					table[nextss] = append(table[nextss], pgs.NumHaplotypes*pos+k-1)
				}
			}
		}
		existingSums = append(existingSums, newSums...)
	}
	return table
}

func backtrack(ctx decimal.Context, path []uint16, sum *decimal.Big, table map[string][]uint16, weights map[uint16]*decimal.Big) [][]uint16 {
	if sum.Sign() == 0 {
		return [][]uint16{path}
	}
	output := make([][]uint16, 0)
	for _, ptr := range table[sum.String()] {
		if locusAlreadyExists(ptr, path) || (len(path) > 0 && ptr > path[len(path)-1]) {
			continue
		}
		newState := make([]uint16, len(path)+1)
		copy(newState, path)
		newState[len(path)] = ptr
		newSum := decimal.WithContext(ctx)
		newSum.Copy(weights[ptr/2])
		if ptr%2 == 1 {
			newSum.Add(newSum, weights[ptr/2])
		}
		ctx.Sub(newSum, sum, newSum)
		if res := backtrack(ctx, newState, newSum, table, weights); res != nil {
			output = append(output, res...)
		}
	}
	return output
}

func buildModuloMap(ctx decimal.Context, modulo *decimal.Big, table map[string][]uint16) map[int32][]string {
	moduloMap := make(map[int32][]string)
	sumDec := decimal.WithContext(ctx)
	var reduced int32
	var preReduced int64
	var ok bool
	for sumStr := range table {
		sumDec.SetUint64(0)
		ctx.SetString(sumDec, sumStr)
		preReduced, ok = tools.BigMod(ctx, sumDec, modulo).Int64()
		if !ok {
			log.Fatalf("Failed to reduce %s to int64", sumStr)
		}
		reduced = int32(preReduced)
		if _, ok := moduloMap[reduced]; !ok {
			moduloMap[reduced] = make([]string, 0)
		}
		moduloMap[reduced] = append(moduloMap[reduced], sumStr)
	}
	return moduloMap
}

func allSumCombinations(ctx decimal.Context, modSum, umodulo int32, modTables []map[int32][]string, fullTables []map[string][]uint16,
	betas []map[uint16]*decimal.Big) map[string][][]uint16 {
	combinations := make(map[string][][]uint16)
	var valueEntry, valueExit int32
	tmpSum, leftSum, rightSum := decimal.WithContext(ctx), decimal.WithContext(ctx), decimal.WithContext(ctx)
	var ok bool
	var s string
	for valueEntry = range modTables[0] {
		valueExit = tools.SubMod(modSum, valueEntry, umodulo)
		if _, ok = modTables[1][valueExit]; !ok {
			continue
		}
		leftSolComb := backtrackedSolutions(ctx, modTables[0][valueEntry], fullTables[0], betas[0])
		rightSolComb := backtrackedSolutions(ctx, modTables[1][valueExit], fullTables[1], betas[1])
		for nonModSumL, lSols := range leftSolComb {
			leftSum.SetUint64(0)
			ctx.SetString(leftSum, nonModSumL)
			for nonModSumR, rSols := range rightSolComb {
				rightSum.SetUint64(0)
				ctx.SetString(rightSum, nonModSumR)
				ctx.Add(tmpSum, leftSum, rightSum)
				for i := range lSols {
					for j := range rSols {
						s = tmpSum.String()
						if _, ok = combinations[s]; !ok {
							combinations[s] = make([][]uint16, 0)
						}
						combinations[s] = append(combinations[s], append(lSols[i], rSols[j]...))
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

func backtrackedSolutions(ctx decimal.Context, sums []string, table map[string][]uint16, betas map[uint16]*decimal.Big) map[string][][]uint16 {
	solutions := make(map[string][][]uint16)
	for _, sum := range sums {
		sumBig, ok := decimal.WithContext(ctx).SetString(sum)
		if !ok {
			log.Fatalf("Failed to convert %s to decimal", sum)
		}
		input := make([]uint16, 0)
		solutions[sum] = backtrack(ctx, input, sumBig, table, betas)
	}
	return solutions
}

func scaleWeights(ctx decimal.Context, weights []*decimal.Big, multiplier *decimal.Big) []*decimal.Big {
	scaledWeights := make([]*decimal.Big, len(weights))
	for i := range weights {
		scaledWeights[i] = decimal.WithContext(ctx)
		scaledWeights[i].Copy(weights[i])
		scaledWeights[i].Mul(scaledWeights[i], multiplier)
	}
	return scaledWeights
}

func makeBetaMap(betas []*decimal.Big, start, end int) map[uint16]*decimal.Big {
	bmap := make(map[uint16]*decimal.Big)
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

func getMaxTotal(ctx decimal.Context, values []*decimal.Big) (*decimal.Big, *decimal.Big) {
	positive, negative := decimal.WithContext(ctx), decimal.WithContext(ctx)
	for _, v := range values {
		if v.Sign() > 0 {
			ctx.Add(positive, positive, v)
			ctx.Add(positive, positive, v)
		} else {
			ctx.Add(negative, negative, v)
			ctx.Add(negative, negative, v)
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
