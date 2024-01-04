package solver

import (
	"container/heap"
	"fmt"
	"log"
	"sort"

	"github.com/ericlagergren/decimal"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

type TwoSplitDP struct {
	target  *decimal.Big
	p       *pgs.PGS
	rounder *Rounder
}

func NewTwoSplitDP(target *decimal.Big, p *pgs.PGS) *TwoSplitDP {
	s := &TwoSplitDP{
		target: target,
		p:      p,
		rounder: &Rounder{
			RoundedMode:   false,
			Ctx:           p.Context,
			ScaledWeights: make([]*decimal.Big, 0),
			ScaledTarget:  decimal.WithContext(p.Context),
		},
	}
	return s
}

type Product struct {
	Factors []int64
	Sum     int64
}

func (dp *TwoSplitDP) Solve() map[string][]uint8 {
	var roundingErrorLeft int64 = 0
	multiplier := decimal.WithContext(dp.p.Context).SetMantScale(1, -dp.target.Scale())
	if dp.target.Scale() > params.PrecisionsLimit {
		multiplier.SetMantScale(1, -params.PrecisionsLimit)
		roundingErrorLeft = int64(dp.p.VariantCount)
		dp.rounder.RoundedMode = true
		dp.rounder.ScaledWeights = scaleWeights(dp.p.Context, dp.p.Weights, multiplier)
		dp.rounder.ScaledTarget.Mul(dp.target, multiplier)
	}
	weights := bigsToInts(dp.p.Context, dp.p.Weights, multiplier)
	target, ok := decimal.WithContext(dp.p.Context).Mul(dp.target, multiplier).RoundToInt().Int64()
	if !ok {
		log.Fatalf("Failed to convert target big to int64: %dp", dp.target.String())
	}
	fmt.Printf("Target: %d\n", target)
	fmt.Printf("Weights: %v\n", weights)

	numSegments := 4
	splitIdxHalf := len(weights) / 2
	splitIdxs := []int{0, splitIdxHalf / 2, splitIdxHalf, splitIdxHalf + (len(weights)-splitIdxHalf)/2, len(weights)}
	betas := make([]map[uint16]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(weights, splitIdxs[i], splitIdxs[i+1])
	}

	tables := make([]map[int64][]uint16, numSegments)
	maxTotalPositive, maxTotalNegative := getIntMaxTotal(weights)
	upper, lower := target-maxTotalNegative+roundingErrorLeft, target-maxTotalPositive
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumTable(betas[i], upper, lower)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
	}

	modulo := int64(tools.FindNextSmallerPrime(uint64(len(tables[0]))))
	fmt.Printf("Modulo: %d\n", modulo)
	moduloMaps := make([]map[int32][]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		moduloMaps[i] = buildModuloMap(modulo, tables[i])
	}

	targets := []int64{target}
	if dp.rounder.RoundedMode {
		for w := target - 1; w > target-roundingErrorLeft; w-- {
			targets = append(targets, w)
		}
	}
	modTargets := make([]int32, len(targets))
	for i := range targets {
		modTargets[i] = int32(tools.Mod(targets[i], modulo))
	}

	// Do recursion to explore all the combinations
	solutionHeap := &genheap{}
	backtracked := make([]map[int64][]*genotype, numSegments)
	for i := 0; i < numSegments; i++ {
		backtracked[i] = make(map[int64][]*genotype)
	}
	products := make([][]*Product, numSegments/2)
	for i := 0; i < len(products); i++ {
		products[i] = make([]*Product, 0)
	}
	step := int32(modulo) / 10
	if step == 0 {
		step = 1
	}
	var midValue, targetDiff, mt int32
	var sumLR int64
	var lkl float64
	var k, j int
	for midValue = 0; midValue < int32(modulo); midValue++ {
		matches := make([][]int64, 0)
		if midValue%step == 0 {
			fmt.Printf("MidValue: %d\n", midValue)
		}
		products[0] = getMatchingSums(midValue, int32(modulo), moduloMaps[:numSegments/2])
		// No pair adds up to this midValue
		if len(products[0]) == 0 {
			continue
		}
		for j, mt = range modTargets {
			targetDiff = tools.SubMod(mt, midValue, int32(modulo))
			products[1] = getMatchingSums(targetDiff, int32(modulo), moduloMaps[numSegments/2:])
			// No pair adds up to targetDiff
			if len(products[1]) == 0 {
				continue
			}
			for _, leftPair := range products[0] {
				for _, rightPair := range products[1] {
					sumLR = leftPair.Sum + rightPair.Sum
					if sumLR != targets[j] {
						continue
					}
					// Get all the parts of the valid solution
					matches = append(matches, []int64{leftPair.Factors[0], leftPair.Factors[1], rightPair.Factors[0], rightPair.Factors[1]})
				}
			}
		}
		partSols := make([][]*genotype, numSegments)
		for _, match := range matches {
			for k = 0; k < numSegments; k++ {
				if _, ok = backtracked[k][match[k]]; !ok {
					combinations := backtrackFromSum(match[k], tables[k], betas[k])
					backtracked[k][match[k]] = make([]*genotype, 0, len(combinations))
					for l := range combinations {
						lkl = calculateNegativeLikelihood(combinations[l], splitIdxs[k]*pgs.NumHaplotypes, splitIdxs[k+1]*pgs.NumHaplotypes, dp.p)
						backtracked[k][match[k]] = append(backtracked[k][match[k]], newGenotype(combinations[l], lkl))
					}
				}
				partSols[k] = backtracked[k][match[k]]
			}
			combinePartials(0, numSegments, make([]uint16, 0), 0, decimal.WithContext(dp.p.Context), partSols, solutionHeap, dp.rounder)
		}
	}
	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.NumHaplotypes)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func buildModuloMap(modulo int64, table map[int64][]uint16) map[int32][]int64 {
	moduloMap := make(map[int32][]int64)
	var reduced int32
	var preReduced int64
	for sum := range table {
		preReduced = tools.Mod(sum, modulo)
		reduced = int32(preReduced)
		if _, ok := moduloMap[reduced]; !ok {
			moduloMap[reduced] = make([]int64, 0)
		}
		moduloMap[reduced] = append(moduloMap[reduced], sum)
	}
	return moduloMap
}

func getMatchingSums(modSum, modulo int32, modTables []map[int32][]int64) []*Product {
	var valueEntry, valueExit int32
	var ok bool
	matches := make([]*Product, 0)
	for valueEntry = range modTables[0] {
		valueExit = tools.SubMod(modSum, valueEntry, modulo)
		if _, ok = modTables[1][valueExit]; !ok {
			continue
		}
		for _, leftSum := range modTables[0][valueEntry] {
			for _, rightSum := range modTables[1][valueExit] {
				matches = append(matches, &Product{[]int64{leftSum, rightSum}, leftSum + rightSum})
			}
		}
	}
	return matches
}

func scaleWeights(ctx decimal.Context, weights []*decimal.Big, multiplier *decimal.Big) []*decimal.Big {
	scaled := make([]*decimal.Big, len(weights))
	for i := range weights {
		scaled[i] = decimal.WithContext(ctx)
		scaled[i].Copy(weights[i])
		scaled[i].Mul(scaled[i], multiplier)
	}
	return scaled
}

func sortByLikelihood(candidates [][]uint16, totalLen int, p *pgs.PGS) [][]uint16 {
	likelihoods := make([]float64, len(candidates))
	for i, candidate := range candidates {
		likelihoods[i] = calculateNegativeLikelihood(candidate, 0, totalLen, p)
	}
	sortBy(candidates, likelihoods)
	return candidates
}

func selectTopLikelihoodCandidates(candidates [][]uint16, totalLen int, p *pgs.PGS, top int) [][]uint16 {
	if len(candidates) <= top {
		return candidates
	}
	h := &genheap{}
	// Push the first N slices onto the genheap
	for i := 0; i < top; i++ {
		heap.Push(h, genotype{candidates[i], calculateNegativeLikelihood(candidates[i], 0, totalLen, p)})
	}
	var lkl float64
	// Update the genheap with remaining slices
	for i := top; i < len(candidates); i++ {
		lkl = calculateNegativeLikelihood(candidates[i], 0, totalLen, p)
		if lkl < (*h)[0].likelihood {
			// Replace the minimum element with the current slice
			heap.Pop(h)
			heap.Push(h, genotype{candidates[i], lkl})
		}
	}

	// Extract slices from the genheap
	result := make([][]uint16, top)
	for i := top - 1; i >= 0; i-- {
		result[i] = heap.Pop(h).(genotype).mutations
	}

	return result
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

func sortedIndices(values []*decimal.Big) []int {
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
