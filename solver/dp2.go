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
	target *decimal.Big
	p      *pgs.PGS
}

func NewTwoSplitDP(target *decimal.Big, p *pgs.PGS) *TwoSplitDP {
	s := &TwoSplitDP{
		target: target,
		p:      p,
	}
	return s
}

type Product struct {
	Factors []*decimal.Big
	Sum     *decimal.Big
}

func (s *TwoSplitDP) Solve() map[string][]uint8 {
	dctx := s.p.Context
	scale := s.target.Scale()
	multiplier := decimal.WithContext(dctx).SetMantScale(1, -scale)
	fmt.Println("Scale:", scale)
	fmt.Println("Multiplier:", multiplier.String())
	scaledWeights := scaleWeights(dctx, s.p.Weights, multiplier)
	scaledTarget := decimal.WithContext(dctx)
	scaledTarget.Copy(s.target)
	scaledTarget.Mul(scaledTarget, multiplier)
	fmt.Println("Scaled target:", scaledTarget.String())
	fmt.Println("Scaled weights:", scaledWeights)

	numSegments := 4
	splitIdxHalf := len(scaledWeights) / 2
	splitIdxs := []int{0, splitIdxHalf / 2, splitIdxHalf, splitIdxHalf + (len(scaledWeights)-splitIdxHalf)/2, len(scaledWeights)}
	betas := make([]map[uint16]*decimal.Big, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(scaledWeights, splitIdxs[i], splitIdxs[i+1])
	}

	tables := make([]map[string][]uint16, numSegments)
	maxTotalPositive, maxTotalNegative := getMaxTotal(dctx, scaledWeights)
	upper, lower := decimal.WithContext(s.p.Context), decimal.WithContext(dctx)
	dctx.Sub(upper, scaledTarget, maxTotalNegative)
	dctx.Sub(lower, scaledTarget, maxTotalPositive)
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumsTable(dctx, betas[i], upper, lower)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		//fmt.Println(tables[i])
		//fmt.Println(betas[i])
	}

	umodulo := int32(tools.FindNextSmallerPrime(uint64(len(tables[0]))))
	modulo := decimal.WithContext(dctx)
	modulo.SetUint64(uint64(umodulo))
	fmt.Printf("Modulo: %s\n", modulo.String())
	moduloMaps := make([]map[int32][]*decimal.Big, numSegments)
	for i := 0; i < numSegments; i++ {
		moduloMaps[i] = buildModuloMap(dctx, modulo, tables[i])
		//fmt.Println(moduloMaps[i])
	}

	// Combine and backtrack
	tmp, ok := tools.BigMod(dctx, scaledTarget, modulo).Int64()
	if !ok {
		log.Fatalf("Failed to convert %s to int64", scaledTarget.String())
	}
	modTarget := int32(tmp)

	// Do recursion to explore all the combinations
	solutionHeap := &genheap{}
	var combine func(int, []uint16, float64, [][]*genotype)
	combine = func(k int, input []uint16, lkl float64, genotypes [][]*genotype) {
		if k == numSegments {
			addToHeap(solutionHeap, lkl, input, params.BigHeap)
			return
		}
		for _, sol := range genotypes[k] {
			carryover := make([]uint16, len(input)+len(sol.mutations))
			copy(carryover, input)
			copy(carryover[len(input):], sol.mutations)
			combine(k+1, carryover, lkl+sol.likelihood, genotypes)
		}
	}

	//fmt.Printf("Mod target: %d\n", modTarget)
	sumLR := decimal.WithContext(s.p.Context)
	backtracked := make([]map[string][]*genotype, numSegments)
	for i := 0; i < numSegments; i++ {
		backtracked[i] = make(map[string][]*genotype)
	}
	products := make([][]*Product, numSegments/2)
	for i := 0; i < len(products); i++ {
		products[i] = make([]*Product, 0)
	}
	step := umodulo / 10
	if step == 0 {
		step = 1
	}
	var midValue, targetDiff int32
	var k int
	for midValue = 0; midValue < umodulo; midValue++ {
		if midValue%step == 0 {
			fmt.Printf("MidValue: %d\n", midValue)
		}
		products[0] = getMatchingSums(s.p.Context, midValue, umodulo, moduloMaps[:numSegments/2])
		// No pair adds up to this midValue
		if len(products[0]) == 0 {
			continue
		}
		targetDiff = tools.SubMod(modTarget, midValue, umodulo)
		products[1] = getMatchingSums(s.p.Context, targetDiff, umodulo, moduloMaps[numSegments/2:])
		// No pair adds up to targetDiff
		if len(products[1]) == 0 {
			continue
		}
		matches := make([][]*decimal.Big, 0)
		for _, leftPair := range products[0] {
			for _, rigthPair := range products[1] {
				sumLR.SetUint64(0)
				s.p.Context.Add(sumLR, leftPair.Sum, rigthPair.Sum)
				if sumLR.Cmp(scaledTarget) != 0 {
					continue
				}
				// Get all the parts of the valid solution
				matches = append(matches, []*decimal.Big{leftPair.Factors[0], leftPair.Factors[1], rigthPair.Factors[0], rigthPair.Factors[1]})
			}
		}
		genotypes := make([][]*genotype, numSegments)
		for _, match := range matches {
			for k = 0; k < numSegments; k++ {
				if genotypes[k], ok = backtracked[k][match[k].String()]; !ok {
					combinations := backtrackFromSum(dctx, match[k], tables[k], betas[k])
					backtracked[k][match[k].String()] = make([]*genotype, len(combinations))
					for l := range combinations {
						backtracked[k][match[k].String()][l] = newGenotype(combinations[l],
							calculateNegativeLikelihood(combinations[l], splitIdxs[k]*pgs.NumHaplotypes, splitIdxs[k+1]*pgs.NumHaplotypes, s.p))
					}
					genotypes[k] = backtracked[k][match[k].String()]
				}
			}
			combine(0, []uint16{}, 0, genotypes)
		}
	}
	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(s.p.Weights)*pgs.NumHaplotypes)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
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
		//fmt.Printf("Position %d/%d\n", i, len(betas))
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

func buildModuloMap(ctx decimal.Context, modulo *decimal.Big, table map[string][]uint16) map[int32][]*decimal.Big {
	moduloMap := make(map[int32][]*decimal.Big)
	var reduced int32
	var preReduced int64
	var ok bool
	for sumStr := range table {
		sumDec := decimal.WithContext(ctx)
		ctx.SetString(sumDec, sumStr)
		preReduced, ok = tools.BigMod(ctx, sumDec, modulo).Int64()
		if !ok {
			log.Fatalf("Failed to reduce %s to int64", sumStr)
		}
		reduced = int32(preReduced)
		if _, ok := moduloMap[reduced]; !ok {
			moduloMap[reduced] = make([]*decimal.Big, 0)
		}
		moduloMap[reduced] = append(moduloMap[reduced], sumDec)
	}
	return moduloMap
}

func getMatchingSums(ctx decimal.Context, modSum, umodulo int32, modTables []map[int32][]*decimal.Big) []*Product {
	var valueEntry, valueExit int32
	var ok bool
	matches := make([]*Product, 0)
	for valueEntry = range modTables[0] {
		valueExit = tools.SubMod(modSum, valueEntry, umodulo)
		if _, ok = modTables[1][valueExit]; !ok {
			continue
		}
		for _, leftSum := range modTables[0][valueEntry] {
			for _, rightSum := range modTables[1][valueExit] {
				jointSum := decimal.WithContext(ctx)
				ctx.Add(jointSum, leftSum, rightSum)
				matches = append(matches, &Product{[]*decimal.Big{leftSum, rightSum}, jointSum})
			}
		}
	}
	return matches
}

func backtrackFromSum(ctx decimal.Context, sum *decimal.Big, table map[string][]uint16, betas map[uint16]*decimal.Big) [][]uint16 {
	input := make([]uint16, 0)
	return backtrack(ctx, input, sum, table, betas)
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

func makeBetaMap(betas []*decimal.Big, start, end int) map[uint16]*decimal.Big {
	bmap := make(map[uint16]*decimal.Big)
	for i := start; i < end; i++ {
		bmap[uint16(i)] = betas[i]
	}
	return bmap
}

func addToHeap(h *genheap, lkl float64, sol []uint16, heapSize int) bool {
	switch {
	case h.Len() < heapSize:
		heap.Push(h, genotype{sol, lkl})
		return true
	case lkl < (*h)[0].likelihood:
		heap.Pop(h)
		heap.Push(h, genotype{sol, lkl})
		return true
	default:
		return false
	}
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
