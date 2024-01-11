package solver

import (
	"container/heap"
	"fmt"
	"log"
	"sort"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

type DP struct {
	target  *apd.Decimal
	p       *pgs.PGS
	rounder *Rounder
}

type Rounder struct {
	RoundedMode   bool
	ScaledWeights []*apd.BigInt
	ScaledTarget  *apd.BigInt
}

type Product struct {
	Factors []int64
	Sum     int64
}

func NewDP(target *apd.Decimal, p *pgs.PGS) *DP {
	s := &DP{
		target: target,
		p:      p,
		rounder: &Rounder{
			RoundedMode:   false,
			ScaledWeights: make([]*apd.BigInt, 0),
			ScaledTarget:  new(apd.BigInt),
		},
	}
	return s
}

func (dp *DP) Solve(numSegments int) map[string][]uint8 {
	var err error
	var roundingError int64 = 0
	multiplier := apd.New(1, int32(dp.p.WeightPrecision))
	if dp.p.WeightPrecision > params.PrecisionsLimit {
		//roundingError = int64(dp.p.VariantCount) * 5 / 4
		roundingError = int64(dp.p.VariantCount)
		dp.rounder.RoundedMode = true
		dp.rounder.ScaledWeights = scaleWeights(dp.p.Context, dp.p.Weights, multiplier)
		sct := new(apd.Decimal)
		dp.p.Context.Mul(sct, dp.target, multiplier)
		dp.rounder.ScaledTarget.SetString(sct.String(), 10)
		multiplier.SetFinite(1, params.PrecisionsLimit)
		fmt.Printf("Scaled target: %s\n", dp.rounder.ScaledTarget.String())
		fmt.Printf("Scaled weights: %v\n", dp.rounder.ScaledWeights)
	}
	weights := DecimalsToInts(dp.p.Context, dp.p.Weights, multiplier)
	tmp := new(apd.Decimal)
	_, err = dp.p.Context.Mul(tmp, dp.target, multiplier)
	if err != nil {
		log.Fatalf("Failed to multiply target and multiplier: %s", dp.target.String())
	}
	_, err = dp.p.Context.RoundToIntegralValue(tmp, tmp)
	if err != nil {
		log.Fatalf("Failed to round target: %s", tmp.String())
	}
	target, err := tmp.Int64()
	if err != nil {
		log.Fatalf("Failed to convert target decimal to int64: %s", tmp.String())
	}
	fmt.Printf("Target: %d\n", target)
	//fmt.Printf("Weights: %v\n", weights)

	sort.Slice(weights, func(i, j int) bool {
		return weights[i] < weights[j]
	})
	fmt.Printf("Sorted weights: %v\n", weights)

	var splitIdxs []int
	if numSegments == 2 {
		splitIdxs = []int{0, len(weights) / 2, len(weights)}
	} else if numSegments == 4 {
		splitIdxHalf := len(weights) / 2
		splitIdxs = []int{0, splitIdxHalf / 2, splitIdxHalf, splitIdxHalf + (len(weights)-splitIdxHalf)/2, len(weights)}
	}
	betas := make([]map[uint16]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(weights, splitIdxs[i], splitIdxs[i+1])
	}

	maxTotalPositive, maxTotalNegative := GetMaxTotal(weights)
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive
	segmentTargetUpperLimit, segmentTargetLowerLimit := make([]int64, numSegments), make([]int64, numSegments)
	//maxTotalPositive, maxTotalNegative := make([]int64, numSegments), make([]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		segmentTargetUpperLimit[i], segmentTargetLowerLimit[i] = SampleMaxMinScores(splitIdxs[i], splitIdxs[i+1],
			100*(splitIdxs[i+1]-splitIdxs[i]), betas[i], dp.p)
		//maxTotalPositive[i], maxTotalNegative[i] = GetMaxTotal(weights[splitIdxs[i]:splitIdxs[i+1]])
	}

	tables := make([]map[int64][]uint16, numSegments)
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumTable(betas[i], upper, lower)
		//tables[i] = calculateSubsetSumTable(betas[i], segmentTargetUpperLimit[i]-maxTotalNegative[i],
		//	segmentTargetLowerLimit[i]-maxTotalPositive[i])
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
	}

	targets := []int64{target}
	if dp.rounder.RoundedMode {
		for w := target - 1; w > target-roundingError; w-- {
			targets = append(targets, w)
		}
	}

	// Do recursion to explore all the combinations
	solutionHeap := &genheap{}
	if numSegments == 2 {
		dp.oneSplitDP(numSegments, splitIdxs, tables, betas, targets, segmentTargetUpperLimit,
			segmentTargetLowerLimit, solutionHeap)
	} else if numSegments == 4 {
		dp.twoSplitDP(numSegments, splitIdxs, tables, betas, targets, segmentTargetUpperLimit,
			segmentTargetLowerLimit, solutionHeap)
	}

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.NumHaplotypes)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func (dp *DP) oneSplitDP(numSegments int, splitIdxs []int, tables []map[int64][]uint16, betas []map[uint16]int64,
	targets []int64, upperLimits, lowerLimits []int64, solHeap *genheap) {
	halfSums := make([][]int64, numSegments)
	backtracked := make([]map[int64][]*genotype, numSegments)
	for i := 0; i < numSegments; i++ {
		backtracked[i] = make(map[int64][]*genotype)
	}
	step := len(tables[0]) / 10
	if step == 0 {
		step = 1
	}
	var i, s int
	var ok bool
	var lkl float64
	var leftSum int64
	var combinations [][]uint16
	halfSols := make([][]*genotype, numSegments)
	for leftSum = range tables[0] {
		if s++; s%step == 0 {
			fmt.Printf("Progress: %d%%\n", s*10/step)
		}
		if leftSum > upperLimits[0] || leftSum < lowerLimits[0] {
			continue
		}
		halfSums[1] = make([]int64, 0)
		for _, t := range targets {
			if t-leftSum > upperLimits[1] || t-leftSum < lowerLimits[1] {
				continue
			}
			if _, ok = tables[1][t-leftSum]; ok {
				halfSums[1] = append(halfSums[1], t-leftSum)
			}
		}
		if len(halfSums[1]) == 0 {
			continue
		}
		halfSums[0] = []int64{leftSum}
		// backtrack partial solutions
		for i = 0; i < numSegments; i++ {
			halfSols[i] = make([]*genotype, 0)
			for _, halfSum := range halfSums[i] {
				combinations = backtrackFromSum(halfSum, tables[i], betas[i])
				for j := range combinations {
					lkl = calculateNegativeLikelihood(combinations[j], splitIdxs[i]*pgs.NumHaplotypes, splitIdxs[i+1]*pgs.NumHaplotypes, dp.p)
					halfSols[i] = append(halfSols[i], newGenotype(combinations[j], lkl))
				}
			}
			//halfSols[i] = make([]*genotype, 0)
			//for _, halfSum := range halfSums[i] {
			//	if _, ok = backtracked[i][halfSum]; !ok {
			//		combinations := backtrackFromSum(halfSum, tables[i], betas[i])
			//		backtracked[i][halfSum] = make([]*genotype, 0, len(combinations))
			//		for l := range combinations {
			//			lkl = calculateNegativeLikelihood(combinations[l], splitIdxs[i]*pgs.NumHaplotypes, splitIdxs[i+1]*pgs.NumHaplotypes, dp.p)
			//			backtracked[i][halfSum] = append(backtracked[i][halfSum], newGenotype(combinations[l], lkl))
			//		}
			//	}
			//	halfSols[i] = append(halfSols[i], backtracked[i][halfSum]...)
			//}
		}
		// combine partial solutions
		combinePartials(0, numSegments, make([]uint16, 0), 0, apd.NewBigInt(0), halfSols, solHeap, dp.rounder)
	}
}

func (dp *DP) twoSplitDP(numSegments int, splitIdxs []int, tables []map[int64][]uint16, betas []map[uint16]int64,
	targets []int64, upperLimits, lowerLimits []int64, solHeap *genheap) {
	var ok bool
	modulo := int64(tools.FindNextSmallerPrime(uint64(len(tables[0]))))
	fmt.Printf("Modulo: %d\n", modulo)
	moduloMaps := make([]map[int32][]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		moduloMaps[i] = buildModuloMap(modulo, tables[i])
	}
	modTargets := make([]int32, len(targets))
	for i := range targets {
		modTargets[i] = int32(tools.Mod(targets[i], modulo))
	}

	// Do recursion to explore all the combinations
	backtracked := make([]map[int64][]*genotype, numSegments)
	for i := 0; i < numSegments; i++ {
		backtracked[i] = make(map[int64][]*genotype)
	}
	step := int32(modulo) / 100
	if step == 0 {
		step = 1
	}
	products := make([][]*Product, numSegments/2)
	var midValue, targetDiff, mt int32
	var sumLR int64
	var lkl float64
	var k, j int
	var matches [][]int64
	for midValue = 0; midValue < int32(modulo); midValue++ {
		matches = make([][]int64, 0)
		if midValue%step == 0 {
			fmt.Printf("MidValue: %d\n", midValue)
		}
		products[0] = getMatchingSums(midValue, int32(modulo), moduloMaps[:numSegments/2], upperLimits[:numSegments/2],
			lowerLimits[:numSegments/2])
		// No pair adds up to this midValue
		if len(products[0]) == 0 {
			continue
		}
		for j, mt = range modTargets {
			targetDiff = tools.SubMod(mt, midValue, int32(modulo))
			products[1] = getMatchingSums(targetDiff, int32(modulo), moduloMaps[numSegments/2:],
				upperLimits[numSegments/2:], lowerLimits[numSegments/2:])
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
			combinePartials(0, numSegments, make([]uint16, 0, len(dp.p.Weights)*pgs.NumHaplotypes), 0, apd.NewBigInt(0), partSols, solHeap, dp.rounder)
		}
	}
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumTable(betas map[uint16]int64, upperBound, lowerBound int64) map[int64][]uint16 {
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
	existingSums := make([]int64, 1)
	existingSums[0] = 0
	var k, i uint16
	var prevSum, nextSum, weight int64
	for _, pos := range indices {
		i++
		//fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[pos] > 0 {
			lowerBound += pgs.NumHaplotypes * betas[pos]
		} else {
			upperBound += pgs.NumHaplotypes * betas[pos]
		}
		newSums := make([]int64, 0)
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.NumHaplotypes; k++ {
				weight = betas[pos] * int64(k)
				nextSum = prevSum + weight
				if nextSum >= lowerBound && nextSum <= upperBound {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
						newSums = append(newSums, nextSum)
					}
					table[nextSum] = append(table[nextSum], pgs.NumHaplotypes*pos+k-1)
				}
			}
		}
		existingSums = append(existingSums, newSums...)
	}
	return table
}

func backtrackFromSum(sum int64, table map[int64][]uint16, betas map[uint16]int64) [][]uint16 {
	input := make([]uint16, 0)
	output := make([][]uint16, 0)
	backtrack(input, &output, sum, table, betas)
	return output
}

func backtrack(path []uint16, result *[][]uint16, sum int64, table map[int64][]uint16, weights map[uint16]int64) {
	if sum == 0 {
		sol := make([]uint16, len(path))
		copy(sol, path)
		*result = append(*result, sol)
	}
	for _, ptr := range table[sum] {
		if locusAlreadyExists(ptr, path) || (len(path) > 0 && ptr > path[len(path)-1]) {
			continue
		}
		path = append(path, ptr)
		newSum := sum - weights[ptr/2]*int64(ptr%2+1)
		backtrack(path, result, newSum, table, weights)
		path = path[:len(path)-1]
	}
}

func combinePartials(segmentNum, totalSegments int, input []uint16, lkl float64, score *apd.BigInt, genotypes [][]*genotype, solHeap *genheap, rnd *Rounder) {
	if segmentNum == totalSegments {
		if rnd.RoundedMode {
			if score.Cmp(rnd.ScaledTarget) != 0 {
				return
			}
		}
		addToHeap(solHeap, lkl, input, params.HeapSize)
		return
	}
	for _, sol := range genotypes[segmentNum] {
		input = append(input, sol.mutations...)
		if rnd.RoundedMode {
			newScore := apd.NewBigInt(0)
			newScore.Add(score, lociToScore(sol.mutations, rnd.ScaledWeights))
			combinePartials(segmentNum+1, totalSegments, input, lkl+sol.likelihood, newScore, genotypes, solHeap, rnd)
		} else {
			combinePartials(segmentNum+1, totalSegments, input, lkl+sol.likelihood, score, genotypes, solHeap, rnd)
		}
		input = input[:len(input)-len(sol.mutations)]
	}
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

func getMatchingSums(modSum, modulo int32, modTables []map[int32][]int64, upperLimits, lowerLimits []int64) []*Product {
	var valueEntry, valueExit int32
	var ok bool
	matches := make([]*Product, 0)
	for valueEntry = range modTables[0] {
		valueExit = tools.SubMod(modSum, valueEntry, modulo)
		if _, ok = modTables[1][valueExit]; !ok {
			continue
		}
		for _, leftSum := range modTables[0][valueEntry] {
			if leftSum > upperLimits[0] || leftSum < lowerLimits[0] {
				continue
			}
			for _, rightSum := range modTables[1][valueExit] {
				if rightSum > upperLimits[1] || rightSum < lowerLimits[1] {
					continue
				}
				matches = append(matches, &Product{[]int64{leftSum, rightSum}, leftSum + rightSum})
			}
		}
	}
	return matches
}

func DecimalsToInts(ctx *apd.Context, decimals []*apd.Decimal, multiplier *apd.Decimal) []int64 {
	ints := make([]int64, len(decimals))
	tmp := new(apd.Decimal)
	var err error
	for i, b := range decimals {
		ctx.Mul(tmp, b, multiplier)
		_, err = ctx.RoundToIntegralValue(tmp, tmp)
		if err != nil {
			log.Fatalf("Failed to round decimal: %s", b.String())
		}
		ints[i], err = tmp.Int64()
		if err != nil {
			log.Fatalf("Failed to convert decimal to int64: %s", tmp.String())
		}
	}
	return ints
}

func makeBetaMap(betas []int64, start, end int) map[uint16]int64 {
	bmap := make(map[uint16]int64, end-start)
	for i := start; i < end; i++ {
		bmap[uint16(i)] = betas[i]
	}
	return bmap
}

func scaleWeights(ctx *apd.Context, weights []*apd.Decimal, multiplier *apd.Decimal) []*apd.BigInt {
	var err error
	var ok bool
	scaled := make([]*apd.BigInt, len(weights))
	for i := range weights {
		tmp := new(apd.Decimal)
		_, err = ctx.Mul(tmp, weights[i], multiplier)
		if err != nil {
			log.Fatalf("Error scaling weights: %v", err)
		}
		_, err = ctx.RoundToIntegralValue(tmp, tmp)
		if err != nil {
			log.Fatalf("Error rounding weights: %v", err)
		}
		scaled[i] = new(apd.BigInt)
		_, ok = scaled[i].SetString(tmp.String(), 10)
		if !ok {
			log.Fatalf("Error converting scaled weight to apd.BigInt: %s", tmp.String())
		}
	}
	return scaled
}

func GetMaxTotal(values []int64) (int64, int64) {
	var positive, negative int64 = 0, 0
	for _, v := range values {
		if v > 0 {
			positive += 2 * v
		} else {
			negative += 2 * v
		}
	}
	return positive, negative
}

func SampleMaxMinScores(segmentStart, segmentEnd, numSamples int, betas map[uint16]int64, p *pgs.PGS) (int64, int64) {
	var err error
	var sample []uint8
	var score, maxScore, minScore int64
	for i := 0; i < numSamples; i++ {
		sample, err = p.SampleSegmentFromPopulation(segmentStart, segmentEnd)
		if err != nil {
			log.Fatalf("Error sampling segment: %v", err)
		}
		score = genotypeToScore(segmentStart, segmentEnd, sample, betas)
		if score > maxScore {
			maxScore = score
		}
		if score < minScore {
			minScore = score
		}
	}
	return maxScore, minScore
}

func lociToScore(loci []uint16, weights []*apd.BigInt) *apd.BigInt {
	score := apd.NewBigInt(0)
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
	for _, locus := range loci {
		sol[locus] = 1
		if locus%2 == 1 {
			sol[locus-1] = 1
		}
	}
	return sol
}

func genotypeToScore(start, end int, snps []uint8, betas map[uint16]int64) int64 {
	score := int64(0)
	for i := 0; i < len(snps); i += pgs.NumHaplotypes {
		for j := 0; j < pgs.NumHaplotypes; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				score += betas[uint16(start+i/2)]
			default:
				log.Printf("Invalid alelle value: %d", snps[i+j])
			}
		}
	}
	return score
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
