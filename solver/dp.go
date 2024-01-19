package solver

import (
	"fmt"
	"log"
	"math"
	"sort"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
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

type Node struct {
	pointers      map[uint16]uint16
	topLikelihood float64
	topPointer    uint16
}

func newNode(lkl float64, ptr uint16) *Node {
	return &Node{make(map[uint16]uint16, 1), lkl, ptr}
}

type Update struct {
	sum        int64
	likelihood float64
	forwardPtr uint16
	backPtr    uint16
}

func newUpdate(sum int64, likelihood float64, fptr, bptr uint16) Update {
	return Update{
		sum:        sum,
		likelihood: likelihood,
		forwardPtr: fptr,
		backPtr:    bptr,
	}
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

func (dp *DP) Solve() map[string][]uint8 {
	numSegments := 2
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
		//fmt.Printf("Scaled target: %s\n", dp.rounder.ScaledTarget.String())
		//fmt.Printf("Scaled weights: %v\n", dp.rounder.ScaledWeights)
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
	fmt.Printf("Weights [%d]: %v\n", len(weights), weights)

	var splitIdxs []int
	splitIdxs = []int{0, len(weights) / 2, len(weights)}
	betas := make([]map[uint16]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(weights, splitIdxs[i], splitIdxs[i+1])
	}

	maxTotalPositive, maxTotalNegative := GetMaxTotal(weights)
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive

	targets := []int64{target}
	if dp.rounder.RoundedMode {
		for w := target - 1; w > target-roundingError; w-- {
			targets = append(targets, w)
		}
	}

	solutionHeap := newGenHeap()
	if len(weights) > 50 {
		tables := make([]map[int64]*Node, numSegments)
		for i := 0; i < numSegments; i++ {
			tables[i] = calculateSubsetSumTableWithLikelihood(betas[i], upper, lower, dp.p)
			fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		}
		dp.probabilisticMitM(numSegments, tables, betas, targets, solutionHeap)
	} else {
		tables := make([]map[int64][]uint16, numSegments)
		for i := 0; i < numSegments; i++ {
			tables[i] = calculateSubsetSumTable(betas[i], upper, lower)
			fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		}
		// Do recursion to explore all the combinations
		dp.deterministicMitM(numSegments, splitIdxs, tables, betas, targets, solutionHeap)
	}

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.NumHplt)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func (dp *DP) probabilisticMitM(numSegments int, tables []map[int64]*Node, betas []map[uint16]int64, targets []int64, solHeap *genHeap) {
	matchHeapSize := 1000 * dp.p.VariantCount
	step := len(tables[0]) / 10
	if step == 0 {
		step = 1
	}
	mheap := newMatchHeap()
	var s int
	var ok bool
	var leftSum int64
	for leftSum = range tables[0] {
		if s++; s%step == 0 {
			fmt.Printf("Progress: %d%%\n", s*10/step)
		}
		for _, t := range targets {
			if _, ok = tables[1][t-leftSum]; ok {
				//fmt.Println("==", leftSum, t-leftSum, tables[0][leftSum].topLikelihood+tables[1][t-leftSum].topLikelihood)
				mheap.addToMatchHeap(tables[0][leftSum].topLikelihood+tables[1][t-leftSum].topLikelihood, []int64{leftSum, t - leftSum}, matchHeapSize)
			}
		}
	}
	score := apd.NewBigInt(0)
	for _, mtch := range *mheap {
		var solution []uint16
		for i := 0; i < numSegments; i++ {
			solution = append(solution, backtrackWithPointers(mtch.sums[i], tables[i], betas[i])...)
		}
		if dp.rounder.RoundedMode {
			score.Set(lociToScore(solution, dp.rounder.ScaledWeights))
			if score.Cmp(dp.rounder.ScaledTarget) != 0 {
				continue
			}
		}
		solHeap.addToGenHeap(mtch.likelihood, solution, params.HeapSize)
	}
}

func (dp *DP) deterministicMitM(numSegments int, splitIdxs []int, tables []map[int64][]uint16, betas []map[uint16]int64,
	targets []int64, solHeap *genHeap) {
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
		halfSums[1] = make([]int64, 0)
		for _, t := range targets {
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
					lkl = calculateNegativeLikelihood(combinations[j], splitIdxs[i]*pgs.NumHplt, splitIdxs[i+1]*pgs.NumHplt, dp.p)
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
		combinePartials(0, numSegments, make([]uint16, 0), 0, apd.NewBigInt(0), halfSols, solHeap, params.HeapSize, dp.rounder)
	}
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumTableWithLikelihood(betas map[uint16]int64, upperBound, lowerBound int64, p *pgs.PGS) map[int64]*Node {
	// Fill out the table using dynamic programming
	table := make(map[int64]*Node)
	indices := make([]uint16, 0, len(betas))
	for i := range betas {
		indices = append(indices, i)
	}
	sort.Slice(indices, func(i, j int) bool {
		return indices[i] < indices[j]
	})
	// the likelihood of all the snps being 0
	allZeroLikelihood := calculateNegativeLikelihood([]uint16{}, int(indices[0])*pgs.NumHplt,
		int(indices[len(indices)-1])*pgs.NumHplt, p)
	// add the zero weight
	table[0] = newNode(allZeroLikelihood, math.MaxUint16)
	existingSums := make([]int64, 1)
	existingSums[0] = 0

	var updates []Update
	var ok bool
	var k, nextPtr uint16
	var nextLikelihood float64
	var prevSum, nextSum, weight int64
	for i, pos := range indices {
		fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[pos] > 0 {
			lowerBound += pgs.NumHplt * betas[pos]
		} else {
			upperBound += pgs.NumHplt * betas[pos]
		}
		updates = make([]Update, 0)
		newSums := make([]int64, 0)
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.NumHplt; k++ {
				nextPtr = pgs.NumHplt*pos + k - 1
				weight = betas[pos] * int64(k)
				nextSum = prevSum + weight
				if nextSum < lowerBound || nextSum > upperBound {
					continue
				}
				nextLikelihood = table[prevSum].topLikelihood + float64(k)*mafToLikelihood(p.Maf[pos][1]) -
					float64(k)*mafToLikelihood(p.Maf[pos][0])
				if _, ok = table[nextSum]; !ok {
					table[nextSum] = newNode(nextLikelihood, nextPtr)
					newSums = append(newSums, nextSum)
				}
				table[nextSum].pointers[nextPtr] = table[prevSum].topPointer
				if nextLikelihood < table[nextSum].topLikelihood {
					updates = append(updates, newUpdate(nextSum, nextLikelihood, nextPtr, table[prevSum].topPointer))
				}
			}
		}
		// We postpone the updates to avoid the sums for the same locus stacking up on each other
		for u := range updates {
			table[updates[u].sum].topLikelihood = updates[u].likelihood
			table[updates[u].sum].topPointer = updates[u].forwardPtr
		}
		existingSums = append(existingSums, newSums...)

	}
	return table
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
			lowerBound += pgs.NumHplt * betas[pos]
		} else {
			upperBound += pgs.NumHplt * betas[pos]
		}
		newSums := make([]int64, 0)
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.NumHplt; k++ {
				weight = betas[pos] * int64(k)
				nextSum = prevSum + weight
				if nextSum >= lowerBound && nextSum <= upperBound {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
						newSums = append(newSums, nextSum)
					}
					table[nextSum] = append(table[nextSum], pgs.NumHplt*pos+k-1)
				}
			}
		}
		existingSums = append(existingSums, newSums...)
	}
	return table
}

func backtrackWithPointers(sum int64, table map[int64]*Node, weights map[uint16]int64) []uint16 {
	output := make([]uint16, 0)
	nextPtr := table[sum].topPointer
	newSum := sum
	for sum != 0 && nextPtr != math.MaxUint16 {
		output = append(output, nextPtr)
		newSum -= weights[nextPtr/2] * int64(nextPtr%2+1)
		nextPtr = table[sum].pointers[nextPtr]
		sum = newSum
	}
	return output
}

func combinePartials(segmentNum, totalSegments int, input []uint16, lkl float64, score *apd.BigInt, genotypes [][]*genotype, solHeap *genHeap, heapSize int, rnd *Rounder) {
	if segmentNum == totalSegments {
		if rnd.RoundedMode {
			if score.Cmp(rnd.ScaledTarget) != 0 {
				return
			}
		}
		solHeap.addToGenHeap(lkl, input, heapSize)
		return
	}
	for _, sol := range genotypes[segmentNum] {
		input = append(input, sol.mutations...)
		if rnd.RoundedMode {
			newScore := apd.NewBigInt(0)
			newScore.Add(score, lociToScore(sol.mutations, rnd.ScaledWeights))
			combinePartials(segmentNum+1, totalSegments, input, lkl+sol.likelihood, newScore, genotypes, solHeap, heapSize, rnd)
		} else {
			combinePartials(segmentNum+1, totalSegments, input, lkl+sol.likelihood, score, genotypes, solHeap, heapSize, rnd)
		}
		input = input[:len(input)-len(sol.mutations)]
	}
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
		if locus%pgs.NumHplt == 1 {
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
	for i := 0; i < len(snps); i += pgs.NumHplt {
		for j := 0; j < pgs.NumHplt; j++ {
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
