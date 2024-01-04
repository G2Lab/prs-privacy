package solver

import (
	"fmt"
	"log"
	"sort"

	"github.com/ericlagergren/decimal"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
)

type OneSplitDP struct {
	target  *decimal.Big
	p       *pgs.PGS
	rounder *Rounder
}

func NewOneSplitDP(target *decimal.Big, p *pgs.PGS) *OneSplitDP {
	s := &OneSplitDP{
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

type Rounder struct {
	RoundedMode   bool
	Ctx           decimal.Context
	ScaledWeights []*decimal.Big
	ScaledTarget  *decimal.Big
}

func (dp *OneSplitDP) Solve() map[string][]uint8 {
	var roundingError int64 = 0
	multiplier := decimal.WithContext(dp.p.Context).SetMantScale(1, -dp.target.Scale())
	if dp.target.Scale() > params.PrecisionsLimit {
		multiplier.SetMantScale(1, -params.PrecisionsLimit)
		//roundingError = int64(dp.p.VariantCount) * 5 / 4
		roundingError = int64(dp.p.VariantCount)
		dp.rounder.RoundedMode = true
		dp.rounder.ScaledWeights = scaleWeights(dp.p.Context, dp.p.Weights, multiplier)
		dp.rounder.ScaledTarget.Mul(dp.target, multiplier)
	}
	weights := bigsToInts(dp.p.Context, dp.p.Weights, multiplier)
	target, ok := decimal.WithContext(dp.p.Context).Mul(dp.target, multiplier).RoundToInt().Int64()
	if !ok {
		log.Fatalf("Failed to convert target big to int64: %s", dp.target.String())
	}
	fmt.Printf("Target: %d\n", target)
	fmt.Printf("Weights: %v\n", weights)

	numSegments := 2
	splitIdxs := []int{0, len(weights) / 2, len(weights)}
	betas := make([]map[uint16]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeIntBetaMap(weights, splitIdxs[i], splitIdxs[i+1])
	}

	tables := make([]map[int64][]uint16, numSegments)
	maxTotalPositive, maxTotalNegative := getIntMaxTotal(weights)
	//if dp.rounder.RoundedMode && maxTotalNegative < 0 {
	//	roundingErrorRight = roundingError
	//}
	//upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive-roundingErrorRight
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive

	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumIntTable(betas[i], upper, lower)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		//fmt.Println(tables[i])
		//fmt.Println(betas[i])
	}

	targets := []int64{target}
	if dp.rounder.RoundedMode {
		//for w := target + roundingErrorRight; w > target-roundingError; w-- {
		//	if w == target {
		//		continue
		//	}
		for w := target - 1; w > target-roundingError; w-- {
			targets = append(targets, w)
		}
	}

	// Do recursion to explore all the combinations
	solutionHeap := &genheap{}
	halfSums := make([][]int64, numSegments)
	step := len(tables[0]) / 10
	if step == 0 {
		step = 1
	}
	var i, s int
	var lkl float64
	var leftSum int64
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
		halfSols := make([][]*genotype, numSegments)
		for i = 0; i < numSegments; i++ {
			halfSols[i] = make([]*genotype, 0)
			for _, halfSum := range halfSums[i] {
				combinations := backtrackFromIntSum(halfSum, tables[i], betas[i])
				for _, seq := range combinations {
					lkl = calculateNegativeLikelihood(seq, splitIdxs[i]*pgs.NumHaplotypes, splitIdxs[i+1]*pgs.NumHaplotypes, dp.p)
					halfSols[i] = append(halfSols[i], newGenotype(seq, lkl))
				}
			}
		}
		// combine partial solutions
		combinePartials(0, numSegments, make([]uint16, 0), 0, decimal.WithContext(dp.p.Context), halfSols, solutionHeap, dp.rounder)
	}

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.NumHaplotypes)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumIntTable(betas map[uint16]int64, upperBound, lowerBound int64) map[int64][]uint16 {
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

func backtrackFromIntSum(sum int64, table map[int64][]uint16, betas map[uint16]int64) [][]uint16 {
	input := make([]uint16, 0)
	return backtrackInt(input, sum, table, betas)
}

func backtrackInt(path []uint16, sum int64, table map[int64][]uint16, weights map[uint16]int64) [][]uint16 {
	if sum == 0 {
		return [][]uint16{path}
	}
	output := make([][]uint16, 0)
	for _, ptr := range table[sum] {
		if locusAlreadyExists(ptr, path) || (len(path) > 0 && ptr > path[len(path)-1]) {
			continue
		}
		newState := make([]uint16, len(path)+1)
		copy(newState, path)
		newState[len(path)] = ptr
		newSum := sum - weights[ptr/2]*int64(ptr%2+1)
		if res := backtrackInt(newState, newSum, table, weights); res != nil {
			output = append(output, res...)
		}
	}
	return output
}

func combinePartials(segmentNum, totalSegments int, input []uint16, lkl float64, score *decimal.Big, genotypes [][]*genotype, solHeap *genheap, rnd *Rounder) {
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
		carryover := make([]uint16, len(input)+len(sol.mutations))
		copy(carryover, input)
		copy(carryover[len(input):], sol.mutations)
		newScore := decimal.WithContext(rnd.Ctx).Copy(score)
		if rnd.RoundedMode {
			rnd.Ctx.Add(newScore, newScore, lociToScore(rnd.Ctx, sol.mutations, rnd.ScaledWeights))
		}
		combinePartials(segmentNum+1, totalSegments, carryover, lkl+sol.likelihood, newScore, genotypes, solHeap, rnd)
	}
}

func bigsToInts(ctx decimal.Context, bigs []*decimal.Big, multiplier *decimal.Big) []int64 {
	ints := make([]int64, len(bigs))
	var ok bool
	for i, b := range bigs {
		tmp := decimal.WithContext(ctx)
		ctx.Mul(tmp, b, multiplier)
		ints[i], ok = tmp.RoundToInt().Int64()
		if !ok {
			log.Fatalf("Failed to convert big to int64: %s", b.String())
		}
	}
	return ints
}

func makeIntBetaMap(betas []int64, start, end int) map[uint16]int64 {
	bmap := make(map[uint16]int64)
	for i := start; i < end; i++ {
		bmap[uint16(i)] = betas[i]
	}
	return bmap
}

func getIntMaxTotal(values []int64) (int64, int64) {
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

func lociToScore(ctx decimal.Context, loci []uint16, weights []*decimal.Big) *decimal.Big {
	score := decimal.WithContext(ctx)
	for _, locus := range loci {
		ctx.Add(score, score, weights[locus/2])
		if locus%pgs.NumHaplotypes == 1 {
			ctx.Add(score, score, weights[locus/2])
		}
	}
	return score
}
