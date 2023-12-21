package solver

import (
	"fmt"

	"github.com/ericlagergren/decimal"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
)

type OneSplitDP struct {
	target *decimal.Big
	p      *pgs.PGS
}

func NewOneSplitDP(target *decimal.Big, p *pgs.PGS) *OneSplitDP {
	s := &OneSplitDP{
		target: target,
		p:      p,
	}
	return s
}

func (s *OneSplitDP) Solve() map[string][]uint8 {
	ctx := s.p.Context
	scale := s.target.Scale()
	multiplier := decimal.WithContext(ctx).SetMantScale(1, -scale)
	scaledWeights := scaleWeights(ctx, s.p.Weights, multiplier)
	scaledTarget := decimal.WithContext(ctx).Copy(s.target)
	scaledTarget.Mul(scaledTarget, multiplier)
	fmt.Println("Scale:", scale)
	fmt.Println("Multiplier:", multiplier.String())
	fmt.Println("Scaled target:", scaledTarget.String())
	fmt.Println("Scaled weights:", scaledWeights)

	numSegments := 2
	splitIdxs := []int{0, len(scaledWeights) / 2, len(scaledWeights)}
	betas := make([]map[uint16]*decimal.Big, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(scaledWeights, splitIdxs[i], splitIdxs[i+1])
	}

	tables := make([]map[string][]uint16, numSegments)
	maxTotalPositive, maxTotalNegative := getMaxTotal(ctx, scaledWeights)
	upper, lower := decimal.WithContext(ctx), decimal.WithContext(s.p.Context)
	s.p.Context.Sub(upper, scaledTarget, maxTotalNegative)
	s.p.Context.Sub(lower, scaledTarget, maxTotalPositive)
	fmt.Println("Upper:", upper.String())
	fmt.Println("Lower:", lower.String())
	for i := 0; i < numSegments; i++ {
		tables[i] = calculateSubsetSumsTable(s.p.Context, betas[i], upper, lower)
		fmt.Printf("Table %d len: %d\n", i, len(tables[i]))
		//fmt.Println(tables[i])
		//fmt.Println(betas[i])
	}

	// Do recursion to explore all the combinations
	partHeaps := make([]*genheap, numSegments)
	for i := 0; i < numSegments; i++ {
		partHeaps[i] = &genheap{}
	}
	solutionHeap := &genheap{}
	halfSums := make([]*decimal.Big, numSegments)
	backtracked := make([][][]uint16, numSegments)
	var i int
	var lkl float64
	for sumLeftStr := range tables[0] {
		halfSums[0], halfSums[1] = decimal.WithContext(ctx), decimal.WithContext(ctx)
		ctx.SetString(halfSums[0], sumLeftStr)
		ctx.Sub(halfSums[1], scaledTarget, halfSums[0])
		if _, ok := tables[1][halfSums[1].String()]; !ok {
			continue
		}
		halfSols := make([][]*genotype, numSegments)
		for i = 0; i < numSegments; i++ {
			backtracked[i] = backtrackFromSum(ctx, halfSums[i], tables[i], betas[i])
			halfSols[i] = make([]*genotype, 0)
			for _, seq := range backtracked[i] {
				lkl = calculateNegativeLikelihood(seq, splitIdxs[i]*pgs.NumHaplotypes, splitIdxs[i+1]*pgs.NumHaplotypes, s.p)
				if addToHeap(partHeaps[i], lkl, seq, params.SmallHeap) {
					halfSols[i] = append(halfSols[i], &genotype{seq, lkl})
				}
			}
		}
		for _, leftSol := range halfSols[0] {
			for _, rightSol := range halfSols[1] {
				lkl = leftSol.likelihood + rightSol.likelihood
				addToHeap(solutionHeap, lkl, append(leftSol.mutations, rightSol.mutations...), params.BigHeap)
			}
		}
	}

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(s.p.Weights)*pgs.NumHaplotypes)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}
