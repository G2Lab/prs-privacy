package solver

import (
	"fmt"
	"math"
	"sync"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
)

type DP struct {
	target  *apd.Decimal
	p       *pgs.PGS
	stats   *pgs.Statistics
	ppl     string
	rounder *Rounder
	known   map[int]uint8
}

type Node struct {
	Pointers      map[uint8]uint8
	TopPointer    uint8
	TopLikelihood float32
}

type ProbabilisticState struct {
	Nodes []map[int64]*Node
	Betas []map[uint8]int64
}

type State interface {
	SaveState(*State, string)
}

type DeterministicState struct {
	Indices [][]int
	Tables  []map[int64][]uint8
	Betas   []map[uint8]int64
}

type Update struct {
	sum        int64
	likelihood float32
	forwardPtr uint8
	backPtr    uint8
}

func newProbabilisticState(nodes []map[int64]*Node, betas []map[uint8]int64) *ProbabilisticState {
	return &ProbabilisticState{Nodes: nodes, Betas: betas}
}

func newNode(ptr uint8, lkl float32) *Node {
	return &Node{
		Pointers:      make(map[uint8]uint8, 1),
		TopPointer:    ptr,
		TopLikelihood: lkl,
	}
}

func newUpdate(sum int64, likelihood float32, fptr, bptr uint8) Update {
	return Update{
		sum:        sum,
		likelihood: likelihood,
		forwardPtr: fptr,
		backPtr:    bptr,
	}
}

func NewDP(score *apd.Decimal, p *pgs.PGS, ppl string, precisionLimit uint32, known map[int]uint8) *DP {
	s := &DP{
		//target:  score,
		target:  ScoreToTarget(score, p),
		p:       p,
		stats:   p.PopulationStats[ppl],
		ppl:     ppl,
		rounder: newRounder(precisionLimit),
		known:   known,
	}
	return s
}

func (dp *DP) SolveProbabilistic() map[string][]uint8 {
	fmt.Println("Solving probabilistically")
	numSegments := 2
	_, target, roundingError := dp.getTargetAndWeightsAsInts()

	state := dp.BuildProbabilisticState(numSegments)
	targets := make([]int64, 0)
	for w := target; w > target-roundingError-1; w-- {
		targets = append(targets, w)
	}

	solutionHeap := newGenHeap()
	dp.probabilisticMitM(numSegments, state.Nodes, state.Betas, targets, solutionHeap)

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.Ploidy, dp.p.EffectAlleles, dp.known)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func (dp *DP) SolveDeterministic() map[string][]uint8 {
	fmt.Println("Solving deterministically")
	numSegments := 2
	weights, target, roundingError := dp.getTargetAndWeightsAsInts()
	indices := splitIndices(dp.p.Loci, dp.known, numSegments)
	betas := make([]map[uint8]int64, numSegments)
	for i := 0; i < numSegments; i++ {
		betas[i] = makeBetaMap(weights, indices[i])
	}

	maxTotalPositive, maxTotalNegative := GetMaxTotal(betas)
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive-roundingError
	tables := make([]map[int64][]uint8, numSegments)
	var wg sync.WaitGroup
	for i := 0; i < numSegments; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			tables[i] = calculateSubsetSumTable(betas[i], indices[i], upper, lower)
		}(i)
	}
	wg.Wait()

	targets := make([]int64, 0)
	for w := target; w > target-roundingError-1; w-- {
		targets = append(targets, w)
	}

	solutionHeap := newGenHeap()
	dp.deterministicMitM(numSegments, indices, tables, betas, targets, solutionHeap)

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.Ploidy, dp.p.EffectAlleles, dp.known)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func (dp *DP) PrepareData(numSeg int) ([]int64, int64, int64, [][]int, []map[uint8]int64) {
	weights, target, roundingError := dp.getTargetAndWeightsAsInts()
	indices := splitIndices(dp.p.Loci, dp.known, numSeg)
	betas := make([]map[uint8]int64, numSeg)
	for i := 0; i < numSeg; i++ {
		betas[i] = makeBetaMap(weights, indices[i])
	}
	return weights, target, roundingError, indices, betas
}

func (dp *DP) BuildProbabilisticState(numSegments int) *ProbabilisticState {
	_, target, roundingError, indices, betas := dp.PrepareData(numSegments)
	maxTotalPositive, maxTotalNegative := GetMaxTotal(betas)
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive-roundingError
	nodes := make([]map[int64]*Node, numSegments)
	var wg sync.WaitGroup
	for i := 0; i < numSegments; i++ {
		wg.Add(1)
		go func(i int) {
			defer wg.Done()
			nodes[i] = calculateSubsetSumTableWithLikelihood(betas[i], indices[i], upper, lower, dp.stats, dp.p.EffectAlleles)
		}(i)
	}
	wg.Wait()
	return newProbabilisticState(nodes, betas)
}

func (dp *DP) probabilisticMitM(numSegments int, tables []map[int64]*Node, betas []map[uint8]int64, targets []int64,
	solHeap *genHeap) {
	fmt.Printf("Table 0 length: %d\n", len(tables[0]))
	fmt.Printf("Table 1 length: %d\n", len(tables[1]))
	matchHeapSize := 1000 * dp.p.NumVariants
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
				mheap.addToMatchHeap(tables[0][leftSum].TopLikelihood+tables[1][t-leftSum].TopLikelihood, []int64{leftSum, t - leftSum}, matchHeapSize)
			}
		}
	}

	score := apd.NewBigInt(0)
	for _, mtch := range *mheap {
		var solution []uint8
		for i := 0; i < numSegments; i++ {
			solution = append(solution, backtrackWithPointers(mtch.sums[i], tables[i], betas[i])...)
		}
		if dp.rounder.RoundedMode {
			score.SetUint64(0)
			score.Set(lociToScore(solution, dp.rounder.ScaledWeights))
			if score.Cmp(dp.rounder.ScaledTarget) != 0 {
				continue
			}
		}
		solHeap.addToGenHeap(mtch.likelihood, solution, params.HeapSize)
	}
}

func (dp *DP) deterministicMitM(numSegments int, indices [][]int, tables []map[int64][]uint8, betas []map[uint8]int64,
	targets []int64, solHeap *genHeap) {
	fmt.Printf("Table 0 length: %d\n", len(tables[0]))
	fmt.Printf("Table 1 length: %d\n", len(tables[1]))
	step := len(tables[0]) / 10
	if step == 0 {
		step = 1
	}
	var i, s int
	var ok bool
	var leftSum int64
	tmp := make([][]uint8, 0)
	// we need to save only the backtracking results on the right,
	// because each left sum is used only once.
	backtracked := make(map[int64][][]uint8)
	halfSums := make([][]int64, numSegments)
	halfSols := make([][][]uint8, numSegments)
	for i = 0; i < numSegments; i++ {
		halfSums[i] = make([]int64, 0)
		halfSols[i] = make([][]uint8, 0)
	}
	halfSums[0] = append(halfSums[0], 0)
	for leftSum = range tables[0] {
		if s++; s%step == 0 {
			fmt.Printf("Progress: %d%%\n", s*10/step)
		}
		halfSums[1] = halfSums[1][:0]
		for _, t := range targets {
			if _, ok = tables[1][t-leftSum]; ok {
				halfSums[1] = append(halfSums[1], t-leftSum)
			}
		}
		if len(halfSums[1]) == 0 {
			continue
		}
		halfSums[0][0] = leftSum
		// backtrack partial solutions
		for i = 0; i < numSegments; i++ {
			halfSols[i] = halfSols[i][:0]
			for _, halfSum := range halfSums[i] {
				if i == 1 && dp.rounder.RoundedMode {
					if _, ok = backtracked[halfSum]; ok {
						halfSols[i] = append(halfSols[i], backtracked[halfSum]...)
						continue
					}
				}
				tmp = backtrackFromSum(halfSum, tables[i], betas[i])
				if i == 1 && dp.rounder.RoundedMode {
					backtracked[halfSum] = make([][]uint8, 0)
					for _, sol := range tmp {
						backtracked[halfSum] = append(backtracked[halfSum], make([]uint8, len(sol)))
						copy(backtracked[halfSum][len(backtracked[halfSum])-1], sol)
					}
				}
				halfSols[i] = append(halfSols[i], tmp...)
			}
		}
		// combine partial solutions
		dp.combinePartials(0, numSegments, indices, make([]uint8, 0), 0, apd.NewBigInt(0), halfSols, solHeap,
			params.HeapSize)
	}
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumTableWithLikelihood(betas map[uint8]int64, indices []int, upperBound, lowerBound int64,
	stats *pgs.Statistics, effectAlleles []uint8) map[int64]*Node {
	// Fill out the table using dynamic programming
	table := make(map[int64]*Node)
	// the fitness of all the snps being 0
	allNonEffectLikelihood := calculateLociLikelihood([]uint8{}, indices, stats.AF, effectAlleles)
	//fmt.Printf("Indices %v\n", indices)
	//fmt.Printf("All zero likelihood: %f\n", allNonEffectLikelihood)
	table[0] = newNode(math.MaxUint8, allNonEffectLikelihood)
	existingSums := make([]int64, 1)
	existingSums[0] = 0
	updates := make([]Update, 0)
	newSums := make([]int64, 0)

	var ok bool
	var k, nextPtr uint8
	//var nextBinIdx, nextBinCount uint8
	var effectAllele, referenceAllele uint8
	var nextLikelihood float32
	//var nextChi float32
	var prevSum, nextSum, weight int64
	for i, pos := range indices {
		fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[uint8(pos)] > 0 {
			lowerBound += pgs.Ploidy * betas[uint8(pos)]
		} else {
			upperBound += pgs.Ploidy * betas[uint8(pos)]
		}
		updates = updates[:0]
		newSums = newSums[:0]
		effectAllele = effectAlleles[pos]
		referenceAllele = ^effectAllele & 1
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.Ploidy; k++ {
				nextPtr = pgs.Ploidy*uint8(pos) + k - 1
				weight = betas[uint8(pos)] * int64(k)
				nextSum = prevSum + weight
				if nextSum < lowerBound || nextSum > upperBound {
					continue
				}
				nextLikelihood = table[prevSum].TopLikelihood + float32(k)*AfToLikelihood(stats.AF[pos][effectAllele]) -
					float32(k)*AfToLikelihood(stats.AF[pos][referenceAllele]) + float32(pgs.Ploidy-k)*AfToLikelihood(pgs.Ploidy)
				if _, ok = table[nextSum]; !ok {
					table[nextSum] = newNode(nextPtr, nextLikelihood)
					newSums = append(newSums, nextSum)
				}
				table[nextSum].Pointers[nextPtr] = table[prevSum].TopPointer
				if nextLikelihood < table[nextSum].TopLikelihood {
					updates = append(updates, newUpdate(nextSum, nextLikelihood, nextPtr, table[prevSum].TopPointer))
				}
			}
		}
		// We postpone the updates to avoid the sums for the same locus stacking up on each other
		for u := range updates {
			// If we had updated with something better before, sum1 + 1k = sum2 +2k, which we do not want
			if table[updates[u].sum].TopLikelihood < updates[u].likelihood {
				continue
			}
			table[updates[u].sum].TopLikelihood = updates[u].likelihood
			table[updates[u].sum].TopPointer = updates[u].forwardPtr
		}
		existingSums = append(existingSums, newSums...)

	}
	return table
}

// We assume that the weights are sorted in ascending order
func calculateSubsetSumTable(betas map[uint8]int64, indices []int, upperBound, lowerBound int64) map[int64][]uint8 {
	// Fill out the table using dynamic programming
	table := make(map[int64][]uint8)
	// add the zero weight
	table[0] = make([]uint8, 0)
	existingSums := make([]int64, 1)
	existingSums[0] = 0
	newSums := make([]int64, 0)
	var k uint8
	var prevSum, nextSum, weight int64
	for _, pos := range indices {
		if betas[uint8(pos)] > 0 {
			lowerBound += pgs.Ploidy * betas[uint8(pos)]
		} else {
			upperBound += pgs.Ploidy * betas[uint8(pos)]
		}
		newSums = newSums[:0]
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.Ploidy; k++ {
				weight = betas[uint8(pos)] * int64(k)
				nextSum = prevSum + weight
				if nextSum < lowerBound || nextSum > upperBound {
					continue
				}
				if _, ok := table[nextSum]; !ok {
					table[nextSum] = make([]uint8, 0)
					newSums = append(newSums, nextSum)
				}
				table[nextSum] = append(table[nextSum], pgs.Ploidy*uint8(pos)+k-1)
			}
		}
		existingSums = append(existingSums, newSums...)
	}
	return table
}

func backtrackWithPointers(sum int64, table map[int64]*Node, weights map[uint8]int64) []uint8 {
	output := make([]uint8, 0)
	nextPtr := table[sum].TopPointer
	newSum := sum
	for sum != 0 && nextPtr != math.MaxUint8 {
		output = append(output, nextPtr)
		newSum -= weights[nextPtr/2] * int64(nextPtr%2+1)
		nextPtr = table[sum].Pointers[nextPtr]
		sum = newSum
	}
	return output
}

func (dp *DP) combinePartials(segmentNum, totalSegments int, indices [][]int, input []uint8, lkl float32, score *apd.BigInt,
	genotypes [][][]uint8, solHeap *genHeap, heapSize int) {
	//var chi, newLkl float32
	var newLkl float32
	if segmentNum == totalSegments {
		if dp.rounder.RoundedMode {
			if score.Cmp(dp.rounder.ScaledTarget) != 0 {
				return
			}
		}
		solHeap.addToGenHeap(lkl, input, heapSize)
		return
	}
	for _, sol := range genotypes[segmentNum] {
		input = append(input, sol...)
		newLkl = calculateLociLikelihood(sol, indices[segmentNum], dp.stats.AF, dp.p.EffectAlleles) + lkl
		if dp.rounder.RoundedMode {
			newScore := apd.NewBigInt(0)
			newScore.Add(score, lociToScore(sol, dp.rounder.ScaledWeights))
			dp.combinePartials(segmentNum+1, totalSegments, indices, input, newLkl, newScore, genotypes, solHeap, heapSize)
		} else {
			dp.combinePartials(segmentNum+1, totalSegments, indices, input, newLkl, score, genotypes, solHeap, heapSize)
		}
		input = input[:len(input)-len(sol)]
	}
}

func backtrackFromSum(sum int64, table map[int64][]uint8, betas map[uint8]int64) [][]uint8 {
	input := make([]uint8, 0)
	output := make([][]uint8, 0)
	backtrack(input, &output, sum, table, betas)
	return output
}

func backtrack(path []uint8, result *[][]uint8, sum int64, table map[int64][]uint8, weights map[uint8]int64) {
	if sum == 0 {
		sol := make([]uint8, len(path))
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

func splitIndices(loci []string, known map[int]uint8, numSegments int) [][]int {
	numUnknownLoci := len(loci) - len(known)
	indices := make([][]int, numSegments)
	for i := 0; i < numSegments; i++ {
		indices[i] = make([]int, 0)
	}
	count := 0
	for i := range loci {
		if _, ok := known[i]; ok {
			continue
		}
		count++
		if count <= numUnknownLoci/numSegments {
			indices[0] = append(indices[0], i)
		} else {
			indices[1] = append(indices[1], i)
		}
	}
	return indices
}

func makeBetaMap(betas []int64, indices []int) map[uint8]int64 {
	bmap := make(map[uint8]int64)
	for _, idx := range indices {
		bmap[uint8(idx)] = betas[idx]
	}
	return bmap
}

func GetMaxTotal(values []map[uint8]int64) (int64, int64) {
	var positive, negative int64 = 0, 0
	for _, row := range values {
		for _, v := range row {
			if v > 0 {
				positive += 2 * v
			} else {
				negative += 2 * v
			}
		}
	}
	return positive, negative
}

func lociToScore(loci []uint8, weights []*apd.BigInt) *apd.BigInt {
	score := apd.NewBigInt(0)
	for _, locus := range loci {
		score.Add(score, weights[locus/pgs.Ploidy])
		if locus%pgs.Ploidy == 1 {
			score.Add(score, weights[locus/pgs.Ploidy])
		}
	}
	return score
}

func lociToGenotype(loci []uint8, total int, efal []uint8, knownAlleles map[int]uint8) []uint8 {
	sol := make([]uint8, total)
	for i := range efal {
		if efal[i] == ReferenceAllele {
			sol[pgs.Ploidy*i] = 1
			sol[pgs.Ploidy*i+1] = 1
		}
		if knwn, ok := knownAlleles[i]; ok {
			switch knwn {
			case 0:
				sol[pgs.Ploidy*i] = 0
				sol[pgs.Ploidy*i+1] = 0
			case 1:
				sol[pgs.Ploidy*i] = 1
				sol[pgs.Ploidy*i+1] = 0
			case 2:
				sol[pgs.Ploidy*i] = 1
				sol[pgs.Ploidy*i+1] = 1
			}
		}
	}
	//
	for _, locus := range loci {
		switch efal[locus/pgs.Ploidy] {
		case AlternativeAllele:
			sol[locus] = 1
			if locus%pgs.Ploidy == 1 {
				sol[locus-1] = 1
			}
		case ReferenceAllele:
			sol[locus] = 0
			if locus%pgs.Ploidy == 1 {
				sol[locus-1] = 0
			}
		}
	}
	return sol
}
