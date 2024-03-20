package solver

import (
	"encoding/json"
	"fmt"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"os"
	"path"
	"sort"

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
}

type Node struct {
	Pointers        map[uint16]uint16
	TopPointer      uint16
	TopLikelihood   float64
	TopChiValue     float64
	CurrentBinIdx   uint8
	CurrentBinCount uint8
}

type State struct {
	Nodes []map[int64]*Node
	Betas []map[uint16]int64
}

type Update struct {
	sum        int64
	likelihood float64
	chi        float64
	binIdx     uint8
	binCount   uint8
	forwardPtr uint16
	backPtr    uint16
}

func newNode(ptr uint16, lkl, chi float64, idx, cnt uint8) *Node {
	return &Node{
		Pointers:        make(map[uint16]uint16, 1),
		TopPointer:      ptr,
		TopLikelihood:   lkl,
		TopChiValue:     chi,
		CurrentBinIdx:   idx,
		CurrentBinCount: cnt,
	}
}

func newUpdate(sum int64, likelihood, chi float64, binIdx, binCount uint8, fptr, bptr uint16) Update {
	return Update{
		sum:        sum,
		likelihood: likelihood,
		chi:        chi,
		binIdx:     binIdx,
		binCount:   binCount,
		forwardPtr: fptr,
		backPtr:    bptr,
	}
}

func NewDP(score *apd.Decimal, p *pgs.PGS, ppl string) *DP {
	s := &DP{
		//target:  ScoreToTarget(score, p),
		target:  score,
		p:       p,
		stats:   p.PopulationStats[ppl],
		ppl:     ppl,
		rounder: newRounder(),
	}
	return s
}

func (dp *DP) SolveFromSaved() map[string][]uint8 {
	numSegments := 2
	_, target, roundingError := getTargetAndWeightsAsInts(dp.p, dp.target, dp.rounder)

	state := loadState(fmt.Sprintf("%s-%s.dp", dp.p.PgsID, dp.ppl), numSegments)
	targets := []int64{target}
	if dp.rounder.RoundedMode {
		for w := target - 1; w > target-roundingError; w-- {
			targets = append(targets, w)
		}
	}

	solutionHeap := newGenHeap()
	dp.probabilisticMitM(numSegments, state.Nodes, state.Betas, targets, solutionHeap)

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.Ploidy, dp.p.EffectAlleles)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
}

func (dp *DP) SolveFromScratch() map[string][]uint8 {
	numSegments := 2
	_, target, roundingError := getTargetAndWeightsAsInts(dp.p, dp.target, dp.rounder)

	state := dp.BuildState(numSegments)
	//state := loadState(fmt.Sprintf("%s-%s.dp", dp.p.PgsID, dp.ppl), numSegments)
	targets := make([]int64, 0)
	for w := target; w > target-roundingError-1; w-- {
		targets = append(targets, w)
	}
	//fmt.Printf("Target: %d\n", target)
	//fmt.Printf("Targets: %d\n", targets)

	solutionHeap := newGenHeap()
	dp.probabilisticMitM(numSegments, state.Nodes, state.Betas, targets, solutionHeap)

	solutions := make(map[string][]uint8)
	for _, sol := range *solutionHeap {
		subset := lociToGenotype(sol.mutations, len(dp.p.Weights)*pgs.Ploidy, dp.p.EffectAlleles)
		solutions[ArrayToString(subset)] = subset
	}
	return solutions
	//	tables := make([]map[int64][]uint16, numSegments)
	//	for i := 0; i < numSegments; i++ {
	//		tables[i] = calculateSubsetSumTable(Betas[i], upper, lower)
	//	}
	//	dp.deterministicMitM(numSegments, freqSortedIndices, splitIdxs, tables, Betas, targets, solutionHeap)
}

func (dp *DP) PrepareData(numSeg int) ([]int64, int64, int64, []int, []int, []map[uint16]int64) {
	weights, target, roundingError := getTargetAndWeightsAsInts(dp.p, dp.target, dp.rounder)
	freqSortedIndices := sortByEffectAlleleFreq(dp.stats.AF, dp.p.EffectAlleles)

	var splitIdxs []int
	splitIdxs = []int{0, len(weights) / 2, len(weights)}
	// Make sure that the frequency bins are not split
	preSplitIdx := splitIdxs[1] - 1
	postSplitIdx := splitIdxs[1]
	prevBin := tools.ValueToBinIdx(dp.stats.AF[freqSortedIndices[preSplitIdx]][dp.p.EffectAlleles[freqSortedIndices[preSplitIdx]]], dp.stats.FreqBinBounds)
	nextBin := tools.ValueToBinIdx(dp.stats.AF[freqSortedIndices[postSplitIdx]][dp.p.EffectAlleles[freqSortedIndices[postSplitIdx]]], dp.stats.FreqBinBounds)
	for prevBin == nextBin && preSplitIdx > 0 && postSplitIdx < len(weights) {
		preSplitIdx--
		prevBin = tools.ValueToBinIdx(dp.stats.AF[freqSortedIndices[preSplitIdx]][dp.p.EffectAlleles[freqSortedIndices[preSplitIdx]]], dp.stats.FreqBinBounds)
		if prevBin != nextBin {
			splitIdxs[1] = preSplitIdx + 1
			break
		}
		postSplitIdx++
		nextBin = tools.ValueToBinIdx(dp.stats.AF[freqSortedIndices[postSplitIdx]][dp.p.EffectAlleles[freqSortedIndices[postSplitIdx]]], dp.stats.FreqBinBounds)
		if prevBin != nextBin {
			splitIdxs[1] = postSplitIdx
			break
		}
	}
	betas := make([]map[uint16]int64, numSeg)
	for i := 0; i < numSeg; i++ {
		betas[i] = makeBetaMap(weights, freqSortedIndices, splitIdxs[i], splitIdxs[i+1])
	}
	return weights, target, roundingError, freqSortedIndices, splitIdxs, betas
}

func (dp *DP) BuildState(numSegments int) *State {
	weights, target, roundingError, freqSortedIndices, splitIdxs, betas := dp.PrepareData(numSegments)
	maxTotalPositive, maxTotalNegative := GetMaxTotal(weights)
	upper, lower := target-maxTotalNegative+roundingError, target-maxTotalPositive-roundingError
	state := new(State)
	state.Betas = betas
	state.Nodes = make([]map[int64]*Node, numSegments)
	for i := 0; i < numSegments; i++ {
		state.Nodes[i] = calculateSubsetSumTableWithLikelihood(betas[i], freqSortedIndices[splitIdxs[i]:splitIdxs[i+1]], upper, lower, dp.stats, dp.p.EffectAlleles)
	}
	return state
}

func SaveState(state *State, filepath string) {
	tf, err := os.OpenFile(filepath, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		log.Fatalf("Error opening tables file: %v", err)
	}
	defer tf.Close()
	encoder := json.NewEncoder(tf)
	if err = encoder.Encode(state); err != nil {
		log.Fatalf("Error encoding tables: %v", err)
	}
}

func loadState(filename string, numSegments int) *State {
	filepath := path.Join(params.DataFolder, filename)
	tf, err := os.Open(filepath)
	if err != nil {
		log.Fatalf("Error opening tables file: %v", err)
	}
	defer tf.Close()
	decoder := json.NewDecoder(tf)
	tables := make([]map[int64]*Node, numSegments)
	betas := make([]map[uint16]int64, numSegments)
	state := &State{Nodes: tables, Betas: betas}
	if err = decoder.Decode(state); err != nil {
		log.Fatalf("Error decoding tables: %v", err)
	}
	return state
}

func (dp *DP) probabilisticMitM(numSegments int, tables []map[int64]*Node, betas []map[uint16]int64, targets []int64,
	solHeap *genHeap) {
	matchHeapSize := 1000 * dp.p.NumVariants
	step := len(tables[0]) / 10
	if step == 0 {
		step = 1
	}
	mheap := newMatchHeap()
	var s int
	var ok bool
	var leftSum int64
	fmt.Printf("Table 0 length: %d\n", len(tables[0]))
	fmt.Printf("Table 1 length: %d\n", len(tables[1]))
	for leftSum = range tables[0] {
		if s++; s%step == 0 {
			fmt.Printf("Progress: %d%%\n", s*10/step)
		}
		for _, t := range targets {
			if _, ok = tables[1][t-leftSum]; ok {
				//mheap.addToMatchHeap(tables[0][leftSum].TopLikelihood+tables[1][t-leftSum].TopLikelihood, []int64{leftSum, t - leftSum}, matchHeapSize)
				//mheap.addToMatchHeap(tables[0][leftSum].TopChiValue+tables[1][t-leftSum].TopChiValue, []int64{leftSum, t - leftSum}, matchHeapSize)
				mheap.addToMatchHeap(CombineLikelihoodAndChiSquared(tables[0][leftSum].TopLikelihood, tables[0][leftSum].TopChiValue)+
					CombineLikelihoodAndChiSquared(tables[1][t-leftSum].TopLikelihood, tables[1][t-leftSum].TopChiValue),
					[]int64{leftSum, t - leftSum}, matchHeapSize)
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
			score.SetUint64(0)
			score.Set(lociToScore(solution, dp.rounder.ScaledWeights))
			if score.Cmp(dp.rounder.ScaledTarget) != 0 {
				continue
			}
		}
		solHeap.addToGenHeap(mtch.likelihood, solution, params.HeapSize)
	}
}

func (dp *DP) deterministicMitM(numSegments int, indices []int, splitIdxs []int, tables []map[int64][]uint16, betas []map[uint16]int64,
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
					lkl = calculateLikelihoodForSelectedIndices(combinations[j], indices[splitIdxs[0]:splitIdxs[i+1]], dp.stats.AF, dp.p.EffectAlleles)
					halfSols[i] = append(halfSols[i], newGenotype(combinations[j], lkl))
				}
			}
			//halfSols[i] = make([]*genotype, 0)
			//for _, halfSum := range halfSums[i] {
			//	if _, ok = backtracked[i][halfSum]; !ok {
			//		combinations := backtrackFromSum(halfSum, tables[i], Betas[i])
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
func calculateSubsetSumTableWithLikelihood(betas map[uint16]int64, indices []int, upperBound, lowerBound int64,
	stats *pgs.Statistics, effectAlleles []uint8) map[int64]*Node {
	// Fill out the table using dynamic programming
	table := make(map[int64]*Node)
	// the fitness of all the snps being 0
	allNonEffectLikelihood := calculateLikelihoodForSelectedIndices([]uint16{}, indices, stats.AF, effectAlleles)
	fmt.Printf("Indices %v\n", indices)
	fmt.Printf("All zero likelihood: %f\n", allNonEffectLikelihood)

	firstBucketIdx := tools.ValueToBinIdx(stats.AF[indices[0]][effectAlleles[indices[0]]], stats.FreqBinBounds)
	lastBucketIdx := tools.ValueToBinIdx(stats.AF[indices[len(indices)-1]][effectAlleles[indices[len(indices)-1]]], stats.FreqBinBounds)
	allNonEffectFreqSpec := make([]float64, lastBucketIdx-firstBucketIdx+1)
	freqSpecSegment := stats.FreqSpectrum[firstBucketIdx : lastBucketIdx+1]
	allNonEffectDistance := CalculateTwoSpectrumDistance(allNonEffectFreqSpec, freqSpecSegment)
	fmt.Printf("First bucket %d, last bucket %d, total %d\n", firstBucketIdx, lastBucketIdx, len(stats.FreqSpectrum))
	fmt.Printf("All zero distance: %f\n", allNonEffectDistance)
	// add the zero weight
	table[0] = newNode(math.MaxUint16, allNonEffectLikelihood, allNonEffectDistance, 0, 0)
	existingSums := make([]int64, 1)
	existingSums[0] = 0

	var updates []Update
	var ok bool
	var k, nextPtr uint16
	var nextBinIdx, nextBinCount uint8
	var effectAllele, referenceAllele uint8
	var nextLikelihood, nextChi float64
	var prevSum, nextSum, weight int64
	for i, pos := range indices {
		fmt.Printf("Position %d/%d\n", i, len(betas))
		if betas[uint16(pos)] > 0 {
			lowerBound += pgs.Ploidy * betas[uint16(pos)]
		} else {
			upperBound += pgs.Ploidy * betas[uint16(pos)]
		}
		updates = make([]Update, 0)
		newSums := make([]int64, 0)
		effectAllele = effectAlleles[pos]
		referenceAllele = ^effectAllele & 1
		nextBinIdx = uint8(tools.ValueToBinIdx(stats.AF[pos][effectAllele], stats.FreqBinBounds))
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.Ploidy; k++ {
				nextPtr = pgs.Ploidy*uint16(pos) + k - 1
				weight = betas[uint16(pos)] * int64(k)
				nextSum = prevSum + weight
				if nextSum < lowerBound || nextSum > upperBound {
					continue
				}
				nextLikelihood = table[prevSum].TopLikelihood + float64(k)*afToLikelihood(stats.AF[pos][effectAllele]) -
					float64(k)*afToLikelihood(stats.AF[pos][referenceAllele]) + float64(pgs.Ploidy-k)*afToLikelihood(pgs.Ploidy)
				nextBinCount = table[prevSum].CurrentBinCount + uint8(k)
				if nextBinIdx != table[prevSum].CurrentBinIdx {
					nextBinCount = uint8(k)
				}
				nextChi = table[prevSum].TopChiValue + IncrementObservedInSpectrum(float64(k), float64(nextBinCount)-1, stats.FreqSpectrum[nextBinIdx])
				if _, ok = table[nextSum]; !ok {
					table[nextSum] = newNode(nextPtr, nextLikelihood, nextChi, nextBinIdx, nextBinCount)
					newSums = append(newSums, nextSum)
				}
				table[nextSum].Pointers[nextPtr] = table[prevSum].TopPointer
				if CombineLikelihoodAndChiSquared(nextLikelihood, nextChi) <
					CombineLikelihoodAndChiSquared(table[nextSum].TopLikelihood, table[nextSum].TopChiValue) {
					//if nextLikelihood < table[nextSum].TopLikelihood {
					//if nextChi < table[nextSum].TopChiValue {
					updates = append(updates, newUpdate(nextSum, nextLikelihood, nextChi, nextBinIdx, nextBinCount, nextPtr, table[prevSum].TopPointer))
				}
			}
		}
		// We postpone the updates to avoid the sums for the same locus stacking up on each other
		for u := range updates {
			// If we have already updated with something better, sum1 + 1k = sum2 +2k
			if CombineLikelihoodAndChiSquared(table[updates[u].sum].TopLikelihood, table[updates[u].sum].TopChiValue) <
				CombineLikelihoodAndChiSquared(updates[u].likelihood, updates[u].chi) {
				//if table[updates[u].sum].TopLikelihood < updates[u].likelihood {
				//if table[updates[u].sum].TopChiValue < updates[u].chi {
				continue
			}
			table[updates[u].sum].TopLikelihood = updates[u].likelihood
			table[updates[u].sum].TopChiValue = updates[u].chi
			table[updates[u].sum].CurrentBinIdx = updates[u].binIdx
			table[updates[u].sum].CurrentBinCount = updates[u].binCount
			table[updates[u].sum].TopPointer = updates[u].forwardPtr
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
		//fmt.Printf("Position %d/%d\n", i, len(Betas))
		if betas[pos] > 0 {
			lowerBound += pgs.Ploidy * betas[pos]
		} else {
			upperBound += pgs.Ploidy * betas[pos]
		}
		newSums := make([]int64, 0)
		for _, prevSum = range existingSums {
			for k = 1; k <= pgs.Ploidy; k++ {
				weight = betas[pos] * int64(k)
				nextSum = prevSum + weight
				if nextSum >= lowerBound && nextSum <= upperBound {
					if _, ok := table[nextSum]; !ok {
						table[nextSum] = make([]uint16, 0)
						newSums = append(newSums, nextSum)
					}
					table[nextSum] = append(table[nextSum], pgs.Ploidy*pos+k-1)
				}
			}
		}
		existingSums = append(existingSums, newSums...)
	}
	return table
}

func backtrackWithPointers(sum int64, table map[int64]*Node, weights map[uint16]int64) []uint16 {
	output := make([]uint16, 0)
	nextPtr := table[sum].TopPointer
	newSum := sum
	for sum != 0 && nextPtr != math.MaxUint16 {
		output = append(output, nextPtr)
		newSum -= weights[nextPtr/2] * int64(nextPtr%2+1)
		nextPtr = table[sum].Pointers[nextPtr]
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
		//if locusAlreadyExists(ptr, path) {
		if locusAlreadyExists(ptr, path) || (len(path) > 0 && ptr > path[len(path)-1]) {
			continue
		}
		path = append(path, ptr)
		newSum := sum - weights[ptr/2]*int64(ptr%2+1)
		backtrack(path, result, newSum, table, weights)
		path = path[:len(path)-1]
	}
}

func makeBetaMap(betas []int64, indices []int, start, end int) map[uint16]int64 {
	bmap := make(map[uint16]int64, end-start)
	for i := start; i < end; i++ {
		bmap[uint16(indices[i])] = betas[indices[i]]
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

func lociToScore(loci []uint16, weights []*apd.BigInt) *apd.BigInt {
	score := apd.NewBigInt(0)
	for _, locus := range loci {
		score.Add(score, weights[locus/pgs.Ploidy])
		if locus%pgs.Ploidy == 1 {
			score.Add(score, weights[locus/pgs.Ploidy])
		}
	}
	return score
}

func lociToGenotype(loci []uint16, total int, efal []uint8) []uint8 {
	sol := make([]uint8, total)
	for i := range efal {
		if efal[i] == ReferenceAllele {
			sol[pgs.Ploidy*i] = 1
			sol[pgs.Ploidy*i+1] = 1
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

func sortByEffectAlleleFreq(m map[int][]float64, efal []uint8) []int {
	indices := make([]int, 0, len(m))
	for ind := range m {
		indices = append(indices, ind)
	}
	sort.Slice(indices, func(i, j int) bool {
		return m[indices[i]][efal[indices[i]]] < m[indices[j]][efal[indices[j]]]
	})
	fmt.Println("Sorted indices by effect allele frequency:")
	for i := range indices {
		fmt.Printf("%d:%f ", indices[i], m[indices[i]][efal[indices[i]]])
	}
	fmt.Println()

	return indices
}

func genotypeToScore(start, end int, snps []uint8, betas map[uint16]int64) int64 {
	score := int64(0)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
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
