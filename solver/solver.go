package solver

import (
	"errors"
	"fmt"
	"log"
	"math"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

const (
	ReferenceAllele   = 0
	AlternativeAllele = 1
)

const (
	UseLikelihood = iota
	UseSpectrum
	UseLikelihoodAndSpectrum
)

//type Solver interface {
//	SolveFromScratchProbabilistic() map[string][]uint8
//	SolveFromSavedProbabilistic() map[string][]uint8
//}

func ScoreToTarget(score *apd.Decimal, p *pgs.PGS) *apd.Decimal {
	target := new(apd.Decimal)
	multiplier := new(apd.Decimal).SetInt64(int64(len(p.Loci) * pgs.Ploidy))
	_, err := p.Context.Mul(target, score, multiplier)
	if err != nil {
		log.Printf("Error caluclating target from the score: %v", err)
		return nil
	}
	return target
}

// Smaller the negative fitness, the more likely the sequence is
func calculateLociLikelihood(mutatedLoci []uint8, indices []int, af map[int][]float32, efal []uint8) float32 {
	var likelihood float32 = 0
	indexed := make(map[uint8]struct{})
	for _, pos := range mutatedLoci {
		indexed[pos] = struct{}{}
	}
	var single, double bool
	for _, j := range indices {
		_, single = indexed[uint8(pgs.Ploidy*j)]
		_, double = indexed[uint8(pgs.Ploidy*j+1)]
		switch {
		case single:
			likelihood += afToLikelihood(pgs.Ploidy) + afToLikelihood(af[j][0]) + afToLikelihood(af[j][1])
		case double:
			likelihood += afToLikelihood(af[j][efal[j]]) * pgs.Ploidy
		default:
			likelihood += afToLikelihood(af[j][^efal[j]&1]) * pgs.Ploidy
		}
	}
	return likelihood
}

func CalculateFullSequenceLikelihood(sequence []uint8, af map[int][]float32, efal []uint8) float32 {
	var likelihood float32 = 0.0
	var effect, other uint8
	for i := 0; i < len(sequence); i += 2 {
		effect = efal[i/pgs.Ploidy]
		other = ^effect & 1
		switch {
		case sequence[i] == effect && sequence[i+1] == effect:
			likelihood += afToLikelihood(af[i/pgs.Ploidy][effect]) * pgs.Ploidy
		case sequence[i] == other && sequence[i+1] == other:
			likelihood += afToLikelihood(af[i/pgs.Ploidy][other]) * pgs.Ploidy
		default:
			likelihood += afToLikelihood(pgs.Ploidy) + afToLikelihood(af[i/pgs.Ploidy][effect]) +
				afToLikelihood(af[i/pgs.Ploidy][other])
		}
	}
	return likelihood
}

func afToLikelihood(af float32) float32 {
	return float32(-math.Log(float64(af)))
}

// SampleFromPopulation samples an individual according to the MAF
func SampleFromPopulation(af map[int][]float32) ([]uint8, error) {
	sample := make([]uint8, len(af)*pgs.Ploidy)
	// Initial sample based on individual priors
	for i := range af {
		for j := 0; j < pgs.Ploidy; j++ {
			ind := tools.SampleFromDistribution(af[i])
			sample[i*pgs.Ploidy+j] = pgs.GENOTYPES[ind]
			if sample[i*pgs.Ploidy+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func SampleSegmentFromPopulation(start, end int, af [][]float32) ([]uint8, error) {
	sample := make([]uint8, (end-start)*pgs.Ploidy)
	// Initial sample based on individual priors
	for i := start; i < end; i++ {
		for j := 0; j < pgs.Ploidy; j++ {
			ind := tools.SampleFromDistribution(af[i])
			sample[(i-start)*pgs.Ploidy+j] = pgs.GENOTYPES[ind]
			if sample[(i-start)*pgs.Ploidy+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func LocusLikelihood(sequence []uint8, i int, af [][]float32) float32 {
	var likelihood float32 = 0.0
	for j := 0; j < pgs.Ploidy; j++ {
		likelihood += afToLikelihood(af[i][sequence[i*pgs.Ploidy+j]])
	}
	return likelihood
}

func AllReferenceAlleleSample(af map[int][]float32) []uint8 {
	sample := make([]uint8, 2*len(af))
	for i := 0; i < len(af); i++ {
		if af[i][0] > 0.5 {
			sample[2*i] = 0
			sample[2*i+1] = 0
		} else {
			sample[2*i] = 1
			sample[2*i+1] = 1
		}
	}
	return sample
}

func locusAlreadyExists(v uint8, array []uint8) bool {
	for _, a := range array {
		if a == v || (v%pgs.Ploidy == 0 && a == v+1) || (v%pgs.Ploidy == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func CalculateDecimalScore(ctx *apd.Context, snps []uint8, weights []*apd.Decimal, efal []uint8) *apd.Decimal {
	score := apd.New(0, 0)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
			switch snps[i+j] == efal[i/pgs.Ploidy] {
			case true:
				ctx.Add(score, score, weights[i/pgs.Ploidy])
			case false:
				continue
			}
		}
	}
	return score
}

func CalculateBigIntScore(snps []uint8, weights []*apd.BigInt) *apd.BigInt {
	score := apd.NewBigInt(0)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				score.Add(score, weights[i/pgs.Ploidy])
			default:
				log.Printf("Invalid alelle value: %d", snps[i+j])
			}
		}
	}
	return score
}

func Accuracy(solution []uint8, target []uint8) float32 {
	if len(solution) != len(target) {
		return 0.0
	}
	var acc float32 = 0.0
	for i := 0; i < len(solution); i += pgs.Ploidy {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.Ploidy / float32(len(solution))
}

func ArrayToString(array []uint8) string {
	str := make([]string, len(array)/2)
	for i := 0; i < len(array); i += pgs.Ploidy {
		str[i/pgs.Ploidy] = fmt.Sprint(array[i] + array[i+1])
	}
	return strings.Join(str, "")
}

func SortByAccuracy(solutions map[string][]uint8, target []uint8) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	accuracies := make([]float32, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		accuracies[i] = Accuracy(solution, target)
		i++
	}
	sortBy(flattened, accuracies, true)
	return flattened
}

func SortByLikelihood(solutions map[string][]uint8, af map[int][]float32, efal []uint8) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	likelihoods := make([]float32, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		likelihoods[i] = CalculateFullSequenceLikelihood(solution, af, efal)
		i++
	}
	sortBy(flattened, likelihoods, false)
	return flattened
}

func CalculateSequenceEASpectrum(sequence []uint8, af map[int][]float32, bins []float32, effectAlleles []uint8) []float32 {
	spectrum := make([]float32, len(bins))
	var binIdx int
	for i := 0; i < len(sequence); i += 2 {
		binIdx = tools.ValueToBinIdx(af[i/pgs.Ploidy][effectAlleles[i/pgs.Ploidy]], bins)
		for j := 0; j < pgs.Ploidy; j++ {
			if sequence[i+j] == effectAlleles[i/pgs.Ploidy] {
				spectrum[binIdx]++
			}
		}
	}
	return spectrum
}

func CalculateLociEASpectrum(mutatedLoci []uint8, af map[int][]float32, bins []float32, effectAlleles []uint8) []float32 {
	spectrum := make([]float32, len(bins))
	var binIdx int
	for _, loc := range mutatedLoci {
		binIdx = tools.ValueToBinIdx(af[int(loc)/pgs.Ploidy][effectAlleles[int(loc)/pgs.Ploidy]], bins)
		spectrum[binIdx] += float32(loc%pgs.Ploidy) + 1
	}
	return spectrum
}

func SortByLikelihoodAndFrequency(solutions map[string][]uint8, stats *pgs.Statistics, effectAlleles []uint8, sorting uint8) [][]uint8 {
	flattened := make([][]uint8, len(solutions))
	fitness := make([]float32, len(solutions))
	i := 0
	for _, solution := range solutions {
		flattened[i] = solution
		switch sorting {
		case UseLikelihood:
			fitness[i] = CalculateFullSequenceLikelihood(solution, stats.AF, effectAlleles)
		case UseSpectrum:
			fitness[i] = CalculateSolutionSpectrumDistance(solution, stats, effectAlleles)
		case UseLikelihoodAndSpectrum:
			fitness[i] = CalculateFullSequenceLikelihood(solution, stats.AF, effectAlleles) +
				CalculateSolutionSpectrumDistance(solution, stats, effectAlleles)
		}
		i++
	}
	sortBy(flattened, fitness, false)
	return flattened
}

func sortBy[T, P params.Ordered](items [][]T, properties []P, reverse bool) {
	for i := 0; i < len(items)-1; i++ {
		for j := i + 1; j < len(items); j++ {
			switch reverse {
			case true:
				if properties[i] < properties[j] {
					items[i], items[j] = items[j], items[i]
					properties[i], properties[j] = properties[j], properties[i]
				}
			case false:
				if properties[i] > properties[j] {
					items[i], items[j] = items[j], items[i]
					properties[i], properties[j] = properties[j], properties[i]
				}
			}
		}
	}
}

func ChiSquaredValue(observed, expected []float32) float32 {
	var chi float32
	for i := 0; i < len(observed); i++ {
		if expected[i] == 0 {
			continue
		}
		chi += float32(math.Pow(float64(observed[i])-float64(expected[i]), 2)) / expected[i]
	}
	return chi
}

func CalculateSolutionSpectrumDistance(solution []uint8, stats *pgs.Statistics, effectAlleles []uint8) float32 {
	return ChiSquaredValue(CalculateSequenceEASpectrum(solution, stats.AF, stats.FreqBinBounds, effectAlleles), stats.FreqSpectrum)
}

func CalculateTwoSpectrumDistance(spectrum1, spectrum2 []float32) float32 {
	return ChiSquaredValue(spectrum1, spectrum2)
}

func IncrementObservedInSpectrum(increment, prevBinCount, expectedBinCount float32) float32 {
	return (2*increment*(prevBinCount-expectedBinCount) + float32(math.Pow(float64(increment), 2))) / expectedBinCount
}

func CombineLikelihoodAndChiSquared(likelihood, chiSquared float32) float32 {
	return likelihood + chiSquared
}

func findAbsMin(values []float64) float64 {
	minV := math.Abs(values[0])
	for _, v := range values {
		if math.Abs(v) < minV {
			minV = math.Abs(v)
		}
	}
	return minV
}

//func SampleMaxMinScores(segmentStart, segmentEnd, numSamples int, Betas map[uint8]int64, af [][]float32) (int64, int64) {
//	var err error
//	var sample []uint8
//	var score, maxScore, minScore int64
//	for i := 0; i < numSamples; i++ {
//		sample, err = SampleSegmentFromPopulation(segmentStart, segmentEnd, af)
//		if err != nil {
//			log.Fatalf("Error sampling segment: %v", err)
//		}
//		score = genotypeToScore(segmentStart, segmentEnd, sample, Betas)
//		if score > maxScore {
//			maxScore = score
//		}
//		if score < minScore {
//			minScore = score
//		}
//	}
//	return maxScore, minScore
//}
