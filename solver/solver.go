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
	ReferenceAllele = 0
	EffectAllele    = 1
)

type Solver interface {
	Solve() map[string][]uint8
}

// Smaller the negative fitness, the more likely the sequence is
func calculateLikelihoodForSelectedIndices(mutatedLoci []uint16, indices []int, af map[int][]float64) float64 {
	var likelihood float64 = 0
	indexed := make(map[uint16]struct{})
	for _, pos := range mutatedLoci {
		indexed[pos] = struct{}{}
	}
	var single, double bool
	for _, j := range indices {
		_, single = indexed[uint16(pgs.Ploidy*j)]
		_, double = indexed[uint16(pgs.Ploidy*j+1)]
		switch {
		case single:
			likelihood += afToLikelihood(pgs.Ploidy) + afToLikelihood(af[j][0]) + afToLikelihood(af[j][1])
		case double:
			likelihood += afToLikelihood(af[j][1]) * pgs.Ploidy
		default:
			likelihood += afToLikelihood(af[j][0]) * pgs.Ploidy
		}
	}
	return likelihood
}

// Smaller the negative fitness, the more likely the sequence is
func calculateNegativeLikelihood(mutatedLoci []uint16, startIdx, endIdx int, af map[int][]float64) float64 {
	var likelihood float64 = 0
	indexed := make(map[uint16]struct{})
	for _, pos := range mutatedLoci {
		indexed[pos] = struct{}{}
	}
	var single, double bool
	for j := startIdx; j < endIdx; j += pgs.Ploidy {
		_, single = indexed[uint16(j)]
		_, double = indexed[uint16(j+1)]
		switch {
		case single:
			likelihood += afToLikelihood(pgs.Ploidy) + afToLikelihood(af[j/2][0]) + afToLikelihood(af[j/2][1])
			//likelihood += afToLikelihood(af[j/2][0])
			//likelihood += afToLikelihood(af[j/2][1])
		case double:
			likelihood += afToLikelihood(af[j/2][1]) * pgs.Ploidy
		default:
			likelihood += afToLikelihood(af[j/2][0]) * pgs.Ploidy
		}
	}
	return likelihood
}

func CalculateFullSequenceLikelihood(sequence []uint8, af map[int][]float64) float64 {
	likelihood := 0.0
	//if len(sequence) != len(p.Weights)*pgs.Ploidy {
	//	fmt.Printf("Error: sequence length %d does not match the number of variants %d\n", len(sequence), len(p.Weights))
	//	fmt.Println(sequence)
	//	return 0.0
	//}
	for i := 0; i < len(sequence); i += 2 {
		switch {
		case sequence[i] == 0 && sequence[i+1] == 0:
			likelihood += afToLikelihood(af[i/pgs.Ploidy][0]) * pgs.Ploidy
		case sequence[i] == 1 && sequence[i+1] == 1:
			likelihood += afToLikelihood(af[i/pgs.Ploidy][1]) * pgs.Ploidy
		default:
			likelihood += afToLikelihood(pgs.Ploidy) + afToLikelihood(af[i/pgs.Ploidy][0]) + afToLikelihood(af[i/pgs.Ploidy][1])
		}
		//for j := 0; j < pgs.Ploidy; j++ {
		//	likelihood += afToLikelihood(af[i/pgs.Ploidy][sequence[i+j]])
		//	//fitness += afToLikelihood(p.StudyEAF[i/pgs.Ploidy][sequence[i+j]])
		//}
	}
	return likelihood
}

func afToLikelihood(af float64) float64 {
	return -math.Log(af)
}

// SampleFromPopulation samples a individual according to the MAF
func SampleFromPopulation(af map[int][]float64) ([]uint8, error) {
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

func SampleSegmentFromPopulation(start, end int, af [][]float64) ([]uint8, error) {
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

func LocusLikelihood(sequence []uint8, i int, af [][]float64) float64 {
	likelihood := 0.0
	for j := 0; j < pgs.Ploidy; j++ {
		likelihood += afToLikelihood(af[i][sequence[i*pgs.Ploidy+j]])
	}
	return likelihood
}

func AllReferenceAlleleSample(af map[int][]float64) []uint8 {
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

func locusAlreadyExists(v uint16, array []uint16) bool {
	for _, a := range array {
		if a == v || (v%pgs.Ploidy == 0 && a == v+1) || (v%pgs.Ploidy == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func CalculateDecimalScore(ctx *apd.Context, snps []uint8, weights []*apd.Decimal) *apd.Decimal {
	score := apd.New(0, 0)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				ctx.Add(score, score, weights[i/pgs.Ploidy])
			default:
				log.Printf("Invalid alelle value: %d", snps[i+j])
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

func Accuracy(solution []uint8, target []uint8) float64 {
	if len(solution) != len(target) {
		return 0.0
	}
	acc := 0.0
	for i := 0; i < len(solution); i += pgs.Ploidy {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.Ploidy / float64(len(solution))
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
	accuracies := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		accuracies[i] = Accuracy(solution, target)
		i++
	}
	sortBy(flattened, accuracies, true)
	return flattened
}

func SortByLikelihood(solutions map[string][]uint8, af map[int][]float64) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	likelihoods := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		likelihoods[i] = CalculateFullSequenceLikelihood(solution, af)
		i++
	}
	sortBy(flattened, likelihoods, false)
	return flattened
}

func CalculateEffectAlleleSpectrum(sequence []uint8, af map[int][]float64, bins []float64) []float64 {
	spectrum := make([]float64, len(bins))
	for i := 0; i < len(sequence); i += 2 {
		for j := 0; j < pgs.Ploidy; j++ {
			if sequence[i+j] == EffectAllele {
				spectrum[tools.ValueToBinIdx(af[i/pgs.Ploidy][EffectAllele], bins)]++
			}
		}
	}
	return spectrum
}

func SortByLikelihoodAndFrequency(solutions map[string][]uint8, stats *pgs.Statistics) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	laf := make([]float64, len(solutions))
	var chi float64
	for _, solution := range solutions {
		flattened[i] = solution
		chi = ChiSquaredValue(CalculateEffectAlleleSpectrum(solution, stats.AF, stats.FreqBinBounds), stats.FreqSpectrum)
		if chi < 1 {
			chi = math.Sqrt(chi)
		}
		laf[i] = CalculateFullSequenceLikelihood(solution, stats.AF) + chi
		i++
	}
	sortBy(flattened, laf, false)
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

func ChiSquaredValue(observed, expected []float64) float64 {
	var chi float64
	for i := 0; i < len(observed); i++ {
		if expected[i] == 0 {
			continue
		}
		chi += math.Pow(observed[i]-expected[i], 2) / expected[i]
	}
	return chi
}

func CalculateSolutionSpectrumDistance(solution []uint8, stats *pgs.Statistics) float64 {
	return ChiSquaredValue(CalculateEffectAlleleSpectrum(solution, stats.AF, stats.FreqBinBounds), stats.FreqSpectrum)
}

func CalculateTwoSpectrumDistance(spectrum1, spectrum2 []float64) float64 {
	return ChiSquaredValue(spectrum1, spectrum2)
}

func IncrementObservedInSpectrum(prevBinCount, expectedBinCount float64) float64 {
	return (2*(prevBinCount-expectedBinCount) + 1) / expectedBinCount
}

func CombineLikelihoodAndChiSquared(likelihood, chiSquared float64) float64 {
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

//func SampleMaxMinScores(segmentStart, segmentEnd, numSamples int, betas map[uint16]int64, af [][]float64) (int64, int64) {
//	var err error
//	var sample []uint8
//	var score, maxScore, minScore int64
//	for i := 0; i < numSamples; i++ {
//		sample, err = SampleSegmentFromPopulation(segmentStart, segmentEnd, af)
//		if err != nil {
//			log.Fatalf("Error sampling segment: %v", err)
//		}
//		score = genotypeToScore(segmentStart, segmentEnd, sample, betas)
//		if score > maxScore {
//			maxScore = score
//		}
//		if score < minScore {
//			minScore = score
//		}
//	}
//	return maxScore, minScore
//}
