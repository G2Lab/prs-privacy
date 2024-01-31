package solver

import (
	"errors"
	"fmt"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
)

type Solver interface {
	Solve() map[string][]uint8
}

// Smaller the negative likelihood, the more likely the sequence is
func calculateNegativeLikelihood(mutatedLoci []uint16, startIdx, endIdx int, af [][]float64) float64 {
	var likelihood float64 = 0
	indexed := make(map[uint16]struct{})
	for _, pos := range mutatedLoci {
		indexed[pos] = struct{}{}
	}
	var single, double bool
	for j := startIdx; j < endIdx; j += pgs.NumHplt {
		_, single = indexed[uint16(j)]
		_, double = indexed[uint16(j+1)]
		switch {
		case single:
			likelihood += afToLikelihood(af[j/2][0])
			likelihood += afToLikelihood(af[j/2][1])
		case double:
			likelihood += afToLikelihood(af[j/2][1]) * pgs.NumHplt
		default:
			likelihood += afToLikelihood(af[j/2][0]) * pgs.NumHplt
		}
	}
	return likelihood
}

func CalculateFullSequenceLikelihood(sequence []uint8, af [][]float64) float64 {
	likelihood := 0.0
	//if len(sequence) != len(p.Weights)*pgs.NumHplt {
	//	fmt.Printf("Error: sequence length %d does not match the number of variants %d\n", len(sequence), len(p.Weights))
	//	fmt.Println(sequence)
	//	return 0.0
	//}
	for i := 0; i < len(sequence); i += 2 {
		for j := 0; j < pgs.NumHplt; j++ {
			likelihood += afToLikelihood(af[i/pgs.NumHplt][sequence[i+j]])
			//likelihood += afToLikelihood(p.CaseEAF[i/pgs.NumHplt][sequence[i+j]])
		}
	}
	return likelihood
}

func afToLikelihood(af float64) float64 {
	return -math.Log(af)
}

// SampleFromPopulation samples a candidate according to the MAF
func SampleFromPopulation(af [][]float64) ([]uint8, error) {
	sample := make([]uint8, len(af)*pgs.NumHplt)
	// Initial sample based on individual priors
	for i := range af {
		for j := 0; j < pgs.NumHplt; j++ {
			ind := tools.SampleFromDistribution(af[i])
			sample[i*pgs.NumHplt+j] = pgs.GENOTYPES[ind]
			if sample[i*pgs.NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func SampleSegmentFromPopulation(start, end int, af [][]float64) ([]uint8, error) {
	sample := make([]uint8, (end-start)*pgs.NumHplt)
	// Initial sample based on individual priors
	for i := start; i < end; i++ {
		for j := 0; j < pgs.NumHplt; j++ {
			ind := tools.SampleFromDistribution(af[i])
			sample[(i-start)*pgs.NumHplt+j] = pgs.GENOTYPES[ind]
			if sample[(i-start)*pgs.NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func SampleUniform(af [][]float64) ([]uint8, error) {
	sample := make([]uint8, len(af)*pgs.NumHplt)
	// Initial sample based on individual priors
	for i := range af {
		for j := 0; j < pgs.NumHplt; j++ {
			sample[i*pgs.NumHplt+j] = tools.SampleUniform(pgs.GENOTYPES)
			if sample[i*pgs.NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func SnpLikelihood(sequence []uint8, i int, af [][]float64) float64 {
	likelihood := 0.0
	for j := 0; j < pgs.NumHplt; j++ {
		likelihood += afToLikelihood(af[i][sequence[i*pgs.NumHplt+j]])
	}
	return likelihood
}

func AllReferenceAlleleSample(af [][]float64) []uint8 {
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

func SampleMaxMinScores(segmentStart, segmentEnd, numSamples int, betas map[uint16]int64, af [][]float64) (int64, int64) {
	var err error
	var sample []uint8
	var score, maxScore, minScore int64
	for i := 0; i < numSamples; i++ {
		sample, err = SampleSegmentFromPopulation(segmentStart, segmentEnd, af)
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

func locusAlreadyExists(v uint16, array []uint16) bool {
	for _, a := range array {
		if a == v || (v%pgs.NumHplt == 0 && a == v+1) || (v%pgs.NumHplt == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func CalculateDecimalScore(ctx *apd.Context, snps []uint8, weights []*apd.Decimal) *apd.Decimal {
	score := apd.New(0, 0)
	for i := 0; i < len(snps); i += pgs.NumHplt {
		for j := 0; j < pgs.NumHplt; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				ctx.Add(score, score, weights[i/2])
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
	for i := 0; i < len(solution); i += pgs.NumHplt {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.NumHplt / float64(len(solution))
}

func ArrayToString(array []uint8) string {
	str := make([]string, len(array)/2)
	for i := 0; i < len(array); i += pgs.NumHplt {
		str[i/pgs.NumHplt] = fmt.Sprint(array[i] + array[i+1])
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

func SortByLikelihood(solutions map[string][]uint8, af [][]float64) [][]uint8 {
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
