package solver

import (
	"fmt"
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
func calculateNegativeLikelihood(mutatedLoci []uint16, startIdx, endIdx int, p *pgs.PGS) float64 {
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
			likelihood += mafToLikelihood(p.Maf[j/2][0])
			likelihood += mafToLikelihood(p.Maf[j/2][1])
		case double:
			likelihood += mafToLikelihood(p.Maf[j/2][1]) * pgs.NumHplt
		default:
			likelihood += mafToLikelihood(p.Maf[j/2][0]) * pgs.NumHplt
		}
	}
	return likelihood
}

func mafToLikelihood(maf float64) float64 {
	return -math.Log(maf)
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
	sortBy(flattened, accuracies)
	return flattened
}

func SortByLikelihood(solutions map[string][]uint8, p *pgs.PGS) [][]uint8 {
	i := 0
	flattened := make([][]uint8, len(solutions))
	likelihoods := make([]float64, len(solutions))
	for _, solution := range solutions {
		flattened[i] = solution
		likelihoods[i] = p.CalculateSequenceLikelihood(solution)
		i++
	}
	sortBy(flattened, likelihoods)
	return flattened
}

func sortBy[T, P params.Ordered](items [][]T, properties []P) {
	for i := 0; i < len(items)-1; i++ {
		for j := i + 1; j < len(items); j++ {
			if properties[i] < properties[j] {
				items[i], items[j] = items[j], items[i]
				properties[i], properties[j] = properties[j], properties[i]
			}
		}
	}
}
