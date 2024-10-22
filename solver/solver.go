package solver

import (
	"fmt"
	"log"
	"math"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/pgs"
)

const (
	ReferenceAllele   = 0
	AlternativeAllele = 1
)

func ScoreToTarget(score *apd.Decimal, p *pgs.PGS) *apd.Decimal {
	target := new(apd.Decimal)
	multiplier := new(apd.Decimal).SetInt64(int64(len(p.Loci) * pgs.Ploidy))
	_, err := p.Context.Mul(target, score, multiplier)
	if err != nil {
		log.Printf("Error calculating target from the score: %v", err)
		return nil
	}
	precision := p.WeightPrecision
	tf, err := target.Float64()
	if err != nil {
		log.Printf("Error converting target to float: %v", err)
		return nil
	}
	if tf != 0 {
		precision += uint32(len(fmt.Sprintf("%.0f", math.Abs(tf))))
	}
	// Rounding the target to the correct precision
	roundingCtx := apd.BaseContext.WithPrecision(precision)
	roundingCtx.Rounding = apd.RoundHalfUp
	_, err = roundingCtx.Round(target, target)
	if err != nil {
		log.Printf("Error rounding the target: %v", err)
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
			likelihood += AfToLikelihood(pgs.Ploidy) + AfToLikelihood(af[j][0]) + AfToLikelihood(af[j][1])
		case double:
			likelihood += AfToLikelihood(af[j][efal[j]]) * pgs.Ploidy
		default:
			likelihood += AfToLikelihood(af[j][^efal[j]&1]) * pgs.Ploidy
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
			likelihood += AfToLikelihood(af[i/pgs.Ploidy][effect]) * pgs.Ploidy
		case sequence[i] == other && sequence[i+1] == other:
			likelihood += AfToLikelihood(af[i/pgs.Ploidy][other]) * pgs.Ploidy
		default:
			likelihood += AfToLikelihood(pgs.Ploidy) + AfToLikelihood(af[i/pgs.Ploidy][effect]) +
				AfToLikelihood(af[i/pgs.Ploidy][other])
		}
	}
	return likelihood
}

func AfToLikelihood(af float32) float32 {
	return float32(-math.Log2(float64(af)))
}

func locusAlreadyExists(v uint8, array []uint8) bool {
	for _, a := range array {
		if a == v || (v%pgs.Ploidy == 0 && a == v+1) || (v%pgs.Ploidy == 1 && a == v-1) {
			return true
		}
	}
	return false
}

func CalculateDecimalSum(ctx *apd.Context, snps []uint8, weights []*apd.Decimal, efal []uint8) *apd.Decimal {
	sum := apd.New(0, 0)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
			switch snps[i+j] == efal[i/pgs.Ploidy] {
			case true:
				ctx.Add(sum, sum, weights[i/pgs.Ploidy])
			case false:
				continue
			}
		}
	}
	return sum
}

func CalculateBigIntSum(snps []uint8, weights []*apd.BigInt, efal []uint8) *apd.BigInt {
	sum := new(apd.BigInt)
	for i := 0; i < len(snps); i += pgs.Ploidy {
		for j := 0; j < pgs.Ploidy; j++ {
			switch snps[i+j] == efal[i/pgs.Ploidy] {
			case true:
				sum.Add(sum, weights[i/pgs.Ploidy])
			case false:
				continue
			}
		}
	}
	return sum
}

func FindMaxAbsoluteBigInt(ints []*apd.BigInt) int64 {
	var maxWeight int64 = 0
	for _, weight := range ints {
		w := weight.Int64()
		if math.Abs(float64(w)) > float64(maxWeight) {
			maxWeight = int64(math.Abs(float64(w)))
		}
	}
	return maxWeight
}

func CalculateDecimalScore(ctx *apd.Context, snps []uint8, weights []*apd.Decimal, efal []uint8) *apd.Decimal {
	score := CalculateDecimalSum(ctx, snps, weights, efal)
	// Normalize the score by dividing by the number of loci and ploidy
	divisor := new(apd.Decimal).SetInt64(int64(len(weights) * pgs.Ploidy))
	_, err := ctx.Quo(score, score, divisor)
	if err != nil {
		log.Println("Error normalizing the score:", err)
		return nil
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

func SortByLikelihood(solutions map[string][]uint8, stats *pgs.Statistics, effectAlleles []uint8) [][]uint8 {
	flattened := make([][]uint8, len(solutions))
	likelihood := make([]float32, len(solutions))
	i := 0
	for _, solution := range solutions {
		flattened[i] = solution
		likelihood[i] = CalculateFullSequenceLikelihood(solution, stats.AF, effectAlleles)
		i++
	}
	sortBy(flattened, likelihood, false)
	return flattened
}

func sortBy(items [][]uint8, properties []float32, reverse bool) {
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
