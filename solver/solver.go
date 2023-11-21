package solver

import (
	"fmt"
	"github.com/nikirill/prs/params"
	"log"
	"strings"

	"github.com/nikirill/prs/pgs"
)

type Solver interface {
	Solve(numThreads int) map[string][]uint8
}

//func CalculateScore(snps []int, weights []float64) float64 {
//	score := 0.0
//	for i := 0; i < len(snps); i++ {
//			switch snps[i] {
//			case 0:
//				continue
//			case 1:
//				score += weights[i]
//			case 2:
//				score += weights[i] + weights[i]
//			default:
//				log.Printf("Invalid alelle value: %d", snps[i])
//			}
//	}
//	return score
//}

func CalculateScore(snps []uint8, weights []float64) float64 {
	score := 0.0
	for i := 0; i < len(snps); i += pgs.NumHaplotypes {
		for j := 0; j < pgs.NumHaplotypes; j++ {
			switch snps[i+j] {
			case 0:
				continue
			case 1:
				score += weights[i/2]
			default:
				log.Printf("Invalid alelle value: %d", snps[i+j])
			}
		}
	}
	return score
}

//func Accuracy(solution []int, target []int) float64 {
//	if len(solution) != len(target) {
//		return 0.0
//	}
//	acc := 0.0
//	for i := 0; i < len(solution); i++ {
//		if solution[i] == target[i] {
//			acc++
//		}
//	}
//	return acc / float64(len(solution))
//}

func Accuracy(solution []uint8, target []uint8) float64 {
	if len(solution) != len(target) {
		return 0.0
	}
	acc := 0.0
	for i := 0; i < len(solution); i += pgs.NumHaplotypes {
		if solution[i]+solution[i+1] == target[i]+target[i+1] {
			acc++
		}
	}
	return acc * pgs.NumHaplotypes / float64(len(solution))
}

func sortInts(positions []int, values []int) {
	for i := 0; i < len(values)-1; i++ {
		for j := i + 1; j < len(values); j++ {
			if values[i] > values[j] {
				values[i], values[j] = values[j], values[i]
				positions[i], positions[j] = positions[j], positions[i]
			}
		}
	}
}

func getTriplets(nums []int) [][]int {
	triplets := make([][]int, 0)

	for i := 0; i < len(nums)-2; i++ {
		for j := i + 1; j < len(nums)-1; j++ {
			for k := j + 1; k < len(nums); k++ {
				triplet := []int{nums[i], nums[j], nums[k]}
				triplets = append(triplets, triplet)
			}
		}
	}

	return triplets
}

func ArrayToString(array []uint8) string {
	str := make([]string, len(array)/2)
	for i := 0; i < len(array); i += pgs.NumHaplotypes {
		str[i/2] = fmt.Sprint(array[i] + array[i+1])
	}
	return strings.Join(str, "")
}

func arrayToStringDiploid(array []int) string {
	str := make([]string, len(array))
	for i := 0; i < len(array); i++ {
		str[i] = fmt.Sprint(array[i])
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

func diploidToSum(diploid []int) []int {
	sum := make([]int, len(diploid)/pgs.NumHaplotypes)
	for i := 0; i < len(diploid); i += pgs.NumHaplotypes {
		sum[i/pgs.NumHaplotypes] = diploid[i] + diploid[i+1]
	}
	return sum
}

func findMin(values []float64) float64 {
	minV := values[0]
	for _, v := range values {
		if v < minV {
			minV = v
		}
	}
	return minV
}
