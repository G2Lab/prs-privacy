package main

import (
	"encoding/csv"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"time"

	"github.com/nikirill/prs/utils"
)

const (
	SCORE      = "score"
	ITERATIONS = 100
	MARGIN     = 0.0001
)

func main() {
	pgs := NewPGS()
	err := pgs.LoadCatalogFile("PGS000073_hmPOS_GRCh38.txt")
	if err != nil {
		log.Println("Error:", err)
		return
	}
	pgs.LoadPriors()

	cohort := NewCohort()
	cohort.CalculatePRS(pgs)
	solution := FindSolution(cohort["NA20543"][SCORE], pgs)
	fmt.Println(pgs.Weights)
	fmt.Println("Guessed:", solution)
	fmt.Print("True:     ")
	for _, location := range pgs.Locations {
		fmt.Printf("%.0f ", cohort["NA20543"][location])
	}
	fmt.Printf("\nGuessed score:%f\n", calculateScore(solution, pgs.Weights))
	fmt.Printf("True score:%f", cohort["NA20543"][SCORE])
	//sortedInd := cohort.SortByScore()
	//err = cohort.SaveScores(sortedInd)
	//if err != nil {
	//	log.Println("Error saving scores:", err)
	//	return
	//}
	//for _, ind := range sortedInd {
	//	fmt.Println(ind, c[ind][SCORE])
	//}
}

type Cohort map[string]map[string]float64

func NewCohort() Cohort {
	c := make(Cohort)
	return c
}

func (c Cohort) CalculatePRS(pgs *PGS) {
	for _, variant := range pgs.Variants {
		chr := variant.GetHmChr()
		position := variant.GetHmPos()
		query, args := utils.IndividualSnpsQuery(chr, position)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			fmt.Println("Error executing bcftools command:", err)
			continue
		}

		samples := strings.Split(string(output), "\t")
		for _, sample := range samples {
			fields := strings.Split(sample, "=")
			if len(fields) != 2 {
				//fmt.Println("Error splitting sample:", sample)
				continue
			}
			individual := fields[0]
			snp := fields[1]
			value, err := utils.SnpToValue(snp)
			if err != nil {
				fmt.Println(err)
				continue
			}
			if _, ok := c[individual]; !ok {
				c[individual] = make(map[string]float64)
			}
			key := fmt.Sprintf("%s:%s", chr, position)
			c[individual][key] = value
			c[individual][SCORE] += value * variant.GetWeight()
		}
	}
}

func (c Cohort) SortByScore() []string {
	sortedInd := make([]string, 0, len(c))
	for ind := range c {
		sortedInd = append(sortedInd, ind)
	}
	for i := 0; i < len(sortedInd)-1; i++ {
		minIndex := i
		for j := i + 1; j < len(sortedInd); j++ {
			if c[sortedInd[j]][SCORE] < c[sortedInd[minIndex]][SCORE] {
				minIndex = j
			}
		}
		if minIndex != i {
			sortedInd[i], sortedInd[minIndex] = sortedInd[minIndex], sortedInd[i]
		}
	}
	return sortedInd
}

func (c Cohort) SaveScores(sortedInd []string) error {
	file, err := os.Create("scores.csv")
	if err != nil {
		return err
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	for _, ind := range sortedInd {
		writer.Write([]string{ind, fmt.Sprintf("%0.17f", c[ind][SCORE])})
	}
	writer.Flush()
	return nil
}

func (c Cohort) LoadScores(filename string) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	reader := csv.NewReader(file)
	for {
		record, err := reader.Read()
		if err != nil {
			return err
		}
		individual := record[0]
		score := record[1]
		c[individual] = make(map[string]float64)
		if v, err := strconv.ParseFloat(score, 64); err == nil {
			c[individual][SCORE] = v
		} else {
			return err
		}
	}
	return nil
}

func FindSolution(targetScore float64, pgs *PGS) []int {
	SNPS := [...]int{0, 1, 2}
	candidates := make([][]int, 10)
	weights := pgs.Weights

	// Initialize candidate solutions according to the SNPs likelihood in the population
	for i := 0; i < len(candidates); i++ {
		candidates[i] = make([]int, len(pgs.Variants))
		for j, loc := range pgs.Locations {
			candidates[i][j] = sample(pgs.GetVariantPriors(loc))
			if candidates[i][j] == -1 {
				fmt.Println("Error sampling candidate")
				return nil
			}
		}
	}

	// Evaluate candidates
	var delta float64
	//mutated := make([][]int, len(candidates))
	for k := 0; k < ITERATIONS; k++ {
		for j, candidate := range candidates {
			delta = calculateScore(candidate, weights) - targetScore
			//fmt.Println(candidate, delta)
			if math.Abs(delta) < MARGIN {
				fmt.Println("Found solution:", candidate)
				return candidate
			}
			// the weights that cover the delta better, get higher probability of being selected
			// big delta -> bigger weights get higher probability.
			// probability = 1 / abs( delta - weight * snp_old + weight * snp_new )
			possibleMutations := make([]int, 0, len(candidate)*(len(SNPS)-1))
			probs := make([]float64, 0)
			for i, snp := range candidate {
				tmp := delta
				tmp -= weights[i] * float64(snp)
				for _, v := range SNPS {
					if v == snp {
						continue
					}
					possibleMutations = append(possibleMutations, v)
					if math.Abs(tmp+weights[i]*float64(v)) < MARGIN {
						//fmt.Println("Found solution:", candidate)
						//fmt.Printf("%d -> %d\n", i, v)
						candidate[i] = v
						return candidate
					}
					probs = append(probs, 1/math.Abs(tmp+weights[i]*float64(v)))
				}
			}
			//mutated[j] = pgs.MutateVariant(candidate, possibleMutations, probs)
			candidates[j] = pgs.MutateVariant(candidate, possibleMutations, probs)
		}
	}
	return nil
}

func sample(distribution map[int]float64) int {
	rand.NewSource(time.Now().UnixNano())
	cumulative := 0.0
	for _, p := range distribution {
		cumulative += p
	}
	r := rand.Float64() * cumulative
	for k, v := range distribution {
		r -= v
		if r <= 0.0 {
			return k
		}
	}
	return -1
}

func calculateScore(snps []int, weights []float64) float64 {
	score := 0.0
	for i, snp := range snps {
		score += float64(snp) * weights[i]
	}
	return score
}
