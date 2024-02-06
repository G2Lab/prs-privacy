package solver

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"time"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
	"github.com/nikirill/prs/tools"
)

type Genetic struct {
	ogTarget      *apd.Decimal
	p             *pgs.PGS
	target        float64
	weights       []float64
	minWeight     float64
	roundingError float64
	af            [][]float64
	freqSpec      []float64
	rounder       *Rounder
	solHeap       *individualHeap
}

func NewGenetic(target *apd.Decimal, p *pgs.PGS, ppl string) *Genetic {
	s := &Genetic{
		ogTarget: target,
		p:        p,
		af:       p.PopulationEAF[ppl],
		freqSpec: p.FreqSpec[ppl],
		rounder:  newRounder(),
		solHeap:  newIndividualHeap(),
	}
	return s
}

func (g *Genetic) Solve() map[string][]uint8 {
	rand.NewSource(time.Now().UnixNano())
	g.weights, g.target, g.roundingError = getTargetAndWeightsAsFloats(g.p, g.ogTarget, g.rounder)
	g.minWeight = findAbsMin(g.weights)
	poolSize := g.p.NumVariants * 10
	//poolSize := g.p.NumVariants
	candidates := make([][]uint8, poolSize, int(float64(3/2)*float64(poolSize)))
	// Initialize mutations solutions according to the allele frequency and frequency spectrum in the population
	var err error
	for i := 0; i < len(candidates); i++ {
		candidates[i], err = SampleFromPopulation(g.af)
		if err != nil {
			fmt.Println(err)
			return nil
		}
	}
	//Evaluate individuals
	for k := 0; k < params.ITERATIONS; k++ {
		if k%500 == 0 {
			fmt.Printf("Iteration %d/%d\n", k, params.ITERATIONS)
		}
		deltas := g.calculateDeltas(candidates)
		candidates, deltas = g.checkForSolutions(candidates, deltas)
		children := g.crossover(candidates, iterationToTemperature(k))
		chDeltas := g.calculateDeltas(children)
		children, chDeltas = g.checkForSolutions(children, chDeltas)
		candidates = g.tournament(append(candidates, children...), append(deltas, chDeltas...), poolSize, iterationToTemperature(k))
		candidates = g.mutate(candidates, iterationToTemperature(k))
	}

	solutions := make(map[string][]uint8)
	for _, sol := range g.solHeap.individuals {
		solutions[ArrayToString(sol.sequence)] = sol.sequence
	}
	return solutions
}

func iterationToTemperature(iteration int) float64 {
	return 1 + float64(params.ITERATIONS-iteration)/2000
}

func (g *Genetic) checkForSolutions(population [][]uint8, deltas []float64) ([][]uint8, []float64) {
	var err error
	for i := 0; i < len(deltas); i++ {
		if math.Abs(deltas[i]) <= g.roundingError {
			if g.rounder.RoundedMode {
				score := CalculateBigIntScore(population[i], g.rounder.ScaledWeights)
				if score.Cmp(g.rounder.ScaledTarget) == 0 {
					g.solHeap.PushIfUnseen(population[i], deltaToFitness(deltas[i], 0), params.HeapSize)
				}
			} else {
				g.solHeap.PushIfUnseen(population[i], deltaToFitness(deltas[i], 0), params.HeapSize)
			}
			population[i], err = SampleFromPopulation(g.af)
			if err != nil {
				log.Fatalf("Error resampling in solution check: %v\n", err)
				return nil, nil
			}
			deltas[i] = CalculateFloatScore(population[i], g.weights) - g.target
		}
	}
	return population, deltas
}

func (g *Genetic) calculateDeltas(population [][]uint8) []float64 {
	deltas := make([]float64, len(population))
	for i := range population {
		deltas[i] = CalculateFloatScore(population[i], g.weights) - g.target
	}
	return deltas
}

func (g *Genetic) crossover(population [][]uint8, T float64) [][]uint8 {
	parents := tools.Shuffle(population)
	splice := func(first, second int) []uint8 {
		// Likelihoods/deltas for the possible crossovers in the form of:
		// 0: first second second second ... second second
		// 1: first first second second ... second second
		// 2: first first first second ... second second
		// ...
		// n: first first first first ... first second

		//likelihoods := make([]float64, len(parents[first])/pgs.NumHplt-1)
		//likelihood := CalculateFullSequenceLikelihood(parents[second], g.af)
		//for k := 0; k < len(g.af)-1; k++ {
		//	likelihood += LocusLikelihood(parents[first], k, g.af)
		//	likelihood -= LocusLikelihood(parents[second], k, g.af)
		//	likelihoods[k] = likelihood
		//}
		//fitness := FitnessFromLikelihoods(likelihoods, T)

		deltas := make([]float64, len(parents[first])/pgs.NumHplt-1)
		delta := CalculateFloatScore(parents[second], g.weights) - g.target
		for k := 0; k < len(g.p.Weights)-1; k++ {
			delta += (float64(parents[first][k*pgs.NumHplt]) + float64(parents[first][k*pgs.NumHplt+1])) * g.weights[k]
			delta -= (float64(parents[second][k*pgs.NumHplt]) + float64(parents[second][k*pgs.NumHplt+1])) * g.weights[k]
		}
		fitness := FitnessFromDeltas(deltas, T)

		selectedIndex := tools.SampleFromDistribution(fitness) + 1
		child := make([]uint8, len(parents[first]))
		copy(child[:selectedIndex*pgs.NumHplt], parents[first][:selectedIndex*pgs.NumHplt])
		copy(child[selectedIndex*pgs.NumHplt:], parents[second][selectedIndex*pgs.NumHplt:])
		return child
	}
	offspring := make([][]uint8, 0, len(parents))
	for i := 0; i < (len(parents)/2)*2; i += 2 {
		offspring = append(offspring, splice(i, i+1))
	}
	return offspring
}

func (g *Genetic) tournament(population [][]uint8, deltas []float64, populationSize int, T float64) [][]uint8 {
	fitness := FitnessFromDeltas(deltas, T)
	scores := make([]float64, len(population))
	var opponentIdx int
	for i := range population {
		for k := 0; k < len(population)/10; k++ {
			opponentIdx = rand.Intn(len(population))
			for opponentIdx == i {
				opponentIdx = rand.Intn(len(population))
			}
			if fitness[i] >= fitness[opponentIdx] {
				scores[i]++
			}
		}
	}
	sortBy(population, scores, true)
	return population[:populationSize]
}

func (g *Genetic) mutate(population [][]uint8, T float64) [][]uint8 {
	for j, original := range population {
		population[j] = g.MutateGenome(original, T)
	}
	return population
}

func FitnessFromDeltas(deltas []float64, T float64) []float64 {
	fitness := make([]float64, len(deltas))
	for i := range deltas {
		// smaller the delta, higher the fitness
		//fitness[i] = 1 / math.Exp(math.Pow(deltas[i], 2))
		fitness[i] = deltaToFitness(deltas[i], T)
	}
	return fitness
}

func FitnessFromLikelihoods(likelihoods []float64, T float64) []float64 {
	fitness := make([]float64, len(likelihoods))
	for i, likelihood := range likelihoods {
		// lower the negative likelihood, higher the fitness
		fitness[i] = 1 / likelihood
		//fitness[i] = 1 / math.Exp(likelihood*T)
	}
	return fitness
}

func (g *Genetic) MutateGenome(original []uint8, T float64) []uint8 {
	//indices := shuffleIndicesByLikelihood(original, g.af)
	//if len(indices) != len(original) {
	//	log.Fatalf("Error shuffling indices: %d != %d", len(indices), len(original))
	//}
	indices := make([]int, len(original))
	for i := range indices {
		indices[i] = i
	}
	rand.Shuffle(len(indices), func(i, j int) {
		indices[i], indices[j] = indices[j], indices[i]
	})
	mutations := make([]uint8, len(original))
	probabilities := make([]float64, len(original))
	originalBins := CalculateAlleleFrequency(original, g.af)
	var newIdx, oldIdx int
	var freqChange float64
	for _, i := range indices {
		for _, v := range pgs.GENOTYPES {
			if v == original[i] {
				continue
			}
			//// If new SNP = 0 and original SNP = 1, delta = delta - weight.
			//// If new SNP = 1 and original SNP = 0, delta = delta + weight.
			//newDelta = delta + (float64(v)-float64(original[i]))*g.weights[i/pgs.NumHplt]
			//if math.Abs(newDelta) <= g.roundingError {
			//	precise := false
			//	if g.rounder.RoundedMode {
			//		cp := make([]uint8, len(original))
			//		copy(cp, original)
			//		cp[i] = v
			//		score := CalculateBigIntScore(cp, g.rounder.ScaledWeights)
			//		if score.Cmp(g.rounder.ScaledTarget) == 0 {
			//			precise = true
			//		}
			//	}
			//	if !g.rounder.RoundedMode || precise {
			//		original[i] = v
			//		return original, newDelta
			//	}
			//}
			oldIdx = tools.ValueToBinIdx(g.af[i/pgs.NumHplt][original[i]], g.p.NumSpecBins)
			newIdx = tools.ValueToBinIdx(g.af[i/pgs.NumHplt][v], g.p.NumSpecBins)
			freqChange = g.specShiftFactor(newIdx, float64(v)-float64(original[i]), originalBins) *
				g.specShiftFactor(oldIdx, float64(original[i])-float64(v), originalBins)
			probabilities[i] = 1 / (afToLikelihood(g.af[i/pgs.NumHplt][v]) * freqChange)
			//probabilities[i] = 1 / math.Abs(newDelta)
			//probabilities[i] = 1 / math.Exp(math.Abs(newDelta))
			//// Falling into a local minimum
			//if math.Abs(newDelta) > g.roundingError && math.Abs(newDelta) < g.minWeight {
			//	newDelta = delta
			//}
			//probabilities[i] = deltaToFitness(newDelta, T)
			//fmt.Printf("%f/%f, ", newDelta, probabilities[i])
			mutations[i] = v
		}
	}
	//fmt.Println(probabilities)
	mutationId := tools.SampleFromDistribution(probabilities)
	//newDelta = delta + (float64(mutations[mutationId])-float64(original[mutationId]))*g.weights[mutationId/pgs.NumHplt]
	original[mutationId] = mutations[mutationId]
	return original
}

func (g *Genetic) specShiftFactor(idx int, shift float64, currentSpec []float64) float64 {
	diff := math.Abs(g.freqSpec[idx] - currentSpec[idx])
	newDiff := math.Abs(diff + shift)
	if newDiff < 1 {
		return math.E
	}
	return math.Exp(math.Pow(diff/newDiff, 2))
}

func deltaToFitness(delta float64, T float64) float64 {
	//return 1 / math.Exp(math.Pow(delta, 2)/T)
	return 1 / math.Exp(math.Pow(delta, 2))
}

func CalculateFloatScore(snps []uint8, weights []float64) float64 {
	score := float64(0)
	for i := 0; i < len(snps); i += pgs.NumHplt {
		for j := 0; j < pgs.NumHplt; j++ {
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

// shuffleIndicesByLikelihood shuffles index order by the likelihood of alternative allele.
// It ensures that if there are several solutions away from the current delta, we tend to pick the one with the highest
// likelihood, but it is still randomized.
func shuffleIndicesByLikelihood(original []uint8, eaf [][]float64) []int {
	var lkl float64
	weightedIndices := make([]int, 0)
	for i, freq := range eaf {
		for j := 0; j < pgs.NumHplt; j++ {
			weightedIndices = append(weightedIndices, pgs.NumHplt*i+j)
			lkl = 1 - freq[original[pgs.NumHplt*i+j]]
			for k := 0; k < int(lkl*100); k++ {
				weightedIndices = append(weightedIndices, pgs.NumHplt*i+j)
			}
		}
	}
	rand.Shuffle(len(weightedIndices), func(i, j int) {
		weightedIndices[i], weightedIndices[j] = weightedIndices[j], weightedIndices[i]
	})
	indices := make([]int, 0, len(original))
	unique := make(map[int]struct{}, len(original))
	for _, v := range weightedIndices {
		if _, exists := unique[v]; !exists {
			indices = append(indices, v)
			unique[v] = struct{}{}
		}
	}
	return indices
}
