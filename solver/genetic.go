package solver

//import (
//	"fmt"
//	"log"
//	"math"
//	"math/rand"
//	"time"
//
//	"github.com/nikirill/prs/params"
//	"github.com/nikirill/prs/pgs"
//	"github.com/nikirill/prs/tools"
//)
//
//type Genetic struct {
//	target     float64
//	trueGenome []uint8
//	p          *pgs.PGS
//}
//
//func NewGenetic(target float64, tg []uint8, p *pgs.PGS) *Genetic {
//	return &Genetic{
//		target:     target,
//		trueGenome: tg,
//		p:          p,
//	}
//}
//
//func (g *Genetic) Solve(numThreads int) map[string][]uint8 {
//	var err error
//	rand.NewSource(time.Now().UnixNano())
//	poolSize := len(g.p.Variants) * 10
//	//poolSize := len(g.p.Variants)
//	solutions := make(map[string][]uint8)
//	candidates := make([][]uint8, poolSize, int(3/2)*poolSize)
//	// Initialize candidate solutions according to the SNPs likelihood in the population
//	for i := 0; i < len(candidates); i++ {
//		candidates[i], err = g.p.SampleFromPopulation()
//		//candidates[i], err = g.p.SampleUniform()
//		if err != nil {
//			fmt.Println(err)
//			return nil
//		}
//	}
//	//Evaluate candidates
//	for k := 0; k < params.ITERATIONS; k++ {
//		if k%500 == 0 {
//			fmt.Printf("Iteration %d/%d\n", k, params.ITERATIONS)
//			//acc := make([]float64, 0, len(candidates))
//			//for _, candidate := range candidates {
//			//	acc = append(acc, Accuracy(candidate, g.trueGenome))
//			//}
//			//slices.Sort(acc)
//			//fmt.Printf("Accuracy: %v\n", acc[len(acc)-10:])
//		}
//		deltas := g.calculateDeltas(candidates)
//		candidates, deltas, solutions = g.checkForSolutions(candidates, deltas, solutions)
//		candidates, deltas = g.mutate(candidates, deltas, iterationToTemperature(k))
//		candidates, deltas, solutions = g.checkForSolutions(candidates, deltas, solutions)
//		children := g.crossover(candidates, iterationToTemperature(k))
//		chDeltas := g.calculateDeltas(children)
//		children, chDeltas, solutions = g.checkForSolutions(children, chDeltas, solutions)
//		candidates = g.tournament(append(candidates, children...), append(deltas, chDeltas...), poolSize, iterationToTemperature(k))
//	}
//
//	return solutions
//}
//
//func extend(base map[string][]uint8, extension []uint8) {
//	sx := ArrayToString(extension)
//	if _, ok := base[sx]; !ok {
//		base[sx] = make([]uint8, len(extension))
//		copy(base[sx], extension)
//	}
//}
//
//func iterationToTemperature(iteration int) float64 {
//	return 1 + float64(params.ITERATIONS-iteration)/2000
//}
//
//func (g *Genetic) checkForSolutions(population [][]uint8, deltas []float64, solutions map[string][]uint8) ([][]uint8, []float64, map[string][]uint8) {
//	var err error
//	for i := 0; i < len(deltas); i++ {
//		if math.Abs(deltas[i]) <= pgs.ErrorMargin {
//			extend(solutions, population[i])
//			population[i], err = g.p.SampleFromPopulation()
//			if err != nil {
//				log.Fatalf("Error resampling in fitness calculation: %v\n", err)
//				return nil, nil, nil
//			}
//			deltas[i] = CalculateScore(population[i], g.p.Weights) - g.target
//		}
//	}
//	return population, deltas, solutions
//}
//
//func (g *Genetic) calculateDeltas(population [][]uint8) []float64 {
//	deltas := make([]float64, len(population))
//	for i := range population {
//		deltas[i] = CalculateScore(population[i], g.p.Weights) - g.target
//	}
//	return deltas
//}
//
//func (g *Genetic) crossover(population [][]uint8, T float64) [][]uint8 {
//	parents := tools.Shuffle(population)
//	splice := func(first, second int) []uint8 {
//		// Likelihoods for the possible crossovers in the form of:
//		// 0: first second second second ... second second
//		// 1: first first second second ... second second
//		// 2: first first first second ... second second
//		// ...
//		// n: first first first first ... first second
//
//		likelihoods := make([]float64, len(parents[first])/pgs.NumHaplotypes-1)
//		likelihood := g.p.CalculateSequenceLikelihood(parents[second])
//		for k := 0; k < len(g.p.Maf)-1; k++ {
//			likelihood += g.p.SnpLikelihood(parents[first], k)
//			likelihood -= g.p.SnpLikelihood(parents[second], k)
//			likelihoods[k] = likelihood
//		}
//		fitness := FitnessFromLikelihoods(likelihoods, T)
//
//		//deltas := make([]float64, len(parents[first])/pgs.NumHaplotypes-1)
//		//delta := CalculateScore(parents[second], g.p.Weights) - g.target
//		//for k := 0; k < len(g.p.Weights)-1; k++ {
//		//	delta += (float64(parents[first][k*pgs.NumHaplotypes]) + float64(parents[first][k*pgs.NumHaplotypes+1])) * g.p.Weights[k]
//		//	delta -= (float64(parents[second][k*pgs.NumHaplotypes]) + float64(parents[second][k*pgs.NumHaplotypes+1])) * g.p.Weights[k]
//		//}
//		//fitness := FitnessFromDeltas(deltas, T)
//
//		selectedIndex := tools.SampleFromDistribution(fitness) + 1
//		child := make([]uint8, len(parents[first]))
//		copy(child[:selectedIndex*pgs.NumHaplotypes], parents[first][:selectedIndex*pgs.NumHaplotypes])
//		copy(child[selectedIndex*pgs.NumHaplotypes:], parents[second][selectedIndex*pgs.NumHaplotypes:])
//		return child
//	}
//	offspring := make([][]uint8, 0, len(parents))
//	for i := 0; i < (len(parents)/2)*2; i += 2 {
//		offspring = append(offspring, splice(i, i+1))
//	}
//	return offspring
//}
//
//func (g *Genetic) tournament(population [][]uint8, deltas []float64, populationSize int, T float64) [][]uint8 {
//	fitness := FitnessFromDeltas(deltas, T)
//	scores := make([]float64, len(population))
//	var opponentIdx int
//	for i := range population {
//		for k := 0; k < len(population)/10; k++ {
//			opponentIdx = rand.Intn(len(population))
//			for opponentIdx == i {
//				opponentIdx = rand.Intn(len(population))
//			}
//			if fitness[i] >= fitness[opponentIdx] {
//				scores[i]++
//			}
//		}
//	}
//	sortBy(population, scores)
//	return population[:populationSize]
//}
//
//func (g *Genetic) mutate(population [][]uint8, deltas []float64, T float64) ([][]uint8, []float64) {
//	for j, original := range population {
//		population[j], deltas[j] = g.MutateGenome(original, deltas[j], T)
//	}
//	return population, deltas
//}
//
//func FitnessFromDeltas(deltas []float64, T float64) []float64 {
//	fitness := make([]float64, len(deltas))
//	for i := range deltas {
//		// smaller the delta, higher the fitness
//		//fitness[i] = 1 / math.Exp(math.Pow(deltas[i], 2))
//		fitness[i] = deltaToFitness(deltas[i], T)
//	}
//	return fitness
//}
//
//func FitnessFromLikelihoods(likelihoods []float64, T float64) []float64 {
//	minLikelihood := findMin(likelihoods)
//	fitness := make([]float64, len(likelihoods))
//	for i, likelihood := range likelihoods {
//		// higher the (negative) likelihood, higher the fitness
//		//fitness[i] = 1 / math.Abs(likelihood)
//		fitness[i] = math.Abs(likelihood - minLikelihood - 1)
//		//fitness[i] = 1 / math.Exp(-likelihood*T/minLikelihood)
//	}
//	return fitness
//}
//
//func (g *Genetic) MutateGenome(original []uint8, delta float64, T float64) ([]uint8, float64) {
//	indices := g.p.ShuffleIndicesByLikelihood(original)
//	if len(indices) != len(original) {
//		log.Fatalf("Error shuffling indices: %d != %d", len(indices), len(original))
//	}
//	mutations := make([]uint8, len(original))
//	probabilities := make([]float64, len(original))
//	newDelta := 0.0
//	for _, i := range indices {
//		for _, v := range pgs.GENOTYPES {
//			if v == original[i] {
//				continue
//			}
//			// If new SNP = 0 and original SNP = 1, delta = delta - weight.
//			// If new SNP = 1 and original SNP = 0, delta = delta + weight.
//			newDelta = delta + (float64(v)-float64(original[i]))*g.p.Weights[i/pgs.NumHaplotypes]
//			if math.Abs(newDelta) <= pgs.ErrorMargin {
//				original[i] = v
//				return original, newDelta
//			}
//			//probabilities[i] = g.p.Maf[i/pgs.NumHaplotypes][v]
//			//probabilities[i] = 1 / math.Abs(newDelta)
//			//probabilities[i] = 1 / math.Exp(math.Abs(newDelta))
//			probabilities[i] = deltaToFitness(newDelta, T)
//			mutations[i] = v
//		}
//	}
//	mutationId := tools.SampleFromDistribution(probabilities)
//	newDelta = delta + (float64(mutations[mutationId])-float64(original[mutationId]))*g.p.Weights[mutationId/pgs.NumHaplotypes]
//	original[mutationId] = mutations[mutationId]
//	return original, newDelta
//}
//
//func deltaToFitness(delta float64, T float64) float64 {
//	//return 1 / math.Exp(math.Pow(delta, 2)/T)
//	return 1 / math.Exp(math.Pow(delta, 2))
//}
//
//func findComplements(solutions map[string][]uint8, p *pgs.PGS, numThreads int) map[string][]uint8 {
//	// Find which positions have the same weight, hence can be swapped
//	weightGroups := make(map[float64][]int)
//	for i, weight := range p.Weights {
//		if _, ok := weightGroups[weight]; !ok {
//			weightGroups[weight] = make([]int, 1, 2)
//			weightGroups[weight][0] = i
//		} else {
//			weightGroups[weight] = append(weightGroups[weight], i)
//		}
//	}
//	// Get all identical weight positions in a list
//	weightCopies := make([][]int, 0)
//	for weight, positions := range weightGroups {
//		if len(positions) > 1 {
//			weightCopies = append(weightCopies, positions)
//			fmt.Printf("Weight %f in %d positions: %v\n", weight, len(positions), positions)
//		}
//	}
//
//	mutexed := tools.NewMutexMap(solutions)
//	for _, solution := range solutions {
//		//fmt.Printf("Exploring %s\n", ArrayToString(solution))
//		explore(solution, weightCopies, mutexed)
//	}
//	//solSlice := make([][]int, 0, len(solutions))
//	//for _, solution := range solutions {
//	//	solSlice = append(solSlice, solution)
//	//}
//	//var wg sync.WaitGroup
//	//for thread := 0; thread < numThreads; thread++ {
//	//	go func() {
//	//		wg.Add(1)
//	//		defer wg.Done()
//	//		for _, solution := range solSlice[thread*len(solutions)/numThreads : (thread+1)*len(solutions)/numThreads] {
//	//			//fmt.Printf("Exploring %s\n", ArrayToString(solution))
//	//			explore(solution, weightCopies, mutexed)
//	//		}
//	//	}()
//	//}
//	//wg.Wait()
//	return mutexed.RetrieveMapOnly()
//}
//
//func explore(source []uint8, positions [][]int, saver *tools.MutexMap) {
//	var leftPos, rightPos int
//	var leftVal, rightVal uint8
//	//fmt.Printf("Exploring %s\n", ArrayToString(source))
//	//fmt.Printf("Num positions: %d\n", len(positions[0]))
//	if len(positions) > 1 {
//		// Sending further one unchanged version
//		explore(source, positions[1:], saver)
//	}
//	// pairwise search
//	for i := 0; i < len(positions[0])-1; i++ {
//		leftPos = positions[0][i]
//		leftVal = source[leftPos*pgs.NumHaplotypes] + source[leftPos*pgs.NumHaplotypes+1]
//		for j := i + 1; j < len(positions[0]); j++ {
//			rightPos = positions[0][j]
//			rightVal = source[rightPos*pgs.NumHaplotypes] + source[rightPos*pgs.NumHaplotypes+1]
//			switch {
//			case leftVal != rightVal:
//				// first, simple swapping
//				swapped := make([]uint8, len(source))
//				copy(swapped, source)
//				swapped[leftPos*pgs.NumHaplotypes], swapped[rightPos*pgs.NumHaplotypes] =
//					swapped[rightPos*pgs.NumHaplotypes], swapped[leftPos*pgs.NumHaplotypes]
//				swapped[leftPos*pgs.NumHaplotypes+1], swapped[rightPos*pgs.NumHaplotypes+1] =
//					swapped[rightPos*pgs.NumHaplotypes+1], swapped[leftPos*pgs.NumHaplotypes+1]
//				saver.Put(ArrayToString(swapped), swapped)
//				if len(positions) > 1 {
//					explore(swapped, positions[1:], saver)
//				}
//				// second, splitting 2 + 0 into 1 + 1
//				if (leftVal == 2 && rightVal == 0) || (leftVal == 0 && rightVal == 2) {
//					split := make([]uint8, len(source))
//					copy(split, source)
//					split[leftPos*pgs.NumHaplotypes], split[leftPos*pgs.NumHaplotypes+1] = 1, 0
//					split[rightPos*pgs.NumHaplotypes], split[rightPos*pgs.NumHaplotypes+1] = 0, 1
//					saver.Put(ArrayToString(split), split)
//					if len(positions) > 1 {
//						explore(split, positions[1:], saver)
//					}
//				}
//			//	converting 1 + 1 into 2 + 0 and 0 + 2
//			case leftVal == 1 && rightVal == 1:
//				leftPushed, rightPushed := make([]uint8, len(source)), make([]uint8, len(source))
//				copy(leftPushed, source)
//				copy(rightPushed, source)
//				leftPushed[leftPos*pgs.NumHaplotypes], leftPushed[leftPos*pgs.NumHaplotypes+1] = 1, 1
//				leftPushed[rightPos*pgs.NumHaplotypes], leftPushed[rightPos*pgs.NumHaplotypes+1] = 0, 0
//				saver.Put(ArrayToString(leftPushed), leftPushed)
//				rightPushed[leftPos*pgs.NumHaplotypes], rightPushed[leftPos*pgs.NumHaplotypes+1] = 0, 0
//				rightPushed[rightPos*pgs.NumHaplotypes], rightPushed[rightPos*pgs.NumHaplotypes+1] = 1, 1
//				saver.Put(ArrayToString(rightPushed), rightPushed)
//				if len(positions) > 1 {
//					explore(leftPushed, positions[1:], saver)
//					explore(rightPushed, positions[1:], saver)
//				}
//			default:
//				continue
//			}
//		}
//	}
//}
//
////func (s *Solve) recursive(numThreads int) map[string][]int {
////	var wg sync.WaitGroup
////	solutions := NewMutexMap(make(map[string][]int))
////	initial := make([]int, 0, len(s.p.Weights))
////	s.branching(initial,0.0, s.target, 0, solutions, numThreads, &wg)
////	wg.Wait()
////	return solutions.RetrieveMapOnly()
////}
//
////func (s *Solve) branching(current []int, sum, target float64, pos int, solutions *MutexMap, threadsLeft int, wg *sync.WaitGroup) {
////	if math.Abs(sum - target) < pgs.ErrorMargin {
////		if len(current) < len(s.p.Weights) {
////			// pad with zeros
////			current = append(current, make([]int, len(s.p.Weights)-len(current))...)
////		}
////		//fmt.Printf("Found solution: %s\n", arrayToStringDiploid(current))
////		solutions.Put(arrayToStringDiploid(current), current)
////		return
////	}
////	if sum > target || pos >= len(s.p.Weights) {
////		return
////	}
////
////	current = append(current, 0)
////	if threadsLeft > 0 {
////		wg.Add(1)
////		threadsLeft -= 1
////		go func() {
////			s.branching(current, sum, target, pos+1, solutions, threadsLeft, wg)
////			wg.Done()
////		}()
////	} else {
////		s.branching(current, sum, target, pos+1, solutions, threadsLeft, wg)
////	}
////	for snp := 1; snp <= 2; snp++ {
////		branched := make([]int, len(current))
////		copy(branched, current)
////		branched[len(branched)-1] = snp
////		if threadsLeft > 0 {
////			wg.Add(1)
////			threadsLeft -= 1
////			go func() {
////				s.branching(branched, sum+float64(snp)*s.p.Weights[pos], target, pos+1, solutions, threadsLeft, wg)
////				wg.Done()
////			}()
////		} else {
////			s.branching(branched, sum+float64(snp)*s.p.Weights[pos], target, pos+1, solutions, threadsLeft, wg)
////		}
////	}
////}
