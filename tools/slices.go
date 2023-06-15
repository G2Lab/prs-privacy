package tools

import (
	"math/rand"
	"time"
)

func SampleFromDistribution(distribution []float64) int {
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

func SampleUniform(values []int) int {
	rand.NewSource(time.Now().UnixNano())
	cumulative := 1.0
	r := rand.Float64() * cumulative
	for k := range values {
		r -= cumulative / float64(len(values))
		if r <= 0.0 {
			return values[k]
		}
	}
	return -1
}

func SampleFromMap(distribution map[int]float64) int {
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

func Shuffle(values [][]int, labels []float64) ([][]int, []float64) {
	rand.NewSource(time.Now().UnixNano())
	for i := len(values) - 1; i > 0; i-- {
		j := rand.Intn(i + 1)
		values[i], values[j] = values[j], values[i]
		labels[i], labels[j] = labels[j], labels[i]
	}
	return values, labels
}
