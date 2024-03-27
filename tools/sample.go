package tools

import (
	"log"
	"math/rand"
)

func SampleFromDistribution(distribution []float32) int {
	var cumulative float32 = 0.0
	for _, p := range distribution {
		cumulative += p
	}
	r := rand.Float32() * cumulative
	for k, v := range distribution {
		r -= v
		if r <= 0.0 {
			return k
		}
	}
	log.Printf("Unsuccessful sampling from distribution: %v", distribution)
	return -1
}

func SampleUniform(values []uint8) uint8 {
	cumulative := 1.0
	r := rand.Float64() * cumulative
	for k := range values {
		r -= cumulative / float64(len(values))
		if r <= 0.0 {
			return values[k]
		}
	}
	return 255
}

func SampleFromMap(distribution map[string]float64) string {
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
	return ""
}

func ShuffleWithLabels(values [][]uint8, labels []float64) ([][]uint8, []float64) {
	for i := len(values) - 1; i > 0; i-- {
		j := rand.Intn(i + 1)
		values[i], values[j] = values[j], values[i]
		labels[i], labels[j] = labels[j], labels[i]
	}
	return values, labels
}

func Shuffle(values [][]uint8) [][]uint8 {
	shuffled := make([][]uint8, len(values))
	copy(shuffled, values)
	for i := len(values) - 1; i > 0; i-- {
		j := rand.Intn(i + 1)
		shuffled[i], shuffled[j] = shuffled[j], shuffled[i]
	}
	return shuffled
}
