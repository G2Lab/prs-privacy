package tools

import (
	"math/rand"
	"time"
)

func SampleFromSlice(distribution []float64) int {
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
