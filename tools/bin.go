package tools

import "math"

func ValueToBinIdx(value float64, bounds []float64) int {
	low := 0
	high := len(bounds) - 1

	for low <= high {
		mid := (low + high) / 2
		if value <= bounds[mid] {
			high = mid - 1
		} else {
			low = mid + 1
		}
	}

	return low
}

func DeriveNumSpectrumBins(l int) int {
	return int(math.Ceil(math.Sqrt(float64(l))))
}
