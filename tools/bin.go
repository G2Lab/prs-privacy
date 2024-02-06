package tools

import "math"

func ValueToBinIdx(value float64, numBins int) int {
	idx := int(value * float64(numBins))
	if idx == numBins {
		idx--
	}
	return idx
}

func DeriveNumSpectrumBins(l int) int {
	return int(math.Ceil(math.Sqrt(float64(l)))) + 1
}
