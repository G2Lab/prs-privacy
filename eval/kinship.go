package main

func kingHomo(ind1, ind2 map[string]uint8, af map[string][]float32) float32 {
	var NAAaa, NAaAa, NAa1, NAa2, H float32 = 0.0, 0.0, 0.0, 0.0, 0.0
	for locus := range ind1 {
		if _, ok := ind2[locus]; !ok {
			continue
		}
		if ind1[locus] == 1 && ind2[locus] == 1 {
			NAaAa++
		}
		if (ind1[locus] == 2 && ind2[locus] == 0) || (ind1[locus] == 0 && ind2[locus] == 2) {
			NAAaa++
		}
		if ind1[locus] == 1 {
			NAa1++
		}
		if ind2[locus] == 1 {
			NAa2++
		}
		H += 2 * af[locus][0] * af[locus][1]
	}
	return 1/2 + (NAaAa-2*NAAaa)/(2*H) - (NAa1-NAa2)/(4*H)
}

func kingRobust(ind1, ind2 map[string]uint8) float32 {
	var NAaAa, NAAaa, NAa1, NAa2 float32 = 0.0, 0.0, 0.0, 0.0
	for locus := range ind1 {
		if _, ok := ind2[locus]; !ok {
			continue
		}
		if ind1[locus] == 1 && ind2[locus] == 1 {
			NAaAa++
		}
		if (ind1[locus] == 2 && ind2[locus] == 0) || (ind1[locus] == 0 && ind2[locus] == 2) {
			NAAaa++
		}
		if ind1[locus] == 1 {
			NAa1++
		}
		if ind2[locus] == 1 {
			NAa2++
		}
	}
	return (NAaAa - 2*NAAaa) / (NAa1 + NAa2)
}

func kingRobustPair(ind1, ind2 map[string]uint8) []int {
	var NAaAa, NAAaa, NAa1, NAa2 int = 0, 0, 0, 0
	for locus := range ind1 {
		if _, ok := ind2[locus]; !ok {
			continue
		}
		if ind1[locus] == 1 && ind2[locus] == 1 {
			NAaAa++
		}
		if (ind1[locus] == 2 && ind2[locus] == 0) || (ind1[locus] == 0 && ind2[locus] == 2) {
			NAAaa++
		}
		if ind1[locus] == 1 {
			NAa1++
		}
		if ind2[locus] == 1 {
			NAa2++
		}
	}
	return []int{NAaAa - 2*NAAaa, NAa1 + NAa2}
}

//func kingRobustBetween(ind1, ind2 map[string]uint8) float32 {
//	var NAaAa, NAAaa, NAa1, NAa2, NAaMin float32 = 0.0, 0.0, 0.0, 0.0, 0.0
//	for locus := range ind1 {
//		if _, ok := ind2[locus]; !ok {
//			continue
//		}
//		if ind1[locus] == 1 && ind2[locus] == 1 {
//			NAaAa++
//		}
//		if (ind1[locus] == 2 && ind2[locus] == 0) || (ind1[locus] == 0 && ind2[locus] == 2) {
//			NAAaa++
//		}
//		if ind1[locus] == 1 {
//			NAa1++
//		}
//		if ind2[locus] == 1 {
//			NAa2++
//		}
//	}
//	if NAa1 > NAa2 {
//		NAaMin = NAa2
//	} else {
//		NAaMin = NAa1
//	}
//	return 1/2 + (2*NAaAa-4*NAAaa-NAa1-NAa2)/(4*NAaMin)
//}
