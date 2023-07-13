package tools

import (
	"strconv"
	"strings"
)

func ParseLocus(locus string) (int, int, error) {
	chr, err := strconv.Atoi(strings.Split(locus, ":")[0])
	if err != nil {
		return 0, 0, err
	}
	pos, err := strconv.Atoi(strings.Split(locus, ":")[1])
	if err != nil {
		return 0, 0, err
	}
	return chr, pos, nil
}

func SplitLocus(locus string) (string, string) {
	parts := strings.Split(locus, ":")
	return parts[0], parts[1]
}
