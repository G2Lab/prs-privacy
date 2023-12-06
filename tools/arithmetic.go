package tools

import (
	"github.com/nikirill/prs/params"
	"math"
)

// Mod function
func Mod[T params.Integer](n, modulus T) T {
	result := n % modulus
	if result < 0 {
		result += modulus
	}
	return result
}

// Modular addition
func AddMod[T params.Integer](a, b, modulus T) T {
	return Mod(a+b, modulus)
}

// Modular subtraction
func SubMod[T params.Integer](a, b, modulus T) T {
	return Mod(a-b, modulus)
}

func FindNextBiggerPrime(n int64) int64 {
	for !IsPrime(n) {
		n++
	}
	return n
}

func FindNextSmallerPrime(n int64) int64 {
	for !IsPrime(n) {
		n--
	}
	return n
}

func IsPrime(n int64) bool {
	n = int64(math.Abs(float64(n)))
	if n == 1 {
		return false
	}
	if n == 2 {
		return true
	}
	maxDivisor := int64(math.Sqrt(float64(n)))
	var i int64
	for i = 2; i <= maxDivisor; i++ {
		if n%i == 0 {
			return false
		}
	}
	return true
}
