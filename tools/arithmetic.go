package tools

import (
	"math"

	"github.com/ericlagergren/decimal"
	"github.com/nikirill/prs/params"
)

func BigMod(ctx decimal.Context, n, modulus *decimal.Big) *decimal.Big {
	rem := decimal.WithContext(ctx)
	ctx.Rem(rem, n, modulus)
	if rem.Sign() == -1 {
		rem.Add(rem, modulus)
	}
	return rem
}

func BigAddMod(ctx decimal.Context, a, b, modulus *decimal.Big) *decimal.Big {
	sum := decimal.WithContext(ctx)
	ctx.Add(sum, a, b)
	return BigMod(ctx, sum, modulus)
}

func BigSubMod(ctx decimal.Context, a, b, modulus *decimal.Big) *decimal.Big {
	diff := decimal.WithContext(ctx)
	ctx.Sub(diff, a, b)
	return BigMod(ctx, diff, modulus)
}

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

func FindNextBiggerPrime(n uint64) uint64 {
	for !IsPrime(n) {
		n++
	}
	return n
}

func FindNextSmallerPrime(n uint64) uint64 {
	for !IsPrime(n) {
		n--
	}
	return n
}

func IsPrime(n uint64) bool {
	n = uint64(math.Abs(float64(n)))
	if n == 2 {
		return true
	}
	if n == 1 || n%2 == 0 {
		return false
	}
	maxDivisor := uint64(math.Sqrt(float64(n)))
	var i uint64
	for i = 3; i <= maxDivisor; i += 2 {
		if n%i == 0 {
			return false
		}
	}
	return true
}
