package main

import (
	"fmt"
	"github.com/cockroachdb/apd/v3"
	"github.com/ericlagergren/decimal"
	"math/big"
	"testing"
)

func BenchmarkBigRatAdd(b *testing.B) {
	rat1 := new(big.Rat).SetFloat64(0.12345678901234567)
	rat2 := new(big.Rat).SetFloat64(0.98765432109876543)

	for i := 0; i < b.N; i++ {
		result := new(big.Rat).Add(rat1, rat2)
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkDecimalAdd(b *testing.B) {
	ctx := decimal.Context{Precision: 20}
	dec1 := decimal.WithContext(ctx).SetFloat64(0.12345678901234567)
	dec2 := decimal.WithContext(ctx).SetFloat64(0.98765432109876543)

	for i := 0; i < b.N; i++ {
		result := decimal.WithContext(decimal.Context{Precision: 20}).Add(dec1, dec2)
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkAPDAdd(b *testing.B) {
	c := apd.BaseContext.WithPrecision(20)
	dec1 := apd.New(12345678901234567, -17)
	dec2 := apd.New(98765432109876543, -17)

	for i := 0; i < b.N; i++ {
		result := new(apd.Decimal)
		c.Add(result, dec1, dec2)
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkBigIntAdd(b *testing.B) {
	int1 := new(apd.BigInt).SetInt64(12345678901234567)
	int2 := new(apd.BigInt).SetInt64(98765432109876543)

	for i := 0; i < b.N; i++ {
		result := new(apd.BigInt)
		result.Add(int1, int2)
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkGoBigIntAdd(b *testing.B) {
	int1 := new(big.Int).SetInt64(12345678901234567)
	int2 := new(big.Int).SetInt64(98765432109876543)

	for i := 0; i < b.N; i++ {
		result := new(big.Int).Add(int1, int2)
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkInt64Add(b *testing.B) {
	num1 := int64(12345678901234567)
	num2 := int64(98765432109876543)

	for i := 0; i < b.N; i++ {
		result := num1 + num2
		_ = result // Avoid compiler optimization
	}
}

func BenchmarkFloat64Add(b *testing.B) {
	num1 := float64(0.12345678901234567)
	num2 := float64(0.98765432109876543)

	for i := 0; i < b.N; i++ {
		result := num1 + num2
		_ = result // Avoid compiler optimization
	}
}

func main() {
	// Run the benchmarks
	benchmarkResults := testing.Benchmark(BenchmarkBigRatAdd)
	fmt.Println("Big.Rat addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkDecimalAdd)
	fmt.Println("Decimal addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkAPDAdd)
	fmt.Println("APD Decimal addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkBigIntAdd)
	fmt.Println("APD BigInt addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkGoBigIntAdd)
	fmt.Println("Go BigInt addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkInt64Add)
	fmt.Println("Int64 addition:", benchmarkResults)

	benchmarkResults = testing.Benchmark(BenchmarkFloat64Add)
	fmt.Println("Float64 addition:", benchmarkResults)
}
