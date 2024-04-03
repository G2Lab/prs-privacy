package solver

import (
	"log"
	"math"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/pgs"
)

type Rounder struct {
	RoundedMode   bool
	ScaledWeights []*apd.BigInt
	ScaledTarget  *apd.BigInt
}

func newRounder() *Rounder {
	return &Rounder{
		RoundedMode:   false,
		ScaledWeights: make([]*apd.BigInt, 0),
		ScaledTarget:  new(apd.BigInt),
	}
}

func getTargetAndWeightsAsInts(p *pgs.PGS, target *apd.Decimal, rdr *Rounder) ([]int64, int64, int64) {
	var err error
	var roundingError int64 = 0
	multiplier := apd.New(1, int32(p.WeightPrecision))
	if p.WeightPrecision > params.PrecisionsLimit {
		//roundingError = int64(p.NumVariants) * 5 / 4
		roundingError = int64(p.NumVariants)
		rdr.RoundedMode = true
		rdr.ScaledWeights = scaleWeights(p.Context, p.Weights, multiplier)
		rdr.ScaledTarget = DecimalToBigInt(p.Context, target, multiplier)
		multiplier.SetFinite(1, params.PrecisionsLimit)
		//fmt.Printf("Scaled ogTarget: %s\n", rdr.ScaledTarget.String())
		//fmt.Printf("Scaled weights: %v\n", rdr.ScaledWeights)
	}
	weights := DecimalsToInts(p.Context, p.Weights, multiplier)
	tmp := new(apd.Decimal)
	_, err = p.Context.Mul(tmp, target, multiplier)
	if err != nil {
		log.Fatalf("Failed to multiply ogTarget and multiplier: %s", target.String())
	}
	_, err = p.Context.RoundToIntegralValue(tmp, tmp)
	if err != nil {
		log.Fatalf("Failed to round ogTarget: %s", tmp.String())
	}
	intTarget, err := tmp.Int64()
	if err != nil {
		log.Fatalf("Failed to convert ogTarget decimal to int64: %s", tmp.String())
	}
	return weights, intTarget, roundingError
}

func getTargetAndWeightsAsFloats(p *pgs.PGS, target *apd.Decimal, rdr *Rounder) ([]float64, float64, float64) {
	var err error
	var roundingError float64 = 0
	precision := int(p.WeightPrecision)
	if p.WeightPrecision > params.PrecisionsLimit {
		precision = params.PrecisionsLimit
		roundingError = math.Pow10(-precision) * float64(p.NumVariants)
		multiplier := apd.New(1, int32(p.WeightPrecision))
		rdr.RoundedMode = true
		rdr.ScaledWeights = scaleWeights(p.Context, p.Weights, multiplier)
		sct := new(apd.Decimal)
		p.Context.Mul(sct, target, multiplier)
		rdr.ScaledTarget.SetString(sct.String(), 10)
		multiplier.SetFinite(1, params.PrecisionsLimit)
	}
	scale := math.Pow10(precision)
	weights := DecimalsToFloats(p.Weights, scale)
	tf, err := target.Float64()
	if err != nil {
		log.Fatalf("Failed to convert ogTarget decimal to float64: %s", target.String())
	}
	tf = float64(int64(tf*scale)) / scale
	return weights, tf, roundingError
}

func DecimalsToInts(ctx *apd.Context, decimals []*apd.Decimal, multiplier *apd.Decimal) []int64 {
	ints := make([]int64, len(decimals))
	tmp := new(apd.Decimal)
	var err error
	for i, b := range decimals {
		ctx.Mul(tmp, b, multiplier)
		_, err = ctx.RoundToIntegralValue(tmp, tmp)
		if err != nil {
			log.Fatalf("Failed to round decimal: %s", b.String())
		}
		ints[i], err = tmp.Int64()
		if err != nil {
			log.Fatalf("Failed to convert decimal to int64: %s", tmp.String())
		}
	}
	return ints
}

func DecimalsToFloats(decimals []*apd.Decimal, scale float64) []float64 {
	var err error
	var flt float64
	floats := make([]float64, len(decimals))
	for i, b := range decimals {
		flt, err = b.Float64()
		if err != nil {
			log.Fatalf("Failed to convert decimal to float64: %s", b.String())
		}
		floats[i] = float64(int64(flt*scale)) / scale
	}
	return floats
}

func DecimalToBigInt(ctx *apd.Context, d *apd.Decimal, multiplier *apd.Decimal) *apd.BigInt {
	tmp := new(apd.Decimal)
	ctx.Mul(tmp, d, multiplier)
	b := new(apd.BigInt)
	b.Set(&tmp.Coeff)
	b.Mul(b, apd.NewBigInt(int64(math.Pow(10, float64(tmp.Exponent)))))
	if tmp.Negative {
		b.Neg(b)
	}
	return b
}
