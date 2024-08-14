package solver

import (
	"github.com/cockroachdb/apd/v3"
	"log"
)

type Rounder struct {
	RoundedMode    bool
	PrecisionLimit uint32
	ScaledWeights  []*apd.BigInt
	ScaledTarget   *apd.BigInt
}

func newRounder(pl uint32) *Rounder {
	return &Rounder{
		RoundedMode:    false,
		PrecisionLimit: pl,
		ScaledWeights:  make([]*apd.BigInt, 0),
		ScaledTarget:   new(apd.BigInt),
	}
}

func (dp *DP) getTargetAndWeightsAsInts() ([]int64, int64, int64) {
	var err error
	target := new(apd.Decimal)
	target.Set(dp.target)
	if len(dp.known) > 0 {
		for i := range dp.p.Loci {
			if _, ok := dp.known[i]; !ok {
				continue
			}
			if (dp.known[i] == 0 && dp.p.EffectAlleles[i] == 0) || (dp.known[i] == 2 && dp.p.EffectAlleles[i] == 1) {
				_, err = dp.p.Context.Sub(target, target, dp.p.Weights[i])
				_, err = dp.p.Context.Sub(target, target, dp.p.Weights[i])
			}
			if dp.known[i] == 1 {
				_, err = dp.p.Context.Sub(target, target, dp.p.Weights[i])
			}
			if err != nil {
				log.Fatalf("Failed to subtract known locus from target: %s", dp.p.Weights[i].String())
			}
		}
	}

	var roundingError int64 = 0
	multiplier := apd.New(1, int32(dp.p.WeightPrecision))
	if dp.p.WeightPrecision > dp.rounder.PrecisionLimit {
		//roundingError = int64(p.NumVariants) * 5 / 4
		roundingError = int64(dp.p.NumVariants - len(dp.known))
		dp.rounder.RoundedMode = true
		dp.rounder.ScaledWeights = scaleWeights(dp.p.Context, dp.p.Weights, multiplier)
		//for i := range dp.rounder.ScaledWeights {
		//	fmt.Printf("%s ", dp.rounder.ScaledWeights[i].String())
		//}
		dp.rounder.ScaledTarget = DecimalToBigInt(dp.p.Context, target, multiplier)
		multiplier.SetFinite(1, int32(dp.rounder.PrecisionLimit))
	}
	//for i := range dp.p.Weights {
	//	fmt.Printf("%s ", dp.p.Weights[i].String())
	//}
	//fmt.Printf("\n")
	//fmt.Printf("Target: %s\n", target.String())
	weights := DecimalsToInts(dp.p.Context, dp.p.Weights, multiplier)
	tmp := new(apd.Decimal)
	_, err = dp.p.Context.Mul(tmp, target, multiplier)
	if err != nil {
		log.Fatalf("Failed to multiply ogTarget and multiplier: %s", target.String())
	}
	_, err = dp.p.Context.RoundToIntegralValue(tmp, tmp)
	if err != nil {
		log.Fatalf("Failed to round ogTarget: %s", tmp.String())
	}
	intTarget, err := tmp.Int64()
	if err != nil {
		log.Fatalf("Failed to convert ogTarget decimal to int64: %s", tmp.String())
	}
	return weights, intTarget, roundingError
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
	var err error
	tmp := new(apd.Decimal)
	_, err = ctx.Mul(tmp, d, multiplier)
	if err != nil {
		log.Fatalf("Failed to multiply decimal by multiplier: %s", d.String())
	}
	roundingCtx := apd.BaseContext.WithPrecision(ctx.Precision)
	roundingCtx.Rounding = apd.RoundHalfUp
	_, err = roundingCtx.RoundToIntegralValue(tmp, tmp)
	if err != nil {
		log.Fatalf("Failed to round decimal: %s", tmp.String())
	}
	decStr := tmp.Text('f')
	b := apd.NewBigInt(0)
	_, success := b.SetString(decStr, 10)
	if !success {
		log.Fatalf("Failed to convert decimal string %s to big.Int %s", decStr, b.String())
	}
	return b
}

func scaleWeights(ctx *apd.Context, weights []*apd.Decimal, multiplier *apd.Decimal) []*apd.BigInt {
	scaled := make([]*apd.BigInt, len(weights))
	for i := range weights {
		scaled[i] = DecimalToBigInt(ctx, weights[i], multiplier)
	}
	return scaled
}
