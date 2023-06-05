package main

import (
	"fmt"
	"testing"
)

func TestPGSLoad(t *testing.T) {
	p := NewPGS()
	err := p.Load("PGS000073_hmPOS_GRCh38.txt")
	if err != nil {
		fmt.Println("Error:", err)
		return
	}

	// Access the loaded data as needed
	fmt.Println("PGS ID:", p.PgsID)
	fmt.Println("Trait Name:", p.TraitName)
	fmt.Println("Trait EFO:", p.TraitEFO)
	fmt.Println("Genome Build:", p.GenomeBuild)
	fmt.Println("Weight Type:", p.WeightType)
	fmt.Println("HmPOS Build:", p.HmPOSBuild)
	fmt.Println("Variants Number:", p.VariantsNumber)
	fmt.Println("Field Names:", p.Fieldnames)
	fmt.Println("Variants:")
	for _, variant := range p.Variants {
		fmt.Println("ID:", variant.GetID())
		fmt.Println("Weight:", variant.GetWeight())
		fmt.Println("Hm Chr:", variant.GetHmChr())
		fmt.Println("Hm Pos:", variant.GetHmPos())
		fmt.Println("----")
	}
}
