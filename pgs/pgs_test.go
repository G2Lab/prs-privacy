package pgs

import (
	"fmt"
	"path"
	"testing"

	"github.com/nikirill/prs/params"
)

func TestPGSLoad(t *testing.T) {
	p := NewPGS()
	err := p.LoadCatalogFile(path.Join(params.LocalDataFolder, "PGS000073_hmPOS_GRCh37.txt"))
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
	fmt.Println("Variants Number:", p.NumVariants)
	fmt.Println("Field Names:", p.Fieldnames)
	fmt.Println("Variants:")
	for _, variant := range p.Variants {
		fmt.Println("ID:", variant.GetLocus())
		fmt.Println("Hm Chr:", variant.GetHmChr())
		fmt.Println("Hm Pos:", variant.GetHmPos())
		fmt.Println("----")
	}
}
