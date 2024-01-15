package pgs

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"github.com/cockroachdb/apd/v3"
	"github.com/nikirill/prs/params"
	"github.com/nikirill/prs/tools"
)

var ALL_FIELDS = []string{
	"rsID",
	"chr_name",
	"chr_position",
	"effect_allele",
	"other_allele",
	"effect_weight",
	"allelefrequency_effect",
	"locus_name",
	"variant_description",
	"OR",
	"hm_source",
	"hm_rsID",
	"hm_chr",
	"hm_pos",
	"hm_inferOtherAllele",
}

const (
	NumHplt                       = 2
	MissingEffectAlleleLikelihood = 0.0004
)

var GENOTYPES = []uint8{0, 1}

type Variant struct {
	fields map[string]interface{}
}

func NewVariant(fields map[string]interface{}) *Variant {
	v := &Variant{
		fields: make(map[string]interface{}),
	}

	//if weight, ok := fields["effect_weight"].(string); ok {
	//	value := new(big.Rat)
	//	if _, ok := value.SetString(weight); ok {
	//		fields["effect_weight"] = value
	//	} else {
	//		log.Printf("Error parsing weight %s: %s", weight, value.RatString())
	//	}
	//}

	//if frequency, ok := fields["allelefrequency_effect"].(string); ok {
	//	if value, err := strconv.ParseFloat(frequency, 64); err == nil {
	//		fields["allelefrequency_effect"] = value
	//	} else {
	//		log.Printf("Error parsing frequency %s: %s", frequency, err)
	//	}
	//}

	for _, field := range ALL_FIELDS {
		if _, exists := fields[field]; !exists {
			fields[field] = nil
		}
	}

	v.fields = fields
	return v
}

func (v *Variant) GetHmChr() string {
	return v.fields["hm_chr"].(string)
}

func (v *Variant) GetHmPos() string {
	return v.fields["hm_pos"].(string)
}

func (v *Variant) GetLocus() string {
	return fmt.Sprintf("%s:%s", v.GetHmChr(), v.GetHmPos())
}

func (v *Variant) GetWeight(ctx *apd.Context) (*apd.Decimal, error) {
	if value, ok := v.fields["effect_weight"].(string); ok {
		weight, _, err := ctx.NewFromString(value)
		if err != nil {
			return nil, err
		}
		return weight, nil
	}
	return nil, errors.New("error parsing weight")
}

func (v *Variant) GetEffectAlleleFrequency() float64 {
	if eaf, ok := v.fields["allelefrequency_effect"].(float64); ok {
		return eaf
	}
	//log.Fatalf("EAF not found for variant %s", v.GetLocus())
	return 0.0
}

type PGS struct {
	PgsID           string
	TraitName       string
	TraitEFO        string
	GenomeBuild     string
	WeightType      string
	HmPOSBuild      string
	VariantCount    int
	Fieldnames      []string
	Variants        map[string]*Variant
	Loci            []string
	Weights         []*apd.Decimal
	Context         *apd.Context
	WeightPrecision uint32
	Maf             [][]float64 // [major, minor] allele frequency from the population
	//Eaf             [][]float64 // [other, effect] allele frequency from the study / catalogue file
}

func NewPGS() *PGS {
	return &PGS{
		Variants: make(map[string]*Variant),
	}
}

func (p *PGS) LoadCatalogFile(inputFile string) error {
	file, err := os.Open(inputFile)
	if err != nil {
		return err
	}
	defer file.Close()

	headerInProgress := true
	maxPrecision := uint32(0)
	scanner := bufio.NewScanner(file)
scannerLoop:
	for scanner.Scan() {
		line := scanner.Text()

		// If it is the header
		if strings.HasPrefix(line, "#") {
			fields := strings.SplitN(line[1:], "=", 2)
			switch strings.ToLower(fields[0]) {
			case "pgs_id":
				p.PgsID = fields[1]
			case "trait_mapped":
				p.TraitName = fields[1]
			case "trait_efo":
				p.TraitEFO = fields[1]
			case "genome_build":
				p.GenomeBuild = fields[1]
			case "weight_type":
				p.WeightType = fields[1]
			case "hmpos_build":
				p.HmPOSBuild = fields[1]
			case "variants_number":
				if value, err := strconv.Atoi(fields[1]); err == nil {
					p.VariantCount = value
				} else {
					log.Printf("Error parsing variants number %s: %s", fields[1], err)
				}
			}
			continue scannerLoop
		}

		// Specified variant fields
		if headerInProgress {
			p.Fieldnames = strings.Split(line, "\t")
			headerInProgress = false
			continue scannerLoop
		}

		fields := make(map[string]interface{})
		values := strings.Split(line, "\t")
		for i, value := range values {
			if p.Fieldnames[i] == "hm_chr" {
				if strings.Contains(value, "_") {
					value = strings.Split(value, "_")[0]
					if value == "Un" {
						for j := 0; j < len(p.Fieldnames); j++ {
							if p.Fieldnames[j] == "chr_name" {
								value = values[j]
								break
							}
						}
					}
				}
				// If there is no mapping, we skip the variant
				if len(value) == 0 {
					continue scannerLoop
				}
			}
			fields[p.Fieldnames[i]] = value
		}
		if maxPrecision < getPrecision(fields["effect_weight"].(string)) {
			maxPrecision = getPrecision(fields["effect_weight"].(string))
		}
		variant := NewVariant(fields)
		p.Variants[variant.GetLocus()] = variant
	}
	//p.Loci = p.GetUnSortedVariantLoci()
	p.Loci, err = p.GetSortedVariantLoci()
	if err != nil {
		return err
	}

	p.Context = &apd.Context{
		Precision:   maxPrecision + 2,
		Rounding:    apd.RoundFloor,
		MaxExponent: int32(maxPrecision) + 3,
		MinExponent: -int32(maxPrecision) - 3,
		Traps:       apd.DefaultTraps,
	}

	p.Weights = make([]*apd.Decimal, len(p.Loci))
	//p.Eaf = make([][]float64, len(p.Loci))
	for i, loc := range p.Loci {
		p.Weights[i], err = p.Variants[loc].GetWeight(p.Context)
		if err != nil {
			log.Fatalf("Variant %s, %v: %v\n", loc, err, p.Variants[loc].fields["effect_weight"])
		}
		//p.Eaf[i] = []float64{1 - p.Variants[loc].GetEffectAlleleFrequency(), p.Variants[loc].GetEffectAlleleFrequency()}
	}
	p.WeightPrecision = maxPrecision
	fmt.Printf("Weight precision: %d digits\n", p.WeightPrecision)

	if err := scanner.Err(); err != nil {
		return err
	}

	return nil
}

func getPrecision(value string) uint32 {
	if strings.Contains(value, "e") {
		exp := strings.Split(value, "e")[0]
		expLen := len(exp)
		if strings.Contains(exp, ".") {
			expLen = len(strings.Split(exp, ".")[1])
		}
		mntLen, err := strconv.Atoi(strings.Split(value, "e")[1][1:])
		if err != nil {
			log.Printf("Error parsing mantissa %s: %v", value, err)
			return 0
		}
		return uint32(expLen + mntLen)
	}
	if !strings.Contains(value, ".") {
		return 0
	}
	return uint32(len(strings.Split(value, ".")[1]))
}

func (p *PGS) GetSortedVariantLoci() ([]string, error) {
	sortedLoc := make([]string, 0, len(p.Variants))
	for locus := range p.Variants {
		sortedLoc = append(sortedLoc, locus)
	}
	for i := 0; i < len(sortedLoc)-1; i++ {
		minIndex := i
		minChr, minPos, err := tools.ParseLocus(sortedLoc[minIndex])
		if err != nil {
			log.Printf("Error parsing initial locus %s: %v", sortedLoc[minIndex], err)
			return nil, err
		}
		for j := i + 1; j < len(sortedLoc); j++ {
			chr, pos, err := tools.ParseLocus(sortedLoc[j])
			if err != nil {
				log.Printf("Error parsing locus %s: %v", sortedLoc[j], err)
				return nil, err
			}
			if chr < minChr || (chr == minChr && pos < minPos) {
				minIndex = j
				minChr = chr
				minPos = pos
			}
		}
		if minIndex != i {
			sortedLoc[i], sortedLoc[minIndex] = sortedLoc[minIndex], sortedLoc[i]
		}
	}
	return sortedLoc, nil
}

func (p *PGS) GetUnSortedVariantLoci() []string {
	loci := make([]string, 0, len(p.Variants))
	for locus := range p.Variants {
		loci = append(loci, locus)
	}
	return loci
}

func (p *PGS) LoadStats() {
	p.loadMAF()
}

func (p *PGS) loadMAF() {
	filename := fmt.Sprintf("%s/%s.maf", params.DataFolder, p.PgsID)
	_, err := os.Stat(filename)
	if os.IsNotExist(err) {
		p.assembleMAF(filename)
	}
	mafFile, err := os.Open(filename)
	if err != nil {
		log.Printf("Error opening MAF file: %v\n", err)
		return
	}
	defer mafFile.Close()
	reader := csv.NewReader(mafFile)
	// Read header
	_, err = reader.Read()
	if err != nil {
		log.Printf("Error reading header: %v", err)
		return
	}
	p.Maf = make([][]float64, 0, len(p.Loci))
	for {
		row, err := reader.Read()
		if err != nil {
			break // Reached end of file or encountered an error
		}
		effectLikelihood, err := strconv.ParseFloat(row[2], 64)
		if err != nil {
			log.Printf("Error converting %s to float: %v", row[2], err)
			continue
		}
		if effectLikelihood == 0 {
			effectLikelihood = MissingEffectAlleleLikelihood
		}
		p.Maf = append(p.Maf, []float64{1 - effectLikelihood, effectLikelihood})
	}
}

// Retrieve Minor-Allele Frequency for each variant from the database
func (p *PGS) assembleMAF(filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Printf("Error creating priors file: %v", err)
		return
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	writer.Write([]string{"chromosome", "position", "MAF"})
	for _, locus := range p.Loci {
		chr, pos := tools.SplitLocus(locus)
		f, err := os.Open(fmt.Sprintf("data/prior/chr%s.csv", chr))
		if err != nil {
			log.Printf("Error opening file chr%s.csv: %s", chr, err)
			continue
		}
		reader := csv.NewReader(f)
		// Read header
		_, err = reader.Read()
		if err != nil {
			fmt.Printf("Error reading MAF file chr%s.csv: %v", chr, err)
			continue
		}
		found := false
		for {
			row, err := reader.Read()
			if err != nil {
				break // Reached end of file or encountered an error
			}
			if row[0] != pos {
				continue
			}
			writer.Write([]string{chr, pos, row[1]})
			found = true
			break
		}
		if !found {
			writer.Write([]string{chr, pos, "0.00000"})
			log.Printf("No MAF found for locus %s", locus)
		}
		writer.Flush()
		f.Close()
	}
}

// ShuffleIndicesByLikelihood shuffles index order by the likelihood of alternative allele.
// It ensures that if there are several solutions away from the current delta, we tend to pick the one with the highest
// likelihood, but it is still randomized.
func (p *PGS) ShuffleIndicesByLikelihood(original []uint8) []int {
	var lkl float64
	weightedIndices := make([]int, 0)
	for i, maf := range p.Maf {
		for j := 0; j < NumHplt; j++ {
			weightedIndices = append(weightedIndices, NumHplt*i+j)
			lkl = 1 - maf[original[NumHplt*i+j]]
			for k := 0; k < int(lkl*100); k++ {
				weightedIndices = append(weightedIndices, NumHplt*i+j)
			}
		}
	}
	rand.Shuffle(len(weightedIndices), func(i, j int) {
		weightedIndices[i], weightedIndices[j] = weightedIndices[j], weightedIndices[i]
	})
	indices := make([]int, 0, len(original))
	unique := make(map[int]struct{})
	for _, v := range weightedIndices {
		if _, exists := unique[v]; !exists {
			indices = append(indices, v)
			unique[v] = struct{}{}
		}
	}
	return indices
}

// SampleFromPopulation samples a candidate according to the MAF
func (p *PGS) SampleFromPopulation() ([]uint8, error) {
	sample := make([]uint8, p.VariantCount*NumHplt)
	// Initial sample based on individual priors
	for i := range p.Loci {
		for j := 0; j < NumHplt; j++ {
			maf := p.Maf[i]
			ind := tools.SampleFromDistribution(maf)
			sample[i*NumHplt+j] = GENOTYPES[ind]
			if sample[i*NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) SampleSegmentFromPopulation(start, end int) ([]uint8, error) {
	sample := make([]uint8, (end-start)*NumHplt)
	// Initial sample based on individual priors
	for i := start; i < end; i++ {
		for j := 0; j < NumHplt; j++ {
			maf := p.Maf[i]
			ind := tools.SampleFromDistribution(maf)
			sample[(i-start)*NumHplt+j] = GENOTYPES[ind]
			if sample[(i-start)*NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) SampleUniform() ([]uint8, error) {
	sample := make([]uint8, p.VariantCount*NumHplt)
	// Initial sample based on individual priors
	for i := range p.Loci {
		for j := 0; j < NumHplt; j++ {
			sample[i*NumHplt+j] = tools.SampleUniform(GENOTYPES)
			if sample[i*NumHplt+j] == 255 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) CalculateSequenceLikelihood(sequence []uint8) float64 {
	//likelihood := 1.0
	likelihood := 0.0
	for i := range p.Loci {
		for j := 0; j < NumHplt; j++ {
			likelihood += math.Log(p.Maf[i][sequence[i*NumHplt+j]])
			//likelihood += math.Log(p.Eaf[i][sequence[i*NumHplt+j]])
		}
	}
	return likelihood
}

func (p *PGS) SnpLikelihood(sequence []uint8, i int) float64 {
	likelihood := 0.0
	for j := 0; j < NumHplt; j++ {
		likelihood += math.Log(p.Maf[i][sequence[i*NumHplt+j]])
	}
	return likelihood
}

func (p *PGS) AllMajorSample() []uint8 {
	sample := make([]uint8, 2*len(p.Weights))
	for i := 0; i < len(p.Weights); i++ {
		if p.Maf[i][0] > 0.5 {
			sample[2*i] = 0
			sample[2*i+1] = 0
		} else {
			sample[2*i] = 1
			sample[2*i+1] = 1
		}
	}
	return sample
}

func (p *PGS) AllMinorSample() []uint8 {
	sample := make([]uint8, 2*len(p.Weights))
	for i := 0; i < len(p.Weights); i++ {
		if p.Maf[i][0] > 0.5 {
			sample[2*i] = 1
			sample[2*i+1] = 1
		} else {
			sample[2*i] = 0
			sample[2*i+1] = 0
		}
	}
	return sample
}
