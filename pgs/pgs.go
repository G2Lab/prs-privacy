package pgs

import (
	"bufio"
	"encoding/json"
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"os/exec"
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
	NumHplt    = 2
	MissingEAF = 0.0004
	//NumSpectrumBins = 20
)

var (
	GENOTYPES   = []uint8{0, 1}
	POPULATIONS = []string{"AFR", "AMR", "EAS", "EUR", "SAS"}
)

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

	if frequency, ok := fields["allelefrequency_effect"].(string); ok {
		if value, err := strconv.ParseFloat(frequency, 64); err == nil {
			fields["allelefrequency_effect"] = value
		} else {
			log.Printf("Error parsing frequency %s: %s", frequency, err)
		}
	}

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
	//log.Fatalf("PopulationEAF not found for variant %s", v.GetLocus())
	return 0.0
}

type PGS struct {
	PgsID           string
	TraitName       string
	TraitEFO        string
	GenomeBuild     string
	WeightType      string
	HmPOSBuild      string
	NumVariants     int
	Fieldnames      []string
	Variants        map[string]*Variant
	Loci            []string
	Weights         []*apd.Decimal
	Context         *apd.Context
	WeightPrecision uint32
	//Maf             [][]float64 // [major, minor] allele frequency from the population
	PopulationEAF map[string][][]float64
	StudyEAF      [][]float64          // [reference, effect] allele frequency from the study / catalogue file
	FreqSpec      map[string][]float64 // Allele Frequency Spectrum per population
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
					p.NumVariants = value
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
					fmt.Printf("No mapping for variant %s\n", values[0])
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
	p.StudyEAF = make([][]float64, len(p.Loci))
	for i, loc := range p.Loci {
		p.Weights[i], err = p.Variants[loc].GetWeight(p.Context)
		if err != nil {
			log.Fatalf("Variant %s, %v: %v\n", loc, err, p.Variants[loc].fields["effect_weight"])
		}
		p.StudyEAF[i] = []float64{1 - p.Variants[loc].GetEffectAlleleFrequency(), p.Variants[loc].GetEffectAlleleFrequency()}
	}
	p.WeightPrecision = maxPrecision
	fmt.Printf("Weight precision: %d digits\n", p.WeightPrecision)
	p.PopulationEAF = make(map[string][][]float64, len(POPULATIONS))
	p.FreqSpec = make(map[string][]float64, len(POPULATIONS))

	if err = scanner.Err(); err != nil {
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
	p.loadDatasetEAF()
	p.LoadFrequencySpectrumData()
}

func (p *PGS) loadDatasetEAF() {
	filename := fmt.Sprintf("%s/%s.eaf", params.DataFolder, p.PgsID)
	_, err := os.Stat(filename)
	if os.IsNotExist(err) {
		p.extractEAF(filename)
	}
	eafFile, err := os.Open(filename)
	if err != nil {
		log.Printf("Error opening MAF file: %v\n", err)
		return
	}
	defer eafFile.Close()
	decoder := json.NewDecoder(eafFile)
	err = decoder.Decode(&p.PopulationEAF)
	if err != nil {
		log.Printf("Error decoding json PopulationEAF: %v", err)
	}
}

func (p *PGS) extractEAF(filename string) {
	eaf := make(map[string][][]float64, len(POPULATIONS))
	for i := range POPULATIONS {
		eaf[POPULATIONS[i]] = make([][]float64, 0, len(p.Loci))
	}
	populationQ := "%CHROM:%POS-%" + strings.Join(POPULATIONS, "_AF\\t%") + "_AF\n"
	var freq, parsed float64
	for _, locus := range p.Loci {
		chr, pos := tools.SplitLocus(locus)
		query, args := tools.RangeQuery(populationQ, chr, pos, pos)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			log.Printf("Error executing bcftools command: %v", err)
			continue
		}
		lines := strings.Split(string(output), "\n")
		//fmt.Printf("Locus %s: %d lines\n", locus, len(lines))
		if len(lines[:len(lines)-1]) == 0 {
			log.Printf("No data for locus %s, inserting default EAF", locus)
			for i := range POPULATIONS {
				eaf[POPULATIONS[i]] = append(eaf[POPULATIONS[i]], []float64{1 - MissingEAF, MissingEAF})
			}
			continue
		}
		for _, line := range lines[:len(lines)-1] {
			fields := strings.Split(line, "-")
			if fields[0] != locus {
				fmt.Printf("Locus %s does not match %s\n", locus, fields[0])
				continue
			}
			afPerPopulation := strings.Split(fields[1], "\t")
			for i, population := range POPULATIONS {
				altAfs := strings.Split(afPerPopulation[i], ",")
				freq = 0
				for _, altAf := range altAfs {
					parsed, err = strconv.ParseFloat(altAf, 64)
					if err != nil {
						log.Printf("Error parsing %s at %s: %v", afPerPopulation[i], locus, err)
						parsed = MissingEAF
					}
					freq += parsed
				}
				switch freq {
				case 0:
					freq = MissingEAF
				case 1:
					freq = 1 - MissingEAF
				default:
				}
				if freq > 1 || freq < 0 {
					log.Printf("Allele frequency is wrong %f for %s at %s", freq, afPerPopulation[i], locus)
				}
				eaf[population] = append(eaf[population], []float64{1 - freq, freq})
			}
		}
	}
	// save the output
	file, err := os.Create(filename)
	if err != nil {
		log.Printf("Error creating PopulationEAF file: %v", err)
		return
	}
	defer file.Close()
	encoder := json.NewEncoder(file)
	err = encoder.Encode(eaf)
	if err != nil {
		log.Printf("Error encoding json PopulationEAF: %v", err)
	}
}

func (p *PGS) LoadFrequencySpectrumData() {
	filename := fmt.Sprintf("%s/%s.freqspec", params.DataFolder, p.PgsID)
	_, err := os.Stat(filename)
	if os.IsNotExist(err) {
		p.computeFrequencySpectrum(filename)
	}
	freqSpecFile, err := os.Open(filename)
	if err != nil {
		log.Printf("Error opening freqspec file: %v", err)
		return
	}
	defer freqSpecFile.Close()
	decoder := json.NewDecoder(freqSpecFile)
	err = decoder.Decode(&p.FreqSpec)
	if err != nil {
		log.Printf("Error decoding json freqspec: %v", err)
	}
}

func (p *PGS) computeFrequencySpectrum(filename string) {
	numSpectrumBins := deriveNumSpectrumBins(len(p.Loci))
	freqSpec := make(map[string][]float64, len(POPULATIONS))
	for _, population := range POPULATIONS {
		freqSpec[population] = make([]float64, numSpectrumBins)
	}
	ancestry := tools.LoadAncestry()
	ancestrySizes := make(map[string]int, len(ancestry))
	//var popNames []string
	for _, popName := range ancestry {
		if strings.Contains(popName, ",") {
			popName = strings.Split(popName, ",")[0]
		}
		if _, exists := ancestrySizes[popName]; !exists {
			//// If the individual has multiple ancestries, we count all
			//if strings.Contains(popName, ",") {
			//	popNames = strings.Split(popName, ",")
			//	for _, p := range popNames {
			//		if _, exists = ancestrySizes[p]; !exists {
			//			ancestrySizes[p] = 0
			//		}
			//		ancestrySizes[p]++
			//	}
			//	continue
			//}
			ancestrySizes[popName] = 0
		}
		ancestrySizes[popName]++
	}
	var digit int
	var normalized, sample, ancs string
	for k, locus := range p.Loci {
		chr, pos := tools.SplitLocus(locus)
		query, args := tools.RangeGenotypesQuery(chr, pos, pos)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			log.Printf("Error executing bcftools command: %v", err)
			continue
		}
		lines := strings.Split(string(output), "\n")
		if len(lines[:len(lines)-1]) == 0 {
			log.Printf("No data for locus %s, inserting all as reference + one effect", locus)
			for _, ppl := range POPULATIONS {
				// Assume that all the samples have the reference allele, and one "ghost" sample has the effect allele
				freqSpec[ppl][valueToBinIdx(p.PopulationEAF[ppl][k][0], numSpectrumBins)] += float64(ancestrySizes[ppl])
				freqSpec[ppl][valueToBinIdx(p.PopulationEAF[ppl][k][1], numSpectrumBins)] += 1
			}
			continue
		}
		for _, line := range lines[:len(lines)-1] {
			fields := strings.Split(line, "-")
			if fields[0] != locus {
				fmt.Printf("Locus %s does not match %s\n", locus, fields[0])
				continue
			}
			samples := strings.Split(fields[1], "\t")
			for _, sample = range samples[:len(samples)-1] {
				gtp := strings.Split(sample, "=")
				ancs = ancestry[gtp[0]]
				if strings.Contains(ancs, ",") {
					ancs = strings.Split(ancs, ",")[0]
				}
				alleles := strings.Split(gtp[1], "|")
				for _, allele := range alleles {
					normalized, err = tools.NormalizeAllele(allele)
					if err != nil {
						if normalized != "." {
							log.Printf("%v: %s, at %s\n", err, sample, locus)
						}
						continue
					}
					digit, err = strconv.Atoi(normalized)
					if err != nil {
						log.Printf("Error converting %s to int: %v", normalized, err)
						continue
					}
					freqSpec[ancs][valueToBinIdx(p.PopulationEAF[ancs][k][digit], numSpectrumBins)]++
				}
			}
		}
	}
	// Normalize the frequency spectrum
	for _, population := range POPULATIONS {
		for i := 0; i < numSpectrumBins; i++ {
			freqSpec[population][i] /= float64(ancestrySizes[population])
		}
	}
	// save the output
	file, err := os.Create(filename)
	if err != nil {
		log.Printf("Error creating freqspec file: %v", err)
		return
	}
	defer file.Close()
	encoder := json.NewEncoder(file)
	err = encoder.Encode(freqSpec)
	if err != nil {
		log.Printf("Error encoding json freqspec: %v", err)
	}
}

func CalculateAlleleFrequency(sequence []uint8, af [][]float64) []float64 {
	numSpectrumBins := deriveNumSpectrumBins(len(sequence) / 2)
	alfreq := make([]float64, numSpectrumBins)
	for i := 0; i < len(sequence); i += 2 {
		for j := 0; j < NumHplt; j++ {
			alfreq[valueToBinIdx(af[i/NumHplt][sequence[i+j]], numSpectrumBins)]++
		}
	}
	return alfreq
}

func valueToBinIdx(value float64, numBins int) int {
	idx := int(value * float64(numBins))
	if idx == numBins {
		idx--
	}
	return idx
}

func deriveNumSpectrumBins(l int) int {
	return int(math.Ceil(math.Sqrt(float64(l)))) + 1
}

//func (p *PGS) loadMAFOld() {
//	filename := fmt.Sprintf("%s/%s.maf", params.DataFolder, p.PgsID)
//	_, err := os.Stat(filename)
//	if os.IsNotExist(err) {
//		p.assembleMAF(filename)
//	}
//	mafFile, err := os.Open(filename)
//	if err != nil {
//		log.Printf("Error opening MAF file: %v\n", err)
//		return
//	}
//	defer mafFile.Close()
//	reader := csv.NewReader(mafFile)
//	// Read header
//	_, err = reader.Read()
//	if err != nil {
//		log.Printf("Error reading header: %v", err)
//		return
//	}
//	p.Maf = make([][]float64, 0, len(p.Loci))
//	for {
//		row, err := reader.Read()
//		if err != nil {
//			break // Reached end of file or encountered an error
//		}
//		effectLikelihood, err := strconv.ParseFloat(row[2], 64)
//		if err != nil {
//			log.Printf("Error converting %s to float: %v", row[2], err)
//			continue
//		}
//		if effectLikelihood == 0 {
//			effectLikelihood = MissingEAF
//		}
//		p.Maf = append(p.Maf, []float64{1 - effectLikelihood, effectLikelihood})
//	}
//	//for i := 0; i < len(p.Loci); i++ {
//	//	p.Maf[i] = p.StudyEAF[i]
//	//}
//}
//
//// Retrieve Minor-Allele Frequency for each variant from the database
//func (p *PGS) assembleMAF(filename string) {
//	file, err := os.Create(filename)
//	if err != nil {
//		log.Printf("Error creating priors file: %v", err)
//		return
//	}
//	defer file.Close()
//	writer := csv.NewWriter(file)
//	writer.Write([]string{"chromosome", "position", "MAF"})
//	for _, locus := range p.Loci {
//		chr, pos := tools.SplitLocus(locus)
//		f, err := os.Open(fmt.Sprintf("data/prior/chr%s.csv", chr))
//		if err != nil {
//			log.Printf("Error opening file chr%s.csv: %s", chr, err)
//			continue
//		}
//		reader := csv.NewReader(f)
//		// Read header
//		_, err = reader.Read()
//		if err != nil {
//			fmt.Printf("Error reading MAF file chr%s.csv: %v", chr, err)
//			continue
//		}
//		found := false
//		for {
//			row, err := reader.Read()
//			if err != nil {
//				break // Reached end of file or encountered an error
//			}
//			if row[0] != pos {
//				continue
//			}
//			writer.Write([]string{chr, pos, row[1]})
//			found = true
//			break
//		}
//		if !found {
//			writer.Write([]string{chr, pos, strconv.FormatFloat(MissingEAF, 'f', 6, 64)})
//			log.Printf("No MAF found for locus %s", locus)
//		}
//		writer.Flush()
//		f.Close()
//	}
//}
