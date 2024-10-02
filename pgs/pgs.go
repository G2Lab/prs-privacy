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
	"sort"
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
	Ploidy     = 2
	MissingEAF = 0.0002
)

var (
	GENOTYPES     = []uint8{0, 1}
	ANCESTRIES    = []string{"AFR", "AMR", "EAS", "EUR", "SAS"}
	ANCESTRIESALL = []string{"AFR", "AMR", "EAS", "EUR", "SAS", "ALL"}
)

type Variant struct {
	fields map[string]interface{}
}

func NewVariant(fields map[string]interface{}) *Variant {
	v := &Variant{
		fields: make(map[string]interface{}),
	}

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

type Statistics struct {
	AF map[int][]float32 // Effect Allele Frequency
	GF map[int][]float32 // Genotype Frequency
	//FreqSpectrum  []float32         // Allele Frequency Spectrum
	//FreqBinBounds []float32         // Bounds of the frequency spectrum bins
}

type PGS struct {
	PgsID           string
	TraitName       string
	TraitEFO        string
	GenomeBuild     string
	WeightType      string
	HmPOSBuild      string
	PgpID           string
	NumVariants     int
	Fieldnames      []string
	Variants        map[string]*Variant
	Loci            []string
	Weights         []*apd.Decimal
	ScoreAlleles    []string // [effect] alleles
	EffectAlleles   []uint8  // 0 if the effect allele is the reference allele, 1 if it is the alternative allele
	Context         *apd.Context
	WeightPrecision uint32
	MinPrecision    uint32
	StudyEAF        []float64 // effect allele frequency from the study / catalogue file
	PopulationStats map[string]*Statistics
	//NumSpecBins     int
	//Maf             [][]float64 // [major, minor] allele frequency from the population
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
	precisions := make([]uint32, 0)
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
			case "pgp_id":
				p.PgpID = fields[1]
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
			if p.Fieldnames[i] == "hm_pos" || p.Fieldnames[i] == "hm_chr" {
				if len(value) == 0 || value == "Unknown" {
					fmt.Printf("%s: No mapping for variant %s\n", p.PgsID, values[0])
					return errors.New("no mapping for one of the variants")
				}
			}
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
			}
			fields[p.Fieldnames[i]] = value
		}
		precisions = append(precisions, getPrecision(fields["effect_weight"].(string)))
		variant := NewVariant(fields)
		p.Variants[variant.GetLocus()] = variant
	}
	p.Loci, err = p.GetSortedVariantLoci()
	if err != nil {
		return err
	}

	sort.Slice(precisions, func(i, j int) bool {
		return precisions[i] < precisions[j]
	})
	var minPrecision uint32
	maxPrecision := precisions[len(precisions)-1]
	if len(precisions) > 2 {
		minPrecision = precisions[2]
	} else {
		minPrecision = precisions[0]
	}

	p.Context = &apd.Context{
		Precision:   maxPrecision + 3,
		Rounding:    apd.RoundFloor,
		MaxExponent: int32(maxPrecision) + 3,
		MinExponent: -int32(maxPrecision) - 3,
		Traps:       apd.DefaultTraps,
	}

	p.Weights = make([]*apd.Decimal, len(p.Loci))
	p.ScoreAlleles = make([]string, len(p.Loci))
	p.EffectAlleles = make([]uint8, len(p.Loci))
	p.StudyEAF = make([]float64, len(p.Loci))
	for i, loc := range p.Loci {
		p.Weights[i], err = p.Variants[loc].GetWeight(p.Context)
		if err != nil {
			log.Fatalf("Variant %s, %v: %v\n", loc, err, p.Variants[loc].fields["effect_weight"])
		}
		p.ScoreAlleles[i] = p.Variants[loc].fields["effect_allele"].(string)
		p.StudyEAF[i] = p.Variants[loc].GetEffectAlleleFrequency()
	}
	p.WeightPrecision = maxPrecision
	p.MinPrecision = minPrecision
	//fmt.Printf("Weight precision: %d digits\n", p.WeightPrecision)

	//p.NumSpecBins = tools.DeriveNumSpectrumBins(len(p.Loci))
	populationsAndAll := append(ANCESTRIES, "ALL")
	p.PopulationStats = make(map[string]*Statistics, len(populationsAndAll))
	for _, population := range populationsAndAll {
		p.PopulationStats[population] = &Statistics{
			AF: make(map[int][]float32, len(p.Loci)),
			GF: make(map[int][]float32, len(p.Loci)),
			//FreqSpectrum:  make([]float32, p.NumSpecBins),
			//FreqBinBounds: make([]float32, p.NumSpecBins),
		}
	}

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

func (p *PGS) LoadStats(dataset string) error {
	return p.LoadDatasetStats(dataset)
}

func (p *PGS) LoadDatasetStats(dataset string) error {
	var efalFilename, statFilename string
	switch dataset {
	case tools.GG:
		efalFilename = fmt.Sprintf("%s/%s.efal", params.LocalDataFolder, p.PgsID)
		statFilename = fmt.Sprintf("%s/%s.stat", params.LocalDataFolder, p.PgsID)
	case tools.UKB:
		efalFilename = fmt.Sprintf("%s/%s.efal", params.UKBiobankInputFolder, p.PgsID)
		statFilename = fmt.Sprintf("%s/%s.stat", params.UKBiobankInputFolder, p.PgsID)
	default:
		log.Printf("Unknown dataset: %s\n", dataset)
		return errors.New("unknown dataset")
	}
	if _, err := os.Stat(efalFilename); os.IsNotExist(err) {
		err = p.extractPopulationAlleles(dataset)
		if err != nil {
			return err
		}
		//fmt.Println("Extracted population alleles", p.EffectAlleles)
		p.savePopulationAlleles(efalFilename)
	} else {
		efalFile, err := os.Open(efalFilename)
		if err != nil {
			log.Printf("Error opening efal file: %v\n", err)
			return err
		}
		defer efalFile.Close()
		decoder := json.NewDecoder(efalFile)
		err = decoder.Decode(&(p.EffectAlleles))
		if err != nil {
			log.Printf("Error decoding effect alleles json: %v", err)
		}
		//fmt.Println("Loaded population alleles", p.EffectAlleles)
	}

	if _, err := os.Stat(statFilename); os.IsNotExist(err) {
		switch dataset {
		case tools.GG:
			p.extractEAF()
			p.extractGF()
		case tools.UKB:
			p.CalculateEAFAndGF(dataset)
		default:
			log.Printf("Unknown dataset: %s\n", dataset)
			return errors.New("unknown dataset")
		}
		//p.computeFrequencySpectrum()
		p.SaveStats(statFilename)
	} else {
		statsFile, err := os.Open(statFilename)
		if err != nil {
			log.Printf("Error opening stats file: %v\n", err)
			return err
		}
		defer statsFile.Close()
		decoder := json.NewDecoder(statsFile)
		err = decoder.Decode(&(p.PopulationStats))
		if err != nil {
			log.Printf("Error decoding population stats json: %v", err)
		}
	}
	return nil
}

func (p *PGS) extractPopulationAlleles(dataset string) error {
	refAltQ := "%CHROM:%POS-%REF\t%ALT\n"
	var alt, ref string
	var missing, alleleSet bool
	for k, locus := range p.Loci {
		missing, alleleSet = true, false
		chr, pos := tools.SplitLocus(locus)
		query, args := tools.RangeQuery(refAltQ, chr, pos, pos, dataset)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			log.Printf("Error executing bcftools command: %v", err)
			return err
		}
		lines := strings.Split(string(output), "\n")
		for j, line := range lines[:len(lines)-1] {
			fields := strings.Split(line, "-")
			if fields[0] != locus {
				continue
			}
			missing = false
			//fmt.Printf("Locus %s: %s\n", locus, fields[1])
			ref = strings.Split(fields[1], "\t")[0]
			alt = strings.Split(fields[1], "\t")[1]
			if alleleSet {
				continue
			}
			switch {
			case strings.Contains(alt, p.ScoreAlleles[k]):
				p.EffectAlleles[k] = 1
				alleleSet = true
			case strings.Contains(ref, p.ScoreAlleles[k]):
				p.EffectAlleles[k] = 0
				alleleSet = true
			default:
				if len(p.ScoreAlleles[k]) > 1 {
					for _, a := range p.ScoreAlleles[k] {
						if strings.Contains(alt, string(a)) {
							p.EffectAlleles[k] = 1
							alleleSet = true
							break
						}
					}
					if alleleSet {
						continue
					}
				}
				log.Printf("Effect/other and reference/alternative alleles do not match at locus %s: %v vs %v",
					locus, p.ScoreAlleles[k], []string{ref, alt})
				if !alleleSet {
					p.EffectAlleles[k] = 1
				}
				if j == len(lines)-2 {
					return errors.New("effect allele mismatch")
				}
			}
		}
		if missing {
			//	No data for locus, insert alt as the effect
			p.EffectAlleles[k] = 1
		}
	}
	return nil
}

func (p *PGS) savePopulationAlleles(filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Printf("Error creating effect alleles file: %v", err)
		return
	}
	defer file.Close()
	encoder := json.NewEncoder(file)
	err = encoder.Encode(&(p.EffectAlleles))
	if err != nil {
		log.Printf("Error encoding json effect alleles: %v", err)
	}
}

func (p *PGS) extractEAF() {
	//populationQ := "%CHROM:%POS-%" + strings.Join(ANCESTRIES, "_AF\\t%") + "_AF\n"
	populationQ := "%CHROM:%POS-%" + strings.Join(ANCESTRIES, "_AF\\t%") + "_AF\\t" + "%AF\n"
	var freq float32
	var parsed float64
	var missing bool
	populationsAndAll := append(ANCESTRIES, "ALL")
	for k, locus := range p.Loci {
		missing = true
		chr, pos := tools.SplitLocus(locus)
		query, args := tools.RangeQuery(populationQ, chr, pos, pos, tools.GG)
		cmd := exec.Command(query, args...)
		output, err := cmd.Output()
		if err != nil {
			log.Printf("Error executing bcftools command: %v", err)
			continue
		}
		lines := strings.Split(string(output), "\n")
		//fmt.Printf("Locus %s: %d lines\n", locus, len(lines))
		for _, line := range lines[:len(lines)-1] {
			fields := strings.Split(line, "-")
			if fields[0] != locus {
				fmt.Printf("Locus %s does not match %s\n", locus, fields[0])
				continue
			}
			missing = false
			//fmt.Printf("AF Locus %s: %s\n", locus, fields[1])
			afPerPopulation := strings.Split(fields[1], "\t")
			for i, population := range populationsAndAll {
				altAfs := strings.Split(afPerPopulation[i], ",")
				freq = 0
				for _, altAf := range altAfs {
					parsed, err = strconv.ParseFloat(altAf, 32)
					if err != nil {
						log.Printf("Error parsing %s at %s: %v", afPerPopulation[i], locus, err)
						parsed = MissingEAF
					}
					freq += float32(parsed)
				}
				freq /= float32(len(altAfs))
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
				if _, ok := p.PopulationStats[population].AF[k]; !ok {
					p.PopulationStats[population].AF[k] = []float32{0, freq}
				} else {
					if p.PopulationStats[population].AF[k][1]+freq > 1 {
						continue
					}
					p.PopulationStats[population].AF[k][1] += freq
				}
				if p.PopulationStats[population].AF[k][1] > 1 {
					log.Printf("++++++ Allele frequency is wrong %v for %s at %s", p.PopulationStats[population].AF[k], lines[:len(lines)-1], locus)
				}
			}
		}
		//log.Printf("No entry for locus %s, hence all the samples have reference", locus)
		for i := range populationsAndAll {
			if missing { // either read a wrong locus or len(lines[:len(lines)-1]) == 0
				p.PopulationStats[populationsAndAll[i]].AF[k] = []float32{1 - MissingEAF, MissingEAF}
				continue
			}
			// Smoothen the extreme cases
			switch p.PopulationStats[populationsAndAll[i]].AF[k][1] {
			case 0:
				p.PopulationStats[populationsAndAll[i]].AF[k][1] = MissingEAF
			case 1:
				p.PopulationStats[populationsAndAll[i]].AF[k][1] = 1 - MissingEAF
			default:
			}
			p.PopulationStats[populationsAndAll[i]].AF[k][0] = 1 - p.PopulationStats[populationsAndAll[i]].AF[k][1]
		}
	}
}

func (p *PGS) extractGF() {
	var err error
	var found bool
	var chr, position, snps, anc string
	var alleles []uint8
	var output []byte
	var counts map[string][]int
	var total map[string]int
	var positionAndSamples, idvSnps []string
	ancestries := tools.LoadAncestry()
	populationsAndAll := append(ANCESTRIES, "ALL")
	for k, locus := range p.Loci {
		found = false
		chr, position = tools.SplitLocus(locus)
		query, args := tools.IndividualSnpsQuery(chr, position, tools.GG)
		output, err = exec.Command(query, args...).Output()
		if err != nil {
			fmt.Println("Error executing bcftools command:", err)
			continue
		}
		total = make(map[string]int)
		counts = make(map[string][]int)
		for _, ppl := range populationsAndAll {
			total[ppl] = 0
			counts[ppl] = make([]int, 3)
		}
		lines := strings.Split(string(output), "\n")
		for _, line := range lines[:len(lines)-1] {
			positionAndSamples = strings.Split(line, "*")
			if positionAndSamples[0] != locus {
				//fmt.Printf("Locus %s does not match %s\n", locus, positionAndSamples[0])
				continue
			}
			found = true
			samples := strings.Split(positionAndSamples[1], "\t")
			samples = samples[:len(samples)-1]
			for _, sample := range samples {
				idvSnps = strings.Split(sample, "=")
				snps, err = tools.NormalizeSnp(idvSnps[1])
				if err != nil {
					fmt.Printf("Error normalizing %s: %v", snps, err)
					continue
				}
				alleles, err = tools.SnpToPair(snps)
				if err != nil {
					fmt.Printf("Error converting SNPs %s to alleles: %v", snps, err)
					continue
				}
				anc = GetIndividualAncestry(idvSnps[0], ancestries)
				counts[anc][alleles[0]+alleles[1]]++
				total[anc]++
				counts["ALL"][alleles[0]+alleles[1]]++
				total["ALL"]++
			}
		}
		for _, ppl := range populationsAndAll {
			p.PopulationStats[ppl].GF[k] = make([]float32, 3)
			if !found {
				p.PopulationStats[ppl].GF[k][0] = p.PopulationStats[ppl].AF[k][0] * p.PopulationStats[ppl].AF[k][0]
				p.PopulationStats[ppl].GF[k][1] = 2 * p.PopulationStats[ppl].AF[k][0] * p.PopulationStats[ppl].AF[k][1]
				p.PopulationStats[ppl].GF[k][2] = p.PopulationStats[ppl].AF[k][1] * p.PopulationStats[ppl].AF[k][1]
				continue
			}
			for i := 0; i < 3; i++ {
				p.PopulationStats[ppl].GF[k][i] = float32(counts[ppl][i]) / float32(total[ppl])
			}
		}
	}
}

func (p *PGS) CalculateEAFAndGF(dataset string) {
	var err error
	var found bool
	var chr, position, snps string
	var output []byte
	var pair []uint8
	var positionAndSamples, idvSnps []string
	var genotypeCount, alleleCount map[uint8]int
	for k, locus := range p.Loci {
		found = false
		chr, position = tools.SplitLocus(locus)
		query, args := tools.IndividualSnpsQuery(chr, position, dataset)
		output, err = exec.Command(query, args...).Output()
		if err != nil {
			fmt.Println("Error executing bcftools command:", err)
			continue
		}
		genotypeCount = make(map[uint8]int)
		alleleCount = make(map[uint8]int)
		lines := strings.Split(string(output), "\n")
		for _, line := range lines[:len(lines)-1] {
			positionAndSamples = strings.Split(line, "*")
			if positionAndSamples[0] != locus {
				continue
			}
			found = true
			samples := strings.Split(positionAndSamples[1], "\t")
			samples = samples[:len(samples)-1]
			for _, sample := range samples {
				idvSnps = strings.Split(sample, "=")
				snps, err = tools.NormalizeSnp(idvSnps[1])
				if err != nil {
					fmt.Printf("Error normalizing %s: %v", idvSnps[1], err)
					continue
				}
				pair, err = tools.SnpToPair(snps)
				if err != nil {
					fmt.Printf("Error converting SNPs %s to alleles: %v", snps, err)
					continue
				}
				genotypeCount[pair[0]+pair[1]]++
				alleleCount[pair[0]]++
				alleleCount[pair[1]]++
			}
		}
		if !found {
			if _, ok := p.PopulationStats["ALL"].AF[k]; ok {
				// Already recorded before due to possibly misaligned index
				continue
			}
			p.PopulationStats["ALL"].AF[k] = []float32{1 - MissingEAF, MissingEAF}
			p.PopulationStats["ALL"].GF[k] = []float32{(1 - MissingEAF) * (1 - MissingEAF), 2 * MissingEAF * (1 - MissingEAF),
				MissingEAF * MissingEAF}
			continue
		}
		totalAlleles := alleleCount[0] + alleleCount[1]
		p.PopulationStats["ALL"].AF[k] = []float32{float32(alleleCount[0]) / float32(totalAlleles),
			float32(alleleCount[1]) / float32(totalAlleles)}
		totalGenotypes := genotypeCount[0] + genotypeCount[1] + genotypeCount[2]
		p.PopulationStats["ALL"].GF[k] = []float32{float32(genotypeCount[0]) / float32(totalGenotypes),
			float32(genotypeCount[1]) / float32(totalGenotypes), float32(genotypeCount[2]) / float32(totalGenotypes)}
	}
}

//func (p *PGS) computeFrequencySpectrum() {
//	p.assignFreqBins()
//	p.querySampleFrequencies()
//}

//func (p *PGS) assignFreqBins() {
//	for _, ppl := range ANCESTRIES {
//		eaf := make([]float64, len(p.PopulationStats[ppl].AF))
//		for i := range p.PopulationStats[ppl].AF {
//			eaf[i] = float64(p.PopulationStats[ppl].AF[i][p.EffectAlleles[i]])
//		}
//		sort.Float64s(eaf)
//		elementsPerBin := len(eaf) / p.NumSpecBins
//		var endIdx int
//		for i := 0; i < p.NumSpecBins; i++ {
//			endIdx = (i + 1) * elementsPerBin
//			// If we're at the last bin, adjust the endIndex to include all remaining elements
//			if i == p.NumSpecBins-1 {
//				endIdx = len(eaf)
//			}
//			p.PopulationStats[ppl].FreqBinBounds[i] = float32(eaf[endIdx-1])
//		}
//	}
//}

//func (p *PGS) querySampleFrequencies() {
//	ancestry := tools.LoadAncestry()
//	ancestrySizes := make(map[string]int, len(ancestry))
//	var popNames []string
//	for _, popNameStr := range ancestry {
//		// Some samples are assigned multiple comma-separated ancestries
//		popNames = strings.Split(popNameStr, ",")
//		for _, popName := range popNames {
//			if _, exists := ancestrySizes[popName]; !exists {
//				ancestrySizes[popName] = 0
//			}
//			ancestrySizes[popName]++
//		}
//	}
//	var normalized, sample, anc string
//	var ancs []string
//	var val int
//	for k, locus := range p.Loci {
//		chr, pos := tools.SplitLocus(locus)
//		query, args := tools.RangeGenotypesQuery(chr, pos, pos)
//		cmd := exec.Command(query, args...)
//		output, err := cmd.Output()
//		if err != nil {
//			log.Printf("Error executing bcftools command: %v", err)
//			continue
//		}
//		lines := strings.Split(string(output), "\n")
//		if len(lines[:len(lines)-1]) == 0 {
//			log.Printf("No data for locus %s in frequency spectrum", locus)
//			for _, ppl := range ANCESTRIES {
//				// Assume that all the samples have the reference allele, and one "ghost" sample has the effect allele
//				p.PopulationStats[ppl].FreqSpectrum[tools.ValueToBinIdx(p.PopulationStats[ppl].AF[k][1],
//					p.PopulationStats[ppl].FreqBinBounds)] += 1
//			}
//			continue
//		}
//		for _, line := range lines[:len(lines)-1] {
//			fields := strings.Split(line, "-")
//			if fields[0] != locus {
//				fmt.Printf("Locus %s does not match %s\n", locus, fields[0])
//				continue
//			}
//			samples := strings.Split(fields[1], "\t")
//			for _, sample = range samples[:len(samples)-1] {
//				gtp := strings.Split(sample, "=")
//				ancs = strings.Split(ancestry[gtp[0]], ",")
//				alleles := strings.Split(gtp[1], "|")
//				for _, allele := range alleles {
//					normalized, err = tools.NormalizeAllele(allele)
//					if err != nil {
//						if normalized != "." {
//							log.Printf("%v: %s, at %s\n", err, sample, locus)
//						}
//						continue
//					}
//					val, err = strconv.Atoi(normalized)
//					if err != nil {
//						log.Printf("Error converting %s to int: %v", normalized, err)
//						continue
//					}
//					if uint8(val) != p.EffectAlleles[k] {
//						continue
//					}
//					for _, anc = range ancs {
//						p.PopulationStats[anc].FreqSpectrum[tools.ValueToBinIdx(p.PopulationStats[anc].AF[k][p.EffectAlleles[k]],
//							p.PopulationStats[anc].FreqBinBounds)]++
//					}
//				}
//			}
//		}
//	}
//	// Normalize the frequency spectrum
//	for _, ppl := range ANCESTRIES {
//		for i := 0; i < p.NumSpecBins; i++ {
//			p.PopulationStats[ppl].FreqSpectrum[i] /= float32(ancestrySizes[ppl])
//		}
//	}
//}

func (p *PGS) SaveStats(filename string) {
	file, err := os.Create(filename)
	if err != nil {
		log.Printf("Error creating stats file: %v", err)
		return
	}
	defer file.Close()
	encoder := json.NewEncoder(file)
	err = encoder.Encode(p.PopulationStats)
	if err != nil {
		log.Printf("Error encoding json PopulationStats: %v", err)
	}
}

func GetIndividualAncestry(idv string, ancestry map[string]string) string {
	if _, ok := ancestry[idv]; !ok {
		log.Fatalf("Individual %s not found in ancestry file\n", idv)
	}
	idvAnc := ancestry[idv]
	if strings.Contains(idvAnc, ",") {
		idvAnc = strings.Split(idvAnc, ",")[0]
	}
	return idvAnc
}

func (p *PGS) FindMaxAbsoluteWeight() float64 {
	maxWeight := 0.0
	for _, weight := range p.Weights {
		w, err := weight.Float64()
		if err != nil {
			log.Println("Error converting weight to float64:", err)
		}
		if math.Abs(w) > maxWeight {
			maxWeight = math.Abs(w)
		}
	}
	return maxWeight
}

func (p *PGS) FindMinAbsoluteWeight() float64 {
	minWeight := 0.0
	for _, weight := range p.Weights {
		w, err := weight.Float64()
		if err != nil {
			log.Println("Error converting weight to float64:", err)
		}
		if math.Abs(w) < minWeight {
			minWeight = math.Abs(w)
		}
	}
	return minWeight
}

func (p *PGS) FindMinAndMaxScores() (*apd.Decimal, *apd.Decimal, *apd.Decimal, *apd.Decimal) {
	maxScore, minScore := apd.New(0, 0), apd.New(0, 0)
	secondMaxScore, secondMinScore := new(apd.Decimal), new(apd.Decimal)
	var smallestPositiveWeight, smallestNegativeWeight *apd.Decimal
	allPositive, allNegative := true, true
	for _, weight := range p.Weights {
		if weight.Negative {
			p.Context.Add(minScore, minScore, weight)
			p.Context.Add(minScore, minScore, weight)
			allPositive = false
			if smallestNegativeWeight == nil || weight.Cmp(smallestNegativeWeight) > 0 {
				smallestNegativeWeight = weight
			}
			continue
		}
		p.Context.Add(maxScore, maxScore, weight)
		p.Context.Add(maxScore, maxScore, weight)
		allNegative = false
		if smallestPositiveWeight == nil || weight.Cmp(smallestPositiveWeight) < 0 {
			smallestPositiveWeight = weight
		}
	}
	var err error
	if allNegative {
		_, err = p.Context.Add(secondMaxScore, maxScore, smallestNegativeWeight)
	} else {
		_, err = p.Context.Sub(secondMaxScore, maxScore, smallestPositiveWeight)
	}
	if err != nil {
		log.Printf("Error computing second max score %s: %v\n", maxScore.String(), err)
		return nil, nil, nil, nil
	}
	if allPositive {
		secondMinScore.Set(smallestPositiveWeight)
	} else {
		_, err = p.Context.Sub(secondMinScore, minScore, smallestNegativeWeight)
		if err != nil {
			log.Printf("Error computing second max score %s: %v\n", maxScore.String(), err)
			return nil, nil, nil, nil
		}
	}
	divisor := new(apd.Decimal).SetInt64(int64(len(p.Weights) * Ploidy))
	_, err = p.Context.Quo(maxScore, maxScore, divisor)
	if err != nil {
		log.Printf("Error normalizing max score %s: %v\n", maxScore.String(), err)
		return nil, nil, nil, nil
	}
	_, err = p.Context.Quo(minScore, minScore, divisor)
	if err != nil {
		log.Printf("Error normalizing min score %s: %v\n", minScore.String(), err)
		return nil, nil, nil, nil
	}
	_, err = p.Context.Quo(secondMaxScore, secondMaxScore, divisor)
	if err != nil {
		log.Printf("Error normalizing second max score %s: %v\n", secondMaxScore.String(), err)
		return nil, nil, nil, nil
	}
	_, err = p.Context.Quo(secondMinScore, secondMinScore, divisor)
	if err != nil {
		log.Printf("Error normalizing second min score %s: %v\n", secondMinScore.String(), err)
		return nil, nil, nil, nil
	}
	return minScore, maxScore, secondMinScore, secondMaxScore
}

func (p *PGS) EstimateMeanAndStd() (map[string]float64, map[string]float64) {
	means, stds := make(map[string]float64), make(map[string]float64)
	prob, mui := 0.0, 0.0
	for i, weight := range p.Weights {
		w, err := weight.Float64()
		if err != nil {
			log.Println("Error converting weight to float64:", err)
		}
		for anc := range p.PopulationStats {
			if len(p.PopulationStats[anc].AF[i]) == 0 {
				continue
			}
			prob = float64(p.PopulationStats[anc].AF[i][p.EffectAlleles[i]])
			mui = prob*prob*(2*w) + 2*prob*(1-prob)*w
			if _, ok := means[anc]; !ok {
				means[anc] = 0.0
				stds[anc] = 0.0
			}
			means[anc] += mui
			stds[anc] += (1-prob)*(1-prob)*math.Pow(0-mui, 2) + 2*prob*(1-prob)*math.Pow(w-mui, 2) +
				prob*prob*math.Pow(2*w-mui, 2)
		}
	}
	for anc := range p.PopulationStats {
		if _, ok := means[anc]; !ok {
			continue
		}
		means[anc] /= float64(Ploidy * len(p.Weights))
		stds[anc] = math.Sqrt(stds[anc]) / float64(Ploidy*len(p.Weights))
	}
	return means, stds
}
