package main

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"github.com/nikirill/prs/tools"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
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

var GENOTYPES = []int{0, 1, 2}

type Variant struct {
	fields map[string]interface{}
	priors map[int]float64
}

func NewVariant(fields map[string]interface{}) *Variant {
	v := &Variant{
		fields: make(map[string]interface{}),
		priors: make(map[int]float64),
	}

	if weight, ok := fields["effect_weight"].(string); ok {
		if value, err := strconv.ParseFloat(weight, 64); err == nil {
			fields["effect_weight"] = value
		} else {
			log.Printf("Error parsing weight %s: %s", weight, err)
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

func (v *Variant) GetWeight() float64 {
	if weight, ok := v.fields["effect_weight"].(float64); ok {
		return weight
	}
	return 0.0
}

type PGS struct {
	PgsID          string
	TraitName      string
	TraitEFO       string
	GenomeBuild    string
	WeightType     string
	HmPOSBuild     string
	VariantsNumber int
	Fieldnames     []string
	Variants       map[string]*Variant
	Loci           []string
	Weights        []float64
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

	scanner := bufio.NewScanner(file)
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
					p.VariantsNumber = value
				} else {
					log.Printf("Error parsing variants number %s: %s", fields[1], err)
				}
			}
			continue
		}

		// Specified variant fields
		if strings.HasPrefix(line, "rsID") {
			p.Fieldnames = strings.Split(line, "\t")
			continue
		}

		fields := make(map[string]interface{})
		values := strings.Split(line, "\t")
		for i, value := range values {
			fields[p.Fieldnames[i]] = value
		}
		variant := NewVariant(fields)
		p.Variants[variant.GetLocus()] = variant
	}
	p.Loci, err = p.GetSortedVariantLoci()
	if err != nil {
		return err
	}
	p.Weights = make([]float64, len(p.Loci))
	for i, loc := range p.Loci {
		p.Weights[i] = p.Variants[loc].GetWeight()
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	return nil
}

func (p *PGS) GetSortedVariantLoci() ([]string, error) {
	sortedLoc := make([]string, 0, len(p.Variants))
	for locus := range p.Variants {
		sortedLoc = append(sortedLoc, locus)
	}
	for i := 0; i < len(sortedLoc)-1; i++ {
		minIndex := i
		minChr, minPos, err := parseLocus(sortedLoc[minIndex])
		if err != nil {
			log.Printf("Error parsing locus %s: %v, %v", sortedLoc[minIndex], err)
			return nil, err
		}
		for j := i + 1; j < len(sortedLoc); j++ {
			chr, pos, err := parseLocus(sortedLoc[j])
			if err != nil {
				log.Printf("Error parsing locus %s: %v, %v", sortedLoc[minIndex], err)
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

func (p *PGS) LoadPriors() {
	_, err := os.Stat(p.PgsID + ".priors")
	if os.IsNotExist(err) {
		priorsFile, error := os.Create(p.PgsID + ".priors")
		if error != nil {
			log.Printf("Error creating priors file: %s", error)
			return
		}
		defer priorsFile.Close()
		writer := csv.NewWriter(priorsFile)
		writer.Write([]string{"chromosome", "position", "0", "1", "2"})
		for _, variant := range p.Variants {
			chr := variant.GetHmChr()
			pos := variant.GetHmPos()
			f, err := os.Open(fmt.Sprintf("data/prior/chr%s.csv", chr))
			if err != nil {
				log.Printf("Error opening file chr%s.csv: %s", chr, err)
				continue
			}
			reader := csv.NewReader(f)
			// Read header
			_, err = reader.Read()
			if err != nil {
				fmt.Printf("Error reading priors file chr%s.csv: %v", chr, err)
				continue
			}

			for {
				row, err := reader.Read()
				if err != nil {
					break // Reached end of file or encountered an error
				}
				if row[0] != pos {
					continue
				}
				zeroOneProb, err1 := strconv.ParseFloat(row[2], 64)
				oneZeroProb, err2 := strconv.ParseFloat(row[3], 64)
				if err1 != nil || err2 != nil {
					log.Printf("Error converting probabilities from position %s in chr%s.csv: %s", pos, chr, err)
					continue
				}
				writer.Write([]string{chr, pos, row[1], fmt.Sprintf("%.5f", zeroOneProb+oneZeroProb), row[4]})
				break
			}
			f.Close()
		}
		writer.Flush()
		priorsFile.Close()
	}
	priorsFile, err := os.Open(p.PgsID + ".priors")
	if err != nil {
		log.Printf("Error opening priors file: %v", err)
		return
	}
	defer priorsFile.Close()
	reader := csv.NewReader(priorsFile)
	// Read header
	header, err := reader.Read()
	if err != nil {
		log.Printf("Error reading header: %v", err)
		return
	}
	for {
		row, err := reader.Read()
		if err != nil {
			break // Reached end of file or encountered an error
		}
		chr := row[0]
		pos := row[1]
		locus := fmt.Sprintf("%s:%s", chr, pos)
		for i := 2; i < len(row); i++ {
			key, err1 := strconv.Atoi(header[i])
			value, err2 := strconv.ParseFloat(row[i], 64)
			if err1 != nil || err2 != nil {
				log.Printf("Error converting key %s or value %s to int: %v, %v", header[i], row[i], err1, err2)
				continue
			}
			p.Variants[locus].priors[key] = value
		}
	}
}

func (p *PGS) MutateGenome(original []int) []int {
	mutations := make([]int, 0, len(GENOTYPES)*len(original))
	probabilities := make([]float64, 0, len(GENOTYPES)*len(original))
	for i := range original {
		// We range over all the possible genotypes, even the present one, to allow for the possibility of no mutation
		for _, v := range GENOTYPES {
			probabilities = append(probabilities, p.Variants[p.Loci[i]].priors[v])
			mutations = append(mutations, v)
		}
	}
	mutationId := tools.SampleFromDistribution(probabilities)
	original[mutationId/len(GENOTYPES)] = mutations[mutationId]
	return original
}

func (p *PGS) GetVariantPriors(locus string) map[int]float64 {
	return p.Variants[locus].priors
}

func (p *PGS) SampleFromPopulation() ([]int, error) {
	sample := make([]int, len(p.Variants))
	//for i := range p.Loci {
	//	sample[i] = tools.SampleUniform(GENOTYPES)
	for i, loc := range p.Loci {
		sample[i] = tools.SampleFromMap(p.Variants[loc].priors)
		if sample[i] == -1 {
			return nil, errors.New("error in population sampling")
		}
	}
	return sample, nil
}

func parseLocus(locus string) (int, int, error) {
	chr, err := strconv.Atoi(strings.Split(locus, ":")[0])
	if err != nil {
		return 0, 0, err
	}
	pos, err := strconv.Atoi(strings.Split(locus, ":")[1])
	if err != nil {
		return 0, 0, err
	}
	return chr, pos, nil
}

func (p *PGS) CalculateSequenceLikelihood(seq []int) float64 {
	likelihood := 0.0
	for i, locus := range p.Loci {
		likelihood += math.Log(p.Variants[locus].priors[seq[i]])
	}
	return likelihood
}
