package pgs

import (
	"bufio"
	"encoding/csv"
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

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
	NUM_HAPLOTYPES = 2
	NUM_GENOTYPES  = 2
)

var GENOTYPES = []int{0, 1}

type Variant struct {
	fields map[string]interface{}
}

func NewVariant(fields map[string]interface{}) *Variant {
	v := &Variant{
		fields: make(map[string]interface{}),
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
	PgsID        string
	TraitName    string
	TraitEFO     string
	GenomeBuild  string
	WeightType   string
	HmPOSBuild   string
	VariantCount int
	Fieldnames   []string
	Variants     map[string]*Variant
	Loci         []string
	Weights      []float64
	Maf          [][]float64
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
					p.VariantCount = value
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
		minChr, minPos, err := tools.ParseLocus(sortedLoc[minIndex])
		if err != nil {
			log.Printf("Error parsing locus %s: %v", sortedLoc[minIndex], err)
			return nil, err
		}
		for j := i + 1; j < len(sortedLoc); j++ {
			chr, pos, err := tools.ParseLocus(sortedLoc[j])
			if err != nil {
				log.Printf("Error parsing locus %s: %v", sortedLoc[minIndex], err)
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

func (p *PGS) LoadStats() {
	p.LoadMAF()
}

func (p *PGS) LoadMAF() {
	filename := fmt.Sprintf("%s.maf", p.PgsID)
	_, err := os.Stat(filename)
	if os.IsNotExist(err) {
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

			for {
				row, err := reader.Read()
				if err != nil {
					break // Reached end of file or encountered an error
				}
				if row[0] != pos {
					continue
				}
				writer.Write([]string{chr, pos, row[1]})
				break
			}
			f.Close()
			writer.Flush()
		}
		file.Close()
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
		likelihood, err := strconv.ParseFloat(row[2], 64)
		if err != nil {
			log.Printf("Error converting %s to float: %v", row[2], err)
			continue
		}
		p.Maf = append(p.Maf, []float64{likelihood, 1 - likelihood})
	}
}

func (p *PGS) MutateGenome(original []int, delta float64) []int {
	mutations := make([]int, 0, len(original))
	probabilities := make([]float64, 0, len(mutations))
	tmp := 0.0
	for i := 0; i < len(original); i += NUM_HAPLOTYPES {
		for j := 0; j < NUM_HAPLOTYPES; j++ {
			tmp = delta - float64(original[i+j])*p.Weights[i/2]
			for _, v := range GENOTYPES {
				if v == original[i+j] {
					continue
				}
				if tmp+float64(v)*p.Weights[i/2] == 0 {
					original[i+j] = v
					return original
				}
				probabilities = append(probabilities, 1/math.Abs(tmp+float64(v)*p.Weights[i/2]))
				mutations = append(mutations, v)
			}
		}
	}
	//fmt.Println(probabilities)
	mutationId := tools.SampleFromDistribution(probabilities)
	original[mutationId] = mutations[mutationId]
	return original
}

//func (p *PGS) MutateGenome(original [][]int) [][]int {
//	mutations := make([]int, 0, NUM_HAPLOTYPES*len(GENOTYPES)*len(original))
//	probabilities := make([]float64, 0, NUM_HAPLOTYPES*len(GENOTYPES)*len(original))
//	prob := 0.0
//	for i := range original {
//		// We range over all the possible genotypes, even the present one, to allow for the possibility of no mutation
//		for j := 0; j < NUM_HAPLOTYPES; j++ {
//			for _, v := range GENOTYPES {
//				// Get individual prior
//				prob = p.Maf[i][v]
//				probabilities = append(probabilities, prob)
//				mutations = append(mutations, v)
//			}
//		}
//	}
//	mutationId := tools.SampleFromDistribution(probabilities)
//	original[mutationId/(NUM_HAPLOTYPES*len(GENOTYPES))][(mutationId/len(GENOTYPES))%NUM_HAPLOTYPES] = mutations[mutationId]
//	return original
//}

// Sample according to the MAF
func (p *PGS) SampleFromPopulation() ([]int, error) {
	sample := make([]int, p.VariantCount*NUM_HAPLOTYPES)
	// Initial sample based on individual priors
	for i := range p.Loci {
		for j := 0; j < NUM_HAPLOTYPES; j++ {
			maf := p.Maf[i]
			ind := tools.SampleFromDistribution(maf)
			sample[i*NUM_HAPLOTYPES+j] = GENOTYPES[ind]
			if sample[i*NUM_HAPLOTYPES+j] == -1 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) CalculateSequenceLikelihood(sequence []int) float64 {
	likelihood := 1.0
	for i := range p.Loci {
		for j := 0; j < NUM_HAPLOTYPES; j++ {
			likelihood += math.Log(p.Maf[i][sequence[i*NUM_HAPLOTYPES+j]])
		}
	}
	return likelihood
}
