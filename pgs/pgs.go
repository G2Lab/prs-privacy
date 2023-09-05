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
	NumHaplotypes = 2
	ErrorMargin   = 1e-14
	Precision    = 1e-15
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

	headerInProgress := true
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
		if headerInProgress {
			p.Fieldnames = strings.Split(line, "\t")
			headerInProgress = false
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
				writer.Write([]string{chr, pos, "0.99999"})
				log.Printf("No MAF found for locus %s", locus)
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

func (p *PGS) MutateGenome(original []int, delta float64) [][]int {
	mutations := make([]int, len(original))
	probabilities := make([]float64, len(original))
	mutated := make([][]int, 0, 1)
	tmp := 0.0
	for i := 0; i < len(original)/NumHaplotypes; i++ {
		for j := 0; j < NumHaplotypes; j++ {
			tmp = delta - float64(original[NumHaplotypes*i+j])*p.Weights[i]
			for _, v := range GENOTYPES {
				if v == original[NumHaplotypes*i+j] {
					continue
				}
				//if tmp+float64(v)*p.Weights[i] == 0 {
				if math.Abs(tmp+float64(v)*p.Weights[i]) <= ErrorMargin {
					match := make([]int, len(original))
					copy(match, original)
					match[NumHaplotypes*i+j] = v
					mutated = append(mutated, match)
					continue
				}
				//probabilities[NumHaplotypes*i+j] = 1 / math.Abs(tmp+float64(v)*p.Weights[i])
				probabilities[NumHaplotypes*i+j] = p.Maf[i][v]
				mutations[NumHaplotypes*i+j] = v
			}
		}
	}
	if len(mutated) > 0 {
		// we found solutions so no need to sample
		return mutated
	}
	mutationId := tools.SampleFromDistribution(probabilities)
	original[mutationId] = mutations[mutationId]
	mutated = append(mutated, original)
	return mutated
}

//func (p *PGS) MutateGenome(original [][]int) [][]int {
//	mutations := make([]int, 0, NumHaplotypes*len(GENOTYPES)*len(original))
//	probabilities := make([]float64, 0, NumHaplotypes*len(GENOTYPES)*len(original))
//	prob := 0.0
//	for i := range original {
//		// We range over all the possible genotypes, even the present one, to allow for the possibility of no mutation
//		for j := 0; j < NumHaplotypes; j++ {
//			for _, v := range GENOTYPES {
//				// Get individual prior
//				prob = p.Maf[i][v]
//				probabilities = append(probabilities, prob)
//				mutations = append(mutations, v)
//			}
//		}
//	}
//	mutationId := tools.SampleFromDistribution(probabilities)
//	original[mutationId/(NumHaplotypes*len(GENOTYPES))][(mutationId/len(GENOTYPES))%NumHaplotypes] = mutations[mutationId]
//	return original
//}

// Sample according to the MAF
func (p *PGS) SampleFromPopulation() ([]int, error) {
	sample := make([]int, p.VariantCount*NumHaplotypes)
	// Initial sample based on individual priors
	for i := range p.Loci {
		for j := 0; j < NumHaplotypes; j++ {
			maf := p.Maf[i]
			ind := tools.SampleFromDistribution(maf)
			sample[i*NumHaplotypes+j] = GENOTYPES[ind]
			if sample[i*NumHaplotypes+j] == -1 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) SampleUniform() ([]int, error) {
	sample := make([]int, p.VariantCount*NumHaplotypes)
	// Initial sample based on individual priors
	for i := range p.Loci {
		for j := 0; j < NumHaplotypes; j++ {
			sample[i*NumHaplotypes+j] = tools.SampleUniform(GENOTYPES)
			if sample[i*NumHaplotypes+j] == -1 {
				return nil, errors.New("error in population sampling")
			}
		}
	}
	return sample, nil
}

func (p *PGS) CalculateSequenceLikelihood(sequence []int) float64 {
	likelihood := 1.0
	for i := range p.Loci {
		for j := 0; j < NumHaplotypes; j++ {
			likelihood += math.Log(p.Maf[i][sequence[i*NumHaplotypes+j]])
		}
	}
	return likelihood
}
