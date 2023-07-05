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

var GENOTYPES = []string{"0|0", "0|1", "1|0", "1|1"}

type Variant struct {
	fields map[string]interface{}
	priors map[string]float64
}

func NewVariant(fields map[string]interface{}) *Variant {
	v := &Variant{
		fields: make(map[string]interface{}),
		priors: make(map[string]float64),
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

func (v *Variant) GetPriors() map[string]float64 {
	return v.priors
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
	correlations   [][]float64
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
			log.Printf("Error parsing locus %s: %v", sortedLoc[minIndex], err)
			return nil, err
		}
		for j := i + 1; j < len(sortedLoc); j++ {
			chr, pos, err := parseLocus(sortedLoc[j])
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

func (p *PGS) LoadPriors() {
	p.LoadIndividualPriors()
	//p.LoadCorrelations()
}

func (p *PGS) LoadIndividualPriors() {
	_, err := os.Stat(p.PgsID + ".priors")
	if os.IsNotExist(err) {
		priorsFile, err := os.Create(p.PgsID + ".priors")
		if err != nil {
			log.Printf("Error creating priors file: %v", err)
			return
		}
		defer priorsFile.Close()
		writer := csv.NewWriter(priorsFile)
		writer.Write([]string{"chromosome", "position", "0|0", "0|1", "1|0", "1|1"})
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
				writer.Write([]string{chr, pos, row[1], row[2], row[3], row[4]})
				break
			}
			f.Close()
		}
		writer.Flush()
		priorsFile.Close()
	}
	priorsFile, err := os.Open(p.PgsID + ".priors")
	if err != nil {
		log.Printf("Error opening priors file: %v\n", err)
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
			value, err := strconv.ParseFloat(row[i], 64)
			if err != nil {
				log.Printf("Error converting key %s or value %s to int: %v", header[i], row[i], err)
				continue
			}
			p.Variants[locus].priors[header[i]] = value
		}
	}
}

func (p *PGS) LoadCorrelations() {
	file, err := os.Open("data/prior/" + p.PgsID + ".pairwise")
	if err != nil {
		log.Printf("Error opening priors file: %v", err)
		return
	}
	defer file.Close()
	p.correlations = make([][]float64, len(GENOTYPES)*p.VariantsNumber)
	reader := csv.NewReader(file)
	j := 0
	for {
		row, err := reader.Read()
		if err != nil {
			break // Reached end of file or encountered an error
		}
		p.correlations[j] = make([]float64, len(row))
		for i := 0; i < len(row); i++ {
			prob, err := strconv.ParseFloat(row[i], 64)
			if err != nil {
				log.Printf("Error converting probabilitiy to float: %s, %v", row[i], err)
				continue
			}
			p.correlations[j][i] = prob
		}
		j++
	}
}

func (p *PGS) MutateGenome(original []string) []string {
	mutations := make([]string, 0, len(GENOTYPES)*len(original))
	probabilities := make([]float64, 0, len(GENOTYPES)*len(original))
	prob := 0.0
	for i := range original {
		// We range over all the possible genotypes, even the present one, to allow for the possibility of no mutation
		for j, v := range GENOTYPES {
			// Get individual prior
			prob = math.Exp(p.Variants[p.Loci[i]].priors[v])
			// Add conditional probabilities based on the rest of the genome
			for k := range original {
				if k == i {
					continue
				}
				// Given the SNP value at position k, what is the probability of the SNP value at position i being j
				prob += math.Exp(p.correlations[k*len(GENOTYPES)+original[k]][i*len(GENOTYPES)+j])
			}
			probabilities = append(probabilities, prob)
			mutations = append(mutations, v)
		}
	}
	mutationId := tools.SampleFromDistribution(probabilities)
	original[mutationId/len(GENOTYPES)] = mutations[mutationId]
	return original
}

func (p *PGS) GetVariantPriors(locus string) map[string]float64 {
	return p.Variants[locus].priors
}

// We use the Gibbs sampling approach: first sample based only on the individual priors,
// then iterate by taking into account correlations correlations
func (p *PGS) SampleFromPopulation() ([]string, error) {
	sample := make([]string, 2*len(p.Variants))
	// Initial sample based on individual priors
	for i, loc := range p.Loci {
		sample[i] = tools.SampleFromMap(p.Variants[loc].priors)
		if sample[i] == "" {
			return nil, errors.New("error in population sampling")
		}
	}
	//// Iterate to convergence
	//NUM_ITERATIONS := 100
	//for i := 0; i < NUM_ITERATIONS; i++ {
	//
	//}
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
