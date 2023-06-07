package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"
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

func (v *Variant) GetLocation() string {
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
	Locations      []string
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
		p.Variants[variant.GetLocation()] = variant
	}
	p.Locations = p.GetSortedVariantLocations()
	p.Weights = make([]float64, len(p.Locations))
	for i, loc := range p.Locations {
		p.Weights[i] = p.Variants[loc].GetWeight()
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	return nil
}

func (p *PGS) GetSortedVariantLocations() []string {
	sortedLoc := make([]string, 0, len(p.Variants))
	for location := range p.Variants {
		sortedLoc = append(sortedLoc, location)
	}
	for i := 0; i < len(sortedLoc)-1; i++ {
		minIndex := i
		for j := i + 1; j < len(sortedLoc); j++ {
			if sortedLoc[j] < sortedLoc[minIndex] {
				minIndex = j
			}
		}
		if minIndex != i {
			sortedLoc[i], sortedLoc[minIndex] = sortedLoc[minIndex], sortedLoc[i]
		}
	}
	return sortedLoc
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
				fmt.Println(row)
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
		location := fmt.Sprintf("%s:%s", chr, pos)
		for i := 2; i < len(row); i++ {
			key, err1 := strconv.Atoi(header[i])
			value, err2 := strconv.ParseFloat(row[i], 64)
			if err1 != nil || err2 != nil {
				log.Printf("Error converting key %s or value %s to int: %v, %v", header[i], row[i], err1, err2)
				continue
			}
			p.Variants[location].priors[key] = value
		}
	}
}

func (p *PGS) MutateVariant(original, mutations []int, probabilities []float64) []int {
	candidatesPerVariant := len(mutations) / len(original)
	for i := range original {
		priors := p.Variants[p.Locations[i]].priors
		currentVariant := i * candidatesPerVariant
		for j := 0; j < candidatesPerVariant; j++ {
			probabilities[currentVariant+j] = probabilities[currentVariant+j] * priors[mutations[currentVariant+j]]
		}
	}
	mutationId := sampleSlice(probabilities)
	original[mutationId/candidatesPerVariant] = mutations[mutationId]
	return original
}

func (p *PGS) GetVariantPriors(location string) map[int]float64 {
	return p.Variants[location].priors
}

func sampleSlice(distribution []float64) int {
	rand.NewSource(time.Now().UnixNano())
	cumulative := 0.0
	for _, p := range distribution {
		cumulative += p
	}
	r := rand.Float64() * cumulative
	for k, v := range distribution {
		r -= v
		if r <= 0.0 {
			return k
		}
	}
	return -1
}
