package main

import (
	"bufio"
	"log"
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

func (v *Variant) GetID() string {
	return v.fields["rsID"].(string)
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
}

func NewPGS() *PGS {
	return &PGS{
		Variants: make(map[string]*Variant),
	}
}

func (p *PGS) Load(inputFile string) error {
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
		p.Variants[variant.GetID()] = variant
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	return nil
}
