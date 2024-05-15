package tools

import (
	"encoding/json"
	"errors"
	"fmt"
	"log"
	"os"
)

const (
	GG = "1000Genomes"
	RL = "Relatives"
	//UKBB = "UKBB"
)

func GetChromosomeFilepath(chr, dataset string) string {
	switch dataset {
	case GG:
		path := "/gpfs/commons/datasets/1000genomes/GRCh37/"
		return path + "ALL.chr" + chr + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
	case RL:
		path := "/gpfs/commons/datasets/1000genomes/release-20130502-supporting/related_samples_vcf/"
		return path + "ALL.chr" + chr + ".phase3_shapeit2_mvncall_integrated_v5_related_samples.20130502.genotypes.vcf.gz"
	default:
		log.Fatalln("Unknown dataset:", dataset)
		return ""
	}
}

func IndividualSnpsQuery(c, p, dataset string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		"%CHROM:%POS-[%SAMPLE=%GT\t]\n",
		"-r",
		fmt.Sprintf("%s:%s-%s", c, p, p),
		GetChromosomeFilepath(c, dataset),
	}
}

func RangeQuery(searchPattern, c, posBegin, posEnd, dataset string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		searchPattern,
		"-r",
		fmt.Sprintf("%s:%s-%s", c, posBegin, posEnd),
		GetChromosomeFilepath(c, dataset),
	}
}

func RangeSnpValuesQuery(c, posBegin, posEnd, dataset string) (string, []string) {
	return RangeQuery("[%GT\t]\n", c, posBegin, posEnd, dataset)
}

//func RangeGenotypesQuery(c, posBegin, posEnd string) (string, []string) {
//	return RangeQuery("%CHROM:%POS-[%SAMPLE=%GT\t]\n", c, posBegin, posEnd)
//}
//
//func RangeNumSamplesQuery(c, posBegin, posEnd string) (string, []string) {
//	return RangeQuery("%NS\n", c, posBegin, posEnd)
//}
//
//func RangeTotalAllelesQuery(c, posBegin, posEnd string) (string, []string) {
//	return RangeQuery("%AN\n", c, posBegin, posEnd)
//}
//
//func RangePopulationAFQuery(population, c, posBegin, posEnd string) (string, []string) {
//	return RangeQuery(fmt.Sprintf("%%%s_AF\n", population), c, posBegin, posEnd)
//}

func AllChrPositionsQuery(c, dataset string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		"%POS\t",
		GetChromosomeFilepath(c, dataset),
	}
}

//func GetSnpsAtPosition(c string, p string) ([]string, error) {
//	query, args := RangeSnpValuesQuery(c, p, p)
//	cmd := exec.Command(query, args...)
//	output, err := cmd.Output()
//	if err != nil {
//		return nil, err
//	}
//	samples := strings.Split(string(output), "\t")
//	return samples[:len(samples)-1], nil
//}

func NormalizeAllele(allele string) (string, error) {
	switch allele {
	case "0", "1":
		return allele, nil
	// simplify multiallelic cases
	case "2", "3", "4", "5", "6", "7", "8", "9":
		return "1", nil
	case ".":
		return ".", errors.New("unknown allele")
	default:
		return "", errors.New("invalid allele: " + allele)
	}
}

func NormalizeSnp(snp string) (string, error) {
	switch snp {
	case "0|0", "0|1", "1|0", "1|1":
		return snp, nil
	/*
		rare multi-allelic cases
	*/
	case "0|2", "0|3", "0|4", "0|5", "0|6", "0|7", "0|8", "0|9":
		return "0|1", nil
	case "2|0", "3|0", "4|0", "5|0", "6|0", "7|0", "8|0", "9|0":
		return "1|0", nil
	case "1|2", "2|1", "2|2",
		"1|3", "3|1", "2|3", "3|2", "3|3",
		"4|1", "1|4", "2|4", "4|2", "4|3", "3|4", "4|4",
		"5|1", "1|5", "2|5", "5|2", "5|3", "3|5", "4|5", "5|4", "5|5",
		"6|1", "1|6", "2|6", "6|2", "6|3", "3|6", "4|6", "6|4", "5|6", "6|5", "6|6",
		"7|1", "1|7", "2|7", "7|2", "7|3", "3|7", "4|7", "7|4", "5|7", "7|5", "6|7", "7|6", "7|7",
		"8|1", "1|8", "2|8", "8|2", "8|3", "3|8", "4|8", "8|4", "5|8", "8|5", "6|8", "8|6", "7|8", "8|7", "8|8",
		"9|1", "1|9", "2|9", "9|2", "9|3", "3|9", "4|9", "9|4", "5|9", "9|5", "6|9", "9|6", "7|9", "9|7", "9|8", "8|9", "9|9":
		return "1|1", nil
	default:
		return snp, errors.New("unknown snp value")
	}
}

func SnpToPair(snp string) ([]uint8, error) {
	switch snp {
	case "0|0":
		return []uint8{0, 0}, nil
	case "0|1":
		return []uint8{0, 1}, nil
	case "1|0":
		return []uint8{1, 0}, nil
	case "1|1":
		return []uint8{1, 1}, nil
	default:
		return nil, fmt.Errorf("invalid snp value: %s", snp)
	}
}

func SnpToSum(snp string) (uint8, error) {
	switch snp {
	case "0|0", "0,0":
		return 0, nil
	case "0|1", "1|0", "0,1", "1,0":
		return 1, nil
	case "1|1", "1,1":
		return 2, nil
	default:
		return 255, fmt.Errorf("invalid snp value: %s", snp)
	}
}

func LoadAncestry() map[string]string {
	data := make(map[string]string)
	file, err := os.Open("data/superpopulations.json")
	if err != nil {
		fmt.Println("Error opening populations file:", err)
		return data
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	if err = decoder.Decode(&data); err != nil {
		fmt.Println("Error decoding populations file:", err)
	}
	return data
}
