package tools

import (
	"errors"
	"fmt"
	"os/exec"
	"strings"
)

func GetChromosomeFilepath(chr string) string {
	path := "/gpfs/commons/datasets/1000genomes/hg38/"
	return path + "ALL.chr" + chr + ".phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz"
}

func IndividualSnpsQuery(c string, p string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		"[%SAMPLE=%GT\t]",
		"-r",
		fmt.Sprintf("%s:%s-%s", c, p, p),
		GetChromosomeFilepath(c),
	}
}

func RangeSnpValuesQuery(c string, posBegin, posEnd string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		"[%GT\t]\n",
		"-r",
		fmt.Sprintf("%s:%s-%s", c, posBegin, posEnd),
		GetChromosomeFilepath(c),
	}
}

func AllChrPositionsQuery(c string) (string, []string) {
	return "bcftools", []string{
		"query",
		"-f",
		"%POS\t",
		GetChromosomeFilepath(c),
	}
}

func GetSnpsAtPosition(c string, p string) ([]string, error) {
	query, args := RangeSnpValuesQuery(c, p, p)
	cmd := exec.Command(query, args...)
	output, err := cmd.Output()
	if err != nil {
		return nil, err
	}
	samples := strings.Split(string(output), "\t")
	return samples[:len(samples)-1], nil
}

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
		rare multiallelic cases
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

func SnpToSum(snp string) (float64, error) {
	switch snp {
	case "0|0":
		return 0, nil
	case "0|1", "1|0":
		return 1, nil
	case "1|1":
		return 2, nil
	default:
		return -1, fmt.Errorf("invalid snp value: %s", snp)
	}
}
