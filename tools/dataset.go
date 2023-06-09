package tools

import "fmt"

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

func SnpToValue(allele string) (float64, error) {
	switch allele {
	case "0|0":
		return 0, nil
	case "0|1", "1|0":
		return 1, nil
	case "1|1":
		return 2, nil
	case "0|2", "2|0", "1|2", "2|1", "2|2", "0|3", "3|0", "1|3", "3|1", "3|2", "2|3", "3|3":
		return 3, nil
	default:
		return 0, fmt.Errorf("invalid allele value: %s", allele)
	}
}
