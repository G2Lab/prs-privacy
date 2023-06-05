#!/usr/bin/env python3

from collections import defaultdict
import matplotlib.pyplot as plt
import subprocess
from pgs import PGS


def allele_to_value(allele):
    if allele == "0|0":
        return 0
    elif allele == "0|1" or allele == "1|0":
        return 1
    elif allele == "1|1":
        return 2
    else:
        raise ValueError("Invalid allele value: " + allele)


def get_chromosome_filepath(chr):
    path = "/gpfs/commons/datasets/1000genomes/hg38/"
    return path + "ALL.chr" + chr + ".phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz"


def bcftools_query(c, p):
    return [
        "bcftools",
        "query",
        "-f",
        "[%SAMPLE=%GT\t]",
        "-r",
        c + ":" + str(p) + "-" + str(p),
        get_chromosome_filepath(c),
    ]


if __name__ == "__main__":
    pgs = PGS()
    pgs.load("PGS000073_hmPOS_GRCh38.txt")

    individuals = defaultdict(lambda: defaultdict(int))
    for variant in pgs.variants.values():
        chr, position = variant.get_hm_chr(), variant.get_hm_pos()
        try:
            output = subprocess.check_output(bcftools_query(chr, position), universal_newlines=True)
            samples = output.split("\t")
            for sample in samples:
                try:
                    individual, genotype = sample.split("=")
                except ValueError:
                    print("Error splitting sample:", sample)
                value = allele_to_value(genotype)
                individuals[individual][f'{chr}:{position}'] = value
                individuals[individual]['total'] += value * variant.get_weight()
        except subprocess.CalledProcessError as e:
            # Handle any errors that occurred during command execution
            print("Error executing bcftools command:", e)

    sorted_ind = sorted(individuals.keys(), key=lambda k: individuals[k]['total'])
    for ind in sorted_ind:
        print(ind, individuals[ind]['total'])

    # for individual, values in individuals.items():
    #     print(individual, values['total'])

    # # Plotting score distributions
    # scores = [individual['total'] for individual in individuals.values()]
    # plt.hist(scores, bins=100, color='orange')  # Customize the number of bins as needed
    # plt.xlabel('Scores')
    # plt.ylabel('Frequency')
    # plt.title(f'Polygenic-score distribution for {pgs.trait_name}')
    # plt.savefig('score_distribution.png', dpi=300)

    # # Find individual alleles for one score
    # print(find_most_likely_genotype(individuals['HG01438']['total'],
    #                                 [variant.get_weight() for variant in pgs.variants.values()]
    #                                 [1 / pgs.variants_number] * pgs.variants_number))
