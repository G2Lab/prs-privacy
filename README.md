# Private Information Leakage from Polygenic Risk Scores

This repository is a companion to the paper "Private Information Leakage from Polygenic Risk Scores".
We provide the code and the data for running a genotype-recovery demo, for plotting the figures from the paper, and for end-to-end reproducibility of the experiments.

#### Table of contents:
- [Setup](#setup)
  - [Installing dependencies](#installing-dependencies)
  - [Downloading pre-processed 1000 Genomes data](#downloading-pre-processed-1000-genomes-data)
- [Genotype recovery demo](#genotype-recovery-demo)
- [Processing the results from the paper](#processing-the-results-from-the-paper)
    - [Re-plotting the figures](#re-plotting-the-figures)
- [End-to-end reproducibility](#end-to-end-reproducibility)

## Setup
### Installing dependencies

The experiments require Go of version 1.23.1 or above. You can manually install the latest version for your operating system and then run `go get -u ./...` to download the dependencies. Alternatively, run the command below to automatically install it. The automatic script currently supports only MacOS and Linux.
```
make install-go
```
Plotting the figures requires Python and the packages in `requirements.txt`. To automatically configure your setup, run
```
make install-python
```

### Downloading pre-processed 1000 Genomes data

We have extracted the genotypes, allele frequencies and genotype frequencies from the samples of the 1000 Genomes Project for all the variants that we use in our experiments. and have uploaded it, along with relevant PGS metadata, to a separate [repository](https://github.com/G2Lab/prs-privacy-data). All the 1000 Genomes experiments are possible to run without access to the full dataset by using the pre-processed data.
To download it, run the command below.
```
make download-data
```
The archive is ~19 MiB, the script will unpack the data and place the following files into the inputs folder:
- `PGS00XXXX_hmPOS_GRCh37.txt` are the PRS metadata files from the [PGS Catalog](https://www.pgscatalog.org/).
- `PGS00XXXX.json` stores the genotypes for each corresponding variant and the PRS for each individual in the dataset. This is loaded in the code as `cohort`.
- `PGS00XXXX.efal` stores whether the effect allele is the reference allele (0) or the alternative allele (1).
- `PGS00XXXX.stat` stores the allele and genotype frequencies of each dataset ancestry for each SNP in the PRS.
- `PGS00XXXX.scores` stores just the corresponding PRS for each individual in the dataset.  

We _do not_ upload any UK Biobank data, including any statistics, as access to it requires a license. All the experiments that involve the UK Biobank data require a link to the dataset to be reproduced.


## Genotype recovery demo
The demo shows the process of genotype recovery from a set of PRSs for a random individual from the 1000 Genomes dataset.
The demo uses the same PGS selection as the "Patient genotypes can be recovered from publicly available PRSs" section of the paper, but we limit the size of the PGSs to be considered to at most 30 SNPs for efficiency. 
The demo takes approximately a minute to run.
It prints the recovered and true genotypes for each PRS and reports the overall
recovery accuracy at the end.
Run it with:

```
make demo
```

## Processing the results from the paper
We have uploaded the measurements from the paper to the same secondary [repository](https://github.com/G2Lab/prs-privacy-data). To download and uncompress the results, run the command below. The archive is ~200 MiB.
```
make download-results
```
The `results` folder will contain the following:
- `validated_pgs.json` contains the full selection of the PGSs from our experiments with the corresponding size (the number of SNPs) for each.
- `validated_loci.json` contains all the covered loci with a list of PGSs where each locus appears.
- `recoveryOutput/` contains genotype predictions and `recoveryAccuracy/` contains the accuracies of these predictions for the 2535 individuals of the 1000 Genomes dataset. See the "Patient genotypes can be recovered from publicly available PRSs" section and Figure 3 in the paper.
- `linking/` contains the PLINK output from running the KING-robust and the GCTA algorithms for matching the predicted genotypes against the true values. See the "Patients and their relatives can be retrieved from genealogy databases using PRSs" section, Figure 4 and Supplementary Figure 2 in the paper. 
-  `uniqueness/` contains the analytical estimates and dataset measurements for the uniqueness and anonymity guarantees of PRS values. Only the results for the 1000 Genomes dataset are included due to the restrictions of sharing the UK Biobank data. See the "Anonymized genotype databases can be de-anonymized by using a single PRS from known
   patients without the need for genotype prediction" section and Figures 5B-C in the paper.  

### Re-plotting the figures

To reproduce the figures from the paper, run
- `make plot-recovery` plots the genotype recovery accuracy results.
- `make plot-linking` plots the linking results.
- `make plot-uniqueness` plots the score uniqueness and anonymity results. It partially requires access to the UK Biobank data.
- `make plot-rounding` plots the PRS distribution when rounding the effect weights to different decimal-places precision. It requires access to the UK Biobank data. 
- `make plot-all` plots all the above.

The figures will be saved in `analysis/figures`.

## End-to-end reproducibility

To reproduce the experiments end-to-end: 
1. Add the paths to the downloaded 1000 Genomes and UK Biobank data to `datasets.yaml`.
2. Run `make pgs-selection` to retrieve suitable PGSs from the PGS Catalog and form the selection for `validated_pgs.json` and `validated_loci.json`.
3. Perform genotype recovery for all the samples in the 1000 Genomes dataset by running `go run eval/run.go -e=solve CHUNK_NUM CHUNK_SIZE` where `CHUNK_NUM` is a sequential id of the chunk and `CHUNK_SIZE` is the number of samples in the chunk. The program sorts all the samples alphabetically. The start solving position is `CHUNK_NUM * CHUNK_SIZE`. Resolving 298 PRSs for up to 50 SNPs requires 2h40m of computing time and 127 GB of RAM on average per sample. The program relies on bcftools to retrieve the genotype data.
4. Run `go run eval/run.go -e=genfreq` to save genotype frequencies for all the predicted SNPs.
5. Calculate the genotype-recovery accuracy by running `go run eval/run.go -e=accuracy`.
6. Calculate the relatedness results by running `go run eval/run.go -e=linking`. You will need PLINK installed.
7. Calculate the uniqueness and anonymity results by running `go run eval/run.go -e=uniqueness_gg` for the 1000 Genomes data and `go run eval/run.go -e=uniqueness_uk` for the UK Biobank data. The analytical estimation is computationally intensive.
8. Finally, run `go run eval/run.go -e=rounding` to measure the impact of rounding effect weights on the PRS distribution.
9. Proceed with producing figures as detailed in [Re-plotting the figures](#re-plotting-the-figures).