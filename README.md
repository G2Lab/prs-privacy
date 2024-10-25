# Private Information Leakage from Polygenic Risk Scores

This repository contains the code and data for the paper "Private Information Leakage from Polygenic Risk Scores".

### Setup

#### Installing dependencies

Run the following commands to install 
1. Go and the required packages 
2. Python and the packages for making plots

```
make install-go
make install-python
```

#### Downloading pre-processed 1000 Genomes data

The pre-processed data of the 1000 Genomes Project individuals for the experiments 
and the demo are available for download by running the following command:
```
make download
```

### Running a demo
The demo randomly selects an individual from the 1000 Genomes Project dataset and 
recovers their genotypes from the PRSs that use up to 30 variants.
The demo prints the recovered and true genotypes for each PRS and reports the overall
recovery accuracy at the end. The demo takes approximately a minute to run.

```
make demo
```

### Reproducing the experiments