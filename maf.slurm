#!/bin/bash
#SBATCH --job-name=priors                  # Job name
#SBATCH --mem=64G                          # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=prior_%j.log             # Standard output and error log
#SBATCH --cpus-per-task=22                  # Number of CPUs per task
#SBATCH --time=1-00:00:00                    # Time limit hrs:min:sec
#SBATCH --partition=pe2                    # Partition (queue) name
#SBATCH --mail-user=knikitin@nygenome.org  # Email to which notifications will be sent

module load bcftools/1.17
module load htslib/1.17
module load go/1.20.2
go build -o ./data/calculate ./data/
./data/calculate --maf=true