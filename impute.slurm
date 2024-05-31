#!/bin/bash
#SBATCH --job-name=impute                   # Job name
#SBATCH --mem=15G                         # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=solution_%j.log           # Standard output and error log
#SBATCH --cpus-per-task=2                  # Number of CPUs per task
#SBATCH --time=2-0:00:00                    # Time limit hrs:min:sec
#SBATCH --partition=pe2                    # Partition (queue) name
#SBATCH --priority=TOP                     # Job priority
#SBATCH --mail-user=knikitin@nygenome.org  # Email to which notifications will be sent

module load bcftools/1.20
module load minimac4/4.1.2
module load go/1.20.2
go build -o solve ./eval/
./solve $1 $2