#!/bin/bash
#SBATCH --job-name=solve                   # Job name
#SBATCH --mem=30G                         # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=solution_%j.log           # Standard output and error log
#SBATCH --cpus-per-task=2                  # Number of CPUs per task
#SBATCH --time=1-0:00:00                    # Time limit hrs:min:sec
#SBATCH --partition=pe2                    # Partition (queue) name
#SBATCH --priority=TOP                     # Job priority
#SBATCH --mail-user=knikitin@nygenome.org  # Email to which notifications will be sent

module load bcftools/1.17
module load htslib/1.17
module load go/1.20.2
go build -o solve ./eval/
GOMEMLIMIT=28GiB ./solve $1 $2