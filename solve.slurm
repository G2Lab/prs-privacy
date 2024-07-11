#!/bin/bash
#SBATCH --job-name=prs                   # Job name
#SBATCH --mem=250G                         # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=solution_%j.log           # Standard output and error log
#SBATCH --cpus-per-task=4                  # Number of CPUs per task
#SBATCH --time=2-0:00:00                    # Time limit hrs:min:sec
#SBATCH --partition=pe2                    # Partition (queue) name
#SBATCH --priority=4294967293                     # Job priority
#SBATCH --mail-user=knikitin@nygenome.org  # Email to which notifications will be sent

module load bcftools/1.20
module load go/1.20.2
go build -o solve ./eval/
GOMEMLIMIT=240GiB ./solve $1 $2