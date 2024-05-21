#!/bin/bash
#SBATCH --job-name=compression                   # Job name
#SBATCH --mem=20G                          # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=comp_%j.log           # Standard output and error log
#SBATCH --cpus-per-task=4                  # Number of CPUs per task
#SBATCH --time=1-0:00:00                    # Time limit hrs:min:sec
#SBATCH --partition=pe2                    # Partition (queue) name
#SBATCH --priority=TOP                     # Job priority
#SBATCH --mail-user=knikitin@nygenome.org  # Email to which notifications will be sent

bash data/compress_references.sh