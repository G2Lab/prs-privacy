#!/bin/bash

num_chunks="$1"
chunk_size=$((2504 / num_chunks))
#chunk_size=$((50 / num_chunks))

for ((i = 0; i < num_chunks; i++))
#for ((i = 2; i < 4; i++))
do
    sbatch solve.slurm "$i" "$chunk_size"
    sleep 3
done
