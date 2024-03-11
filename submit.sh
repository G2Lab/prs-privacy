#!/bin/bash

num_chunks="$1"
chunk_size=$((2504 / num_chunks))
#chunk_size=$((100 / num_chunks))

for i in $(seq 0 "$num_chunks")
do
    sbatch solve.slurm "$i" "$chunk_size"
    sleep 1
done
