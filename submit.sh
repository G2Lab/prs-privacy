#!/bin/bash

num_chunks=$1
chunk_size=$((2504 / num_chunks + 1))

for i in $(seq 0 $((num_chunks - 1)))
do
    sbatch solve.slurm "$i" "$chunk_size"
    sleep 1
done
