#!/bin/bash

for i in {9..12}
do
    sbatch solve.slurm "$i"
done