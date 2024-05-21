#!/bin/bash
module load minimac4
for i in {1..22}
do
    minimac4 --compress-reference /gpfs/commons/datasets/1000genomes/GRCh37/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > /gpfs/commons/groups/gursoy_lab/knikitin/prs/data/references/"$i".msav
    echo "Chr $i done"
done