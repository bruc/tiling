#!/bin/sh -x

for fa in ../references/{pAV-CMV-GFP,scAAV_CBA_eGFP}.fasta
do
    base=$(basename $fa .fasta)
    ./reformat_viral_plasmid.pl -itr_file=../references/itr-flip.fa \
				-file=$fa \
				>../references/${base}_split.fasta
done


