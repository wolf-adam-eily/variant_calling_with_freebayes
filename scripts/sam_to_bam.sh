#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load samtools

for file in ../alignments/*; do export BASENAME=$(basename $file .sam).bam; samtools sort -@ 1 $file ../binary_alignments/$BASENAME;echo $BASENAME;done;

