#!/bin/bash
#SBATCH --job-name=hisat2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load hisat2

for file in ../trimmed_sample_files/*; do export BASENAME=$(basename $file .fastq).sam; hisat2 -p 2 -x ../index/whitefish -q $file -S ../alignments/$BASENAME; done;

