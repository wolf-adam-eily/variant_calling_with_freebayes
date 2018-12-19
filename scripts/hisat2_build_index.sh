#!/bin/bash
#SBATCH --job-name=hisat2_build_index
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

module load hisat2

hisat2-build -p 2 ../assembly/scaffolds.fasta ../index/whitefish
