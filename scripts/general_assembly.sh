#!/bin/bash
#SBATCH --job-name=assembly
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

mkdir ../assembly

module load SPAdes

cd ../trimmed_sample_files/
 
spades.py --s1 SRR8068137.trimmed.fastq --s2 SRR8068139.trimmed.fastq --s3 SRR8068141.trimmed.fastq --s4 SRR8068144.trimmed.fastq --s5 SRR8068148.trimmed.fastq --s6 SRR8068138.trimmed.fastq --s7 SRR8068140.trimmed.fastq --s8 SRR8068142.trimmed.fastq --s9 SRR8068147.trimmed.fastq -o ../assembly -m 100 -t 8

