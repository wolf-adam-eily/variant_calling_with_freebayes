#!/bin/bash
#SBATCH --job-name=freebayes
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

module load freebayes
module load vcflib

freebayes --fasta-reference ../assembly/K21/final_contigs.fasta \
-b ../alignments/SRR8068137.trimmed.bam \
-b ../alignments/SRR8068139.trimmed.bam \
-b ../alignments/SRR8068141.trimmed.bam \
-b ../alignments/SRR8068144.trimmed.bam \
-b ../alignments/SRR8068148.trimmed.bam \
-b ../alignments/SRR8068138.trimmed.bam \
-b ../alignments/SRR8068140.trimmed.bam \
-b ../alignments/SRR8068142.trimmed.bam \
-b ../alignments/SRR8068147.trimmed.bam | vcffilter -f "QUAL > 20" >../variant_call/results.vcf

