#!/bin/bash
#SBATCH --job-name=fastq_dump
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --mem=10G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
module load sratoolkit

fastq-dump SRR8068137
fastq-dump SRR8068138
fastq-dump SRR8068139
fastq-dump SRR8068140
fastq-dump SRR8068141
fastq-dump SRR8068142
fastq-dump SRR8068144
fastq-dump SRR8068147
fastq-dump SRR8068148
