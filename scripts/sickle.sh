#!/bin/bash
#SBATCH --job-name=sickle
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
module load sickle

for file in ../sample_files/*; do export BASENAME=$(basename $file .fastq).trimmed.fastq; sickle se -f $file -t sanger -o ../trimmed_sample_files/$BASENAME -q 30 -l 50; done;
