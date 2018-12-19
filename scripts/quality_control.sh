#!/bin/bash
#SBATCH --job-name=quality_control
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

module load fastqc
module load MultiQC

mkdir /home/CAM/$USER/tmp/
export TMPDIR=/home/CAM/$USER/tmp/

for file in ../trimmed_sample_files/*; do fastqc $file -o ../quality_control; done;

cd ../quality_control
multiqc .
