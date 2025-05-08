#!/bin/bash
#SBATCH --job-name=checkdepth.sh
#SBATCH --time=01:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Set the folders containing the vcf files
folders=("~/Zymoproj/UKsamples/VCF" "~/Zymoproj/USsamples/VCF" "~/Zymoproj/CHsamples/VCF")

# Loop through each folder
for folder in "${folders[@]}"; do
  # Loop through each strict vcf file in the folder
  for file in "$folder"/*max-m-80.biallelic-only.mac1.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # vcftools check of depth per sample
    vcftools --gzvcf "$file" --depth --out "${folder}/${filename_base}_depth"

    # Calculate the number of SNPs per sample
    vcftools --gzvcf "$file" --counts --out "${folder}/${filename_base}_counts"
  done

  # Loop through each relaxed vcf file in the folder
  for file in "$folder"/*relaxed.mac1.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # vcftools check of depth per sample
    vcftools --gzvcf "$file" --depth --out "${folder}/${filename_base}_depth"

    # Calculate the number of SNPs per sample
    vcftools --gzvcf "$file" --counts --out "${folder}/${filename_base}_counts"
  done
done