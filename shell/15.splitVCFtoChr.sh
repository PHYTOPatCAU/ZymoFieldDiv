#!/bin/bash
#SBATCH --job-name=split_vcf_to_chr
#SBATCH --error=split_vcf_to_chr.err
#SBATCH --output=split_vcf_to_chr.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Set the folders containing the VCF files
folders=("~/Zymoproj/UKsamples/VCF" "~/Zymoproj/USsamples/VCF" "~/Zymoproj/CHsamples/VCF")

# Output directory
output_dir="~/Zymoproj/merged"

# Loop through each folder
for folder in "${folders[@]}"; do
  # Loop through each relaxed VCF file in the folder
  for file in "$folder"/*relaxed.mac1.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # Sort the VCF file by position
    sorted_vcf="${folder}/${filename_base}_sorted.vcf.gz"
    bcftools sort "$file" -Oz -o "$sorted_vcf"

    # Split the sorted VCF file by chromosome
    for i in {1..21}; do
      output_vcf="${output_dir}/${filename_base}_chr${i}_mac1.vcf.gz"
      bcftools view -r $i "$sorted_vcf" -Oz -o "$output_vcf"
    done
  done
done