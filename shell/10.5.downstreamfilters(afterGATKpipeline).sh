#!/bin/bash
#SBATCH --job-name=within_field_filter
#SBATCH --error=within_field_filter.err
#SBATCH --output=within_field_filter.out
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=08:00:00

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Set the folders containing the vcf files
folders=("~/Zymoproj/UKsamples/VCF" "~/Zymoproj/USsamples/VCF" "~/Zymoproj/CHsamples/VCF")

# Loop through each folder
for folder in "${folders[@]}"; do
  # Loop through each vcf file in the folder for the strict set
  for file in "$folder"/*qualityfilter2021.excl.vcf; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"
    filename_base="${filename_base#Zt_}"

    # VCFtools processing steps 
    vcftools --vcf "$file" \
      --recode --recode-INFO-all \
      --max-missing 0.8 \
      --mac 1 \
      --remove-filtered-all --remove-indels \
      --min-alleles 2 --max-alleles 2 \
      --out "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1"

    # Compress the VCF file
    bgzip "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf"
    tabix -p vcf "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf.gz"

    # Subset chr 1-13 (core)
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13 "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf.gz" -O z -o "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.core.vcf.gz"
    tabix -p vcf "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.core.vcf.gz"

    # Subset chr 14-21 (accessory)
    bcftools view -r 14,15,16,17,18,19,20,21 "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf.gz" -O z -o "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.acc.vcf.gz"
    tabix -p vcf "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.acc.vcf.gz"

    # Convert VCF to PLINK format and generate plink files for the whole input file
    plink --vcf "${folder}/Zt_${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf.gz" \
      --make-bed --freq --missing \
      --out "${folder}/Zt_${filename_base}.strict"
  done

  # Loop through each vcf file in the folder for the relaxed set
  for file in "$folder"/*qualityfilter.excl.vcf; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"
    filename_base="${filename_base#Zt_}"

    # VCFtools processing steps for relaxed set
    vcftools --vcf "$file" \
      --recode --recode-INFO-all \
      --mac 1 \
      --out "${folder}/Zt_${filename_base}.relaxed.mac1"

    # Compress the VCF file
    bgzip "${folder}/Zt_${filename_base}.relaxed.mac1.recode.vcf"
    tabix -p vcf "${folder}/Zt_${filename_base}.relaxed.mac1.recode.vcf.gz"

    # Subset chr 1-13 (core)
    bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13 "${folder}/Zt_${filename_base}.relaxed.mac1.recode.vcf.gz" -O z -o "${folder}/Zt_${filename_base}.relaxed.core.vcf.gz"
    tabix -p vcf "${folder}/Zt_${filename_base}.relaxed.core.vcf.gz"

    # Subset chr 14-21 (accessory)
    bcftools view -r 14,15,16,17,18,19,20,21 "${folder}/Zt_${filename_base}.relaxed.mac1.recode.vcf.gz" -O z -o "${folder}/Zt_${filename_base}.relaxed.acc.vcf.gz"
    tabix -p vcf "${folder}/Zt_${filename_base}.relaxed.acc.vcf.gz"

    # Convert VCF to PLINK format and generate plink files for the whole input file
    plink --vcf "${folder}/Zt_${filename_base}.relaxed.mac1.recode.vcf.gz" \
      --make-bed --freq --missing \
      --out "${folder}/Zt_${filename_base}.relaxed"
  done
done