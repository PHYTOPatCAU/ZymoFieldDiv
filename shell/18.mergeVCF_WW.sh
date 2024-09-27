#!/bin/bash
#SBATCH --job-name=merge_vcf_ww
#SBATCH --error=merge_vcf_ww.err
#SBATCH --output=merge_vcf_ww.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Directory containing the VCF files
merged_directory="/work_beegfs/suaph296/Zymoproj/merged"

# Initialize an array to store the VCF files
vcf_files=()

# Loop through each VCF file in the merged directory
for vcf in "${merged_directory}"/*relaxed.mac1.recode.vcf.gz; do
  vcf_files+=("$vcf")
done

# Check if there are any VCF files to merge
if [ ${#vcf_files[@]} -gt 0 ]; then
  output_vcf="${merged_directory}/Zt_WW_relaxed_mac1.vcf.gz"
  bcftools concat "${vcf_files[@]}" -Oz -o "$output_vcf"
  echo "Merged VCF file saved to ${output_vcf}"

  # Index the merged VCF file
  bcftools index "$output_vcf"
  echo "Index for merged VCF file created at ${output_vcf}.csi"

  # Subset the core VCF file
  core_vcf="${merged_directory}/Zt_WW_relaxed_core_v2.vcf.gz"
  bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13 "$output_vcf" -O z -o "$core_vcf"
  echo "Core VCF file saved to ${core_vcf}"

  # Index the core VCF file
  tabix -p vcf "$core_vcf"
  echo "Index for core VCF file created at ${core_vcf}.tbi"
else
  echo "No VCF files found to merge."
fi