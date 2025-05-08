#!/bin/bash
#SBATCH --job-name=concat_eff_per_field
#SBATCH --error=concat_eff_per_field.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Define the expected number of samples for each region
declare -A expected_samples
expected_samples=( ["UK"]=160 ["CH"]=158 ["US"]=97 )

# Base directory containing the region directories
base_directory="~/Zymoproj/merged"

# Regions to process
regions=("UK" "US" "CH")

# Loop over each region
for region in "${regions[@]}"; do
  # Directory containing the VCF files for the current region
  region_directory="${base_directory}/${region}_Effectors"

  # Initialize an array to store valid VCF files
  valid_vcfs=()

  # Loop through each VCF file in the region directory
  for vcf in "${region_directory}"/*.vcf.gz; do
    # Get the number of samples in the VCF file
    num_samples=$(bcftools query -l "$vcf" | wc -l)
    
    # Check if the number of samples matches the expected number
    if [ "$num_samples" -eq "${expected_samples[$region]}" ]; then
      valid_vcfs+=("$vcf")
    else
      echo "Skipping $vcf due to different number of samples ($num_samples)"
    fi
  done

  # Concatenate the valid VCF files
  if [ ${#valid_vcfs[@]} -gt 0 ]; then
    output_vcf="${base_directory}/${region}_eff.vcf.gz"
    bcftools concat "${valid_vcfs[@]}" -Oz -o "$output_vcf"
    echo "Concatenated VCF file for ${region} saved to ${output_vcf}"
  else
    echo "No valid VCF files to concatenate for ${region}."
  fi
done