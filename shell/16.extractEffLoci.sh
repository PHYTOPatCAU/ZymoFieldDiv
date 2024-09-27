#!/bin/bash
#SBATCH --job-name=extract_eff_loci
#SBATCH --error=extract_eff_loci.err
#SBATCH --output=extract_eff_loci.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Define the loci file
LOCI_FILE="effector_loci2.csv"

# Set the folders containing the VCF files
folders=("/work_beegfs/suaph296/Zymoproj/UKsamples/VCF" "/work_beegfs/suaph296/Zymoproj/USsamples/VCF" "/work_beegfs/suaph296/Zymoproj/CHsamples/VCF")

# Loop through each folder
for folder in "${folders[@]}"; do
  # Determine the region from the folder name
  region=$(basename "$folder" | cut -d'/' -f3 | cut -d's' -f1)

  # Create the output directory
  output_dir="/work_beegfs/suaph296/Zymoproj/merged/${region}_Effectors"
  mkdir -p "$output_dir"

  # Loop through each relaxed VCF file in the folder
  for VCF_FILE in "$folder"/*relaxed.mac1.recode.vcf.gz; do
    # Extract the first 5 characters of the VCF file name
    DIR_NAME=$(basename "$VCF_FILE" | cut -c1-5)

    # Loop through each line in effector_loci.csv
    while IFS=';' read -r gene loc; do
      # Skip the header line
      if [ "$gene" == "gene" ]; then
        continue
      fi

      # Debugging: Print the values read from the file
      echo "Read from file - Gene: $gene, Location: $loc"

      # Check if the values are not empty
      if [ -z "$gene" ] || [ -z "$loc" ]; then
        echo "Error: One or more fields are empty. Skipping this line."
        continue
      fi

      # Ensure there are no hidden characters in the region string
      loc=$(echo "$loc" | tr -d '\r')

      # Debugging: Print the region to check its format
      echo "Processing region: $loc for gene: $gene"

      # Print the exact bcftools command being run
      echo "Running command: bcftools view -r $loc \"$VCF_FILE\" -Oz -o \"${output_dir}/${gene}.vcf.gz\""

      # Use bcftools to extract the region for each gene
      bcftools view -r "$loc" "$VCF_FILE" -Oz -o "${output_dir}/${gene}.vcf.gz"

      # Check if bcftools command was successful
      if [ $? -ne 0 ]; then
        echo "Error processing region: $loc for gene: $gene"
      else
        echo "Successfully processed region: $loc for gene: $gene"
      fi
    done < <(tail -n +2 "$LOCI_FILE") # Skip the header line
  done
done