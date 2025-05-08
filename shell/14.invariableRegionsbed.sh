#!/bin/bash
#SBATCH --job-name=invariable_regions_bed
#SBATCH --error=invariable_regions_bed.err
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

# Set the folders containing the vcf files
folders=("~/Zymoproj/UKsamples/VCF" "~/Zymoproj/USsamples/VCF" "~/Zymoproj/CHsamples/VCF")

# Loop through each folder
for folder in "${folders[@]}"; do
  # Loop through each relaxed vcf file in the folder
  for file in "$folder"/*relaxed.mac1.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # Generate BED file for invariant regions
    output_bed="${folder}/${filename_base}_invariable_regions.bed"
    bcftools view -i 'TYPE="snp" && INFO/AC=0' "$file" | \
    bcftools query -f '%CHROM\t%POS0\t%END\n' > "$output_bed"
  done
done