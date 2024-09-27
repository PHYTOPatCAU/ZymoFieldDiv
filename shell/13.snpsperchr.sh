#!/bin/bash
#SBATCH --job-name=snpsperchr_all
#SBATCH --error=snpsperchr_all.err
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
folders=("/work_beegfs/suaph296/Zymoproj/UKsamples/VCF" "/work_beegfs/suaph296/Zymoproj/USsamples/VCF" "/work_beegfs/suaph296/Zymoproj/CHsamples/VCF")

# Loop through each folder
for folder in "${folders[@]}"; do
  # Loop through each strict vcf file in the folder
  for file in "$folder"/*max-m-80.biallelic-only.mac1.maf05.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # Get the SNP count per chromosome and add file type
    output_file="${folder}/snps_per_chr_${filename_base}.txt"
    bcftools query -f '%CHROM\n' "$file" | sort | uniq -c | awk -v fname="$filename_no_ext" -v ftype="strict" '{print fname"\t"$2"\t"$1"\t"ftype}' > "$output_file"
  done

  # Loop through each relaxed vcf file in the folder
  for file in "$folder"/*relaxed.mac1.recode.vcf.gz; do
    # Extract filename without path and extension
    filename_no_ext="${file##*/}"
    filename_base="${filename_no_ext%%.*}"

    # Get the SNP count per chromosome and add file type
    output_file="${folder}/snps_per_chr_${filename_base}.txt"
    bcftools query -f '%CHROM\n' "$file" | sort | uniq -c | awk -v fname="$filename_no_ext" -v ftype="relaxed" '{print fname"\t"$2"\t"$1"\t"ftype}' > "$output_file"
  done
done

# Get a single txt file in /work_beegfs/suaph296/Zymoproj directory with all SNP counts
cat /work_beegfs/suaph296/Zymoproj/*samples/VCF/snps_per_chr_*.txt > /work_beegfs/suaph296/Zymoproj/snp_counts_all.txt