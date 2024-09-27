#!/bin/bash
#SBATCH --job-name=VCFft2
#SBATCH --error=VCFft2.err
#SBATCH --output=VCFft2.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G
#SBATCH --time=24:00:00

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate var_call

# Define region prefixes
regions=("CH" "UK" "US")
vcf_path_base="/work_beegfs/suaph296/Zymoproj"
ref="/work_beegfs/suaph296/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

for region in "${regions[@]}"; do
    vcf_path="${vcf_path_base}/${region}samples/gVCF"
    input_vcf="${vcf_path}/Zt_${region}_genotyped.vcf"
    output_vcf="${vcf_path}/Zt_${region}_genotyped.SNP.vcf"

    if [ ! -f "$input_vcf" ]; then
        echo "Genotyped VCF file for region $region does not exist. Skipping..."
        continue
    fi

    # Command for selecting only SNPs
    gatk --java-options "-Xmx400G" SelectVariants \
        -R $ref \
        -V $input_vcf \
        --select-type-to-include SNP \
        -O $output_vcf

    echo "Variant selection complete for region $region."
done