#!/bin/bash
#SBATCH --job-name=excludeFiltered
#SBATCH --error=excludeFiltered.err
#SBATCH --output=excludeFiltered.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=04:00:00
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
    vcf_path="${vcf_path_base}/${region}samples/VCF"
    input_vcf="${vcf_path}/Zt_${region}_genotyped.SNP.qualityfilter2021.vcf"
    output_vcf="${vcf_path}/Zt_${region}_genotyped.SNP.qualityfilter2021.excl.vcf"

    if [ ! -f "$input_vcf" ]; then
        echo "Quality-filtered VCF file for region $region does not exist. Skipping..."
        continue
    fi

    # Command for excluding filtered variants
    gatk SelectVariants \
        -R $ref \
        -V $input_vcf \
        -O $output_vcf \
        --exclude-filtered true

    echo "Excluding filtered variants complete for region $region."
done