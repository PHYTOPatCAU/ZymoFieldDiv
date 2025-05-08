#!/bin/bash
#SBATCH --job-name=variantfiltev
#SBATCH --error=variantfitlev.err
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=12:00:00

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate var_call

# Define region prefixes
regions=("CH" "UK" "US")
vcf_path_base="~/Zymoproj"
ref="~/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

for region in "${regions[@]}"; do
    vcf_path="${vcf_path_base}/${region}samples/gVCF"
    input_vcf="${vcf_path}/Zt_${region}_genotyped.SNP.vcf"
    output_dir="${vcf_path_base}/${region}samples/VCF"
    output_vcf="${output_dir}/Zt_${region}_genotyped.SNP.qualityfilter.vcf"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    if [ ! -f "$input_vcf" ]; then
        echo "Genotyped SNP VCF file for region $region does not exist. Skipping..."
        continue
    fi

    # Command for variant filtration
    gatk --java-options "-Xmx200G" VariantFiltration \
        -R $ref \
        -V $input_vcf \
        --filter-name "MQ" --filter-expression "MQ < 20.0" \
        --filter-name "QDFilter" --filter-expression "QD < 5.0" \
        --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
        --filter-name "SOR" --filter-expression "SOR > 3.0" \
        --filter-name "FS" --filter-expression "FS > 60.0" \
        --filter-name "Low_depth3" --filter-expression "DP < 10" \
        --filter-name "ReadPosRankSum_lower" --filter-expression "ReadPosRankSum < -2.0" \
        --filter-name "ReadPosRankSum_upper" --filter-expression "ReadPosRankSum > 2.0" \
        --filter-name "MQRankSum_lower" --filter-expression "MQRankSum < -2.0" \
        --filter-name "MQRankSum_upper" --filter-expression "MQRankSum > 2.0" \
        -O $output_vcf

    echo "Variant filtration complete for region $region."
done