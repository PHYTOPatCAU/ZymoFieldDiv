#!/bin/bash
#SBATCH --job-name=subset_bcf
#SBATCH --error=subset_bcf.err
#SBATCH --output=subset_bcf.out
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --time=03:00:00
# Load modules
module load gcc12-env
module load miniconda3

# Activate conda environment
conda activate samtools

# Define region prefixes
regions=("CH" "UK" "US")
vcf_path_base="~/Zymoproj"

for region in "${regions[@]}"; do
    vcf_path="${vcf_path_base}/${region}samples/VCF"
    input_vcf="${vcf_path}/Zt_${region}_genotyped.SNP.qualityfilter2021.excl.vcf"
    sample_names_file="${vcf_path}/${region}_sample_names.txt"
    subset_counts_file="${vcf_path}/subset_counts_${region}.txt"

    cd $vcf_path

    # Extract sample names
    bcftools query -l $input_vcf > $sample_names_file

    # Create output file with header
    printf "Subset Size\tSNP Count\n" > $subset_counts_file

    # Loop for different subset sizes
    for i in {10..160..10}; do
        # Randomly select samples
        shuf -n $((i + 10)) $sample_names_file | head -n $i > ${region}_subset_$i.txt

        # Count the SNPs in the subsets using vcftools
        vcftools \
            --vcf $input_vcf \
            --keep ${region}_subset_$i.txt \
            --mac 1 \
            --recode-INFO-all \
            --remove-filtered-all \
            --out ${region}_subset_$i

        # Append the subset size and SNP count to the output file
        snp_count=$(grep -v "^#" ${region}_subset_$i.recode.vcf | wc -l)
        printf "%d\t%d\n" $i $snp_count >> $subset_counts_file
    done

    echo "Subsampling complete for region $region."
done