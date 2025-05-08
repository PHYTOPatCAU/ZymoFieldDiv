#!/bin/bash
#SBATCH --job-name=admixww_relaxed
#SBATCH --error=admixww_relaxed.err
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --qos=long
#SBATCH --time=5-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Set PATH
vcf_path=~/Zymoproj/merged
output_dir=${vcf_path}/admixture/relaxed

# make output directory
mkdir -p $output_dir

# Change to VCF directory
cd $vcf_path

# Pruning VCF file and making BED file for admixture analysis
plink --vcf Zt_WW_relaxed_mac1.vcf.gz --maf 0.05 --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.2 --make-bed --out ${output_dir}/maf05prunning

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
plink --vcf Zt_WW_relaxed_mac1.vcf.gz --maf 0.05 --double-id --allow-extra-chr --set-missing-var-ids @:# --extract ${output_dir}/maf05prunning.prune.in --make-bed --out ${output_dir}/maf05pruned_data

# Modify the .bim file to replace chromosome names with 0
awk '{$1="0";print $0}' ${output_dir}/maf05pruned_data.bim > ${output_dir}/maf05pruned_data.bim.tmp
mv ${output_dir}/maf05pruned_data.bim.tmp ${output_dir}/maf05pruned_data.bim

# Run Admixture
for K in {1..10}; do
    for repeat in {1..10}; do
        admixture --cv ${output_dir}/maf05pruned_data.bed $K -j24 -s ${repeat} | tee ${output_dir}/log${K}.${repeat}.out
    done
done
# Extract CV errors and save to a text file
grep -h "CV error" ${output_dir}/log*.out | awk '{print FILENAME, $3}' > ${output_dir}/cv_errors.txt
