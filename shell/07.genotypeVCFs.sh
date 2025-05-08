#!/bin/bash

# Region prefixes
regions=("CH" "UK" "US")
vcf_path_base="~/Zymoproj"
ref="~/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

for region in "${regions[@]}"; do
    # Generate a temporary script for the region
    script_path="/tmp/genotype_${region}.sh"
    cat > "$script_path" << EOF
#!/bin/bash
#SBATCH --job-name=genotype_${region}
#SBATCH --error=genotype_${region}.err
#SBATCH --output=genotype_${region}.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --qos=long
#SBATCH --time=12-00:00:00
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

vcf_path="${vcf_path_base}/${region}samples/gVCF"
combined_vcf="\${vcf_path}/Zt_${region}.g.vcf"
output_vcf="\${vcf_path}/Zt_${region}_genotyped.vcf"

if [ ! -f "\$combined_vcf" ]; then
    echo "Combined gVCF file for region \$region does not exist. Skipping..."
    exit 1
fi

gatk --java-options "-Xmx400G -XX:+UseParallelGC -XX:ParallelGCThreads=32" GenotypeGVCFs \
   -V \$combined_vcf \
   -O \$output_vcf \
   -R $ref \
   --max-alternate-alleles 2 \
   --max-genotype-count 6
echo "Genotyping complete for region \$region."
EOF

    # Submit the job
    sbatch "$script_path"
done