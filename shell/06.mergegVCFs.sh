#!/bin/bash
#SBATCH --job-name=mergeVCF%j
#SBATCH --error=mergeVCF%j.err
#SBATCH --out=mergeVCF%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate var_call

# Define the root directories for gVCF
GVCF_ROOT_DIRS=(
    "~/Zymoproj/UKsamples/gVCF"
    "~/Zymoproj/CHsamples/gVCF"
    "~/Zymoproj/USsamples/gVCF"
)

# Set the reference genome file
REF="~/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

# Iterate over each gVCF root directory
for GVCF_ROOT_DIR in "${GVCF_ROOT_DIRS[@]}"; do
    # Make the list of samples
    samples=$(find "$GVCF_ROOT_DIR" -type f -name '*.g.vcf' | sed 's/^/--variant /')

    # Extract the sample region (e.g., UK, CH, US) from the directory path
    region=$(basename "$(dirname "$GVCF_ROOT_DIR")")

    # Make joint VCF file
    gatk --java-options "-Xmx200G" CombineGVCFs \
        $(echo $samples) \
        -O "$GVCF_ROOT_DIR/Zt_${region}.g.vcf" \
        -R "$REF"
done