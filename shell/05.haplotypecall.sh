#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate the conda environment
conda activate var_call

# Define the root directories for BAM
BAM_ROOT_DIRS=(
    "/work_beegfs/suaph296/Zymoproj/UKsamples/BAM"
    "/work_beegfs/suaph296/Zymoproj/CHsamples/BAM"
    "/work_beegfs/suaph296/Zymoproj/USsamples/BAM"
)

# Set the reference genome file
REF=Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa

# Iterate over each BAM root directory
for BAM_ROOT_DIR in "${BAM_ROOT_DIRS[@]}"; do
    # Define the input and output directories
    BAMFolder="$BAM_ROOT_DIR/MD"
    gVCFFolder="${BAM_ROOT_DIR%/BAM}/gVCF"
    
    # Create the gVCF folder if it doesn't exist
    mkdir -p "$gVCFFolder"
    
    # Iterate through each DuplMark.bam file in the MD folder
    for bamfile in "$BAMFolder"/*.DuplMark.bam; do
        # Extract the sample name from the file name
        Sample=$(basename "$bamfile" .DuplMark.bam)

        # Start a separate SLURM job for each DuplMark.bam file
        sbatch <<EOF
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=HaplotypeCaller_${Sample}

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate the conda environment
conda activate var_call

# GATK HaplotypeCaller command
gatk --java-options "-Xmx200G" HaplotypeCaller -R $REF -ploidy 1 --emit-ref-confidence GVCF -I "$bamfile" -O "$gVCFFolder/${Sample}.g.vcf"
EOF

    done
done