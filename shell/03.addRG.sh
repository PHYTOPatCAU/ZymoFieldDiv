#!/bin/bash
#SBATCH --job-name=AddOrReplaceRG_%j
#SBATCH --error=AddOrReplaceRG_%j.err
#SBATCH --output=AddOrReplaceRG_%j.out
#SBATCH --time=20:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
# path to the conda env
conda activate var_call

# Define the root directories for BAM and BAMRgFolders
BAM_ROOT_DIRS=(
    "~/Zymoproj/UKsamples/BAM"
    "~/Zymoproj/CHsamples/BAM"
    "~/Zymoproj/USsamples/BAM"
)

BAM_RG_ROOT_DIRS=(
    "~/Zymoproj/UKsamples/BAMRgFolders"
    "~/Zymoproj/CHsamples/BAMRgFolders"
    "~/Zymoproj/USsamples/BAMRgFolders"
)

# Iterate over each BAM root directory
for BAM_ROOT_DIR in "${BAM_ROOT_DIRS[@]}"; do
    # Iterate over each BAM file in the current BAM root directory
    for BAM_FILE in "$BAM_ROOT_DIR"/*.bam; do
        # Extract the sample name from the BAM file name
        SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
        
        # Define the corresponding BAMRgFolder
        BAM_RG_FOLDER="${BAM_ROOT_DIR/BAM/BAMRgFolders}/$SAMPLE_NAME"
        
        # Create the BAMRgFolder if it doesn't exist
        mkdir -p "$BAM_RG_FOLDER"
        
        # Add read groups to the BAM file and output to the BAMRgFolder
        gatk AddOrReplaceReadGroups \
            I="$BAM_FILE" \
            O="$BAM_RG_FOLDER/${SAMPLE_NAME}_rg.bam" \
            RGID="$SAMPLE_NAME" \
            RGLB=lib1 \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM="$SAMPLE_NAME"
    done
done