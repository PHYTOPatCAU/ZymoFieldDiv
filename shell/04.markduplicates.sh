#!/bin/bash
#SBATCH --job-name=mark_duplicates
#SBATCH --output=mark_duplicates_%j.out
#SBATCH --error=mark_duplicates_%j.err
#SBATCH --nodes=2
#SBATCH --mem=50G
#SBATCH --time=24:00:00
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

# Iterate over each BAM root directory
for BAM_ROOT_DIR in "${BAM_ROOT_DIRS[@]}"; do
    # Define the input and output directories
    BamRGFolder="$BAM_ROOT_DIR/RGBAM"
    MD_FOLDER="$BAM_ROOT_DIR/MD"
    
    # Create the MD folder if it doesn't exist
    mkdir -p "$MD_FOLDER"
    
    # Loop through all .bam files in the input folder
    for bamFile in "$BamRGFolder"/*.bam; do
        # Get the sample name from the file name
        Sample=$(basename "$bamFile" .bam)
        
        # Run picard MarkDuplicates with required parameter values
        gatk MarkDuplicates \
            INPUT="$bamFile" \
            OUTPUT="$MD_FOLDER/$Sample.DuplMark.bam" \
            METRICS_FILE="$MD_FOLDER/metrics_$Sample.txt" \
            CREATE_INDEX=true
    done
done