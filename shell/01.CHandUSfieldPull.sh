#!/bin/bash
#SBATCH --job-name=runtoname_%j
#SBATCH --error=runtoname_%j.err
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=08:00:00
# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
module load curl/7.83.0
# Activate conda environment
conda activate data_download
conda install -c bioconda sra-tools
conda install -c bioconda entrez-direct
 
#search for the list of accessions inside the bioproject as in Feurtey et al. 2023
esearch -db sra -query "PRJNA596434" | efetch -format runinfo | cut -f 1 -d ',' > run_accessions.txt
##esearch -db bioproject -query "PRJNA596434" | efetch -format docsum | xtract -pattern DocumentSummary -element Run@acc > run_accessions.txt
#retrieve the fastq files
while read sra_run; do
  fastq-dump --split-files $sra_run
done < run_accessions.txt
#rename the files
# Activate conda environment if not already activated
conda activate samtools
# Define the input directory where the FASTQ files and response.txt are located
input_dir="/work_beegfs/suaph296/Zymoproj/DEsamples/fastq/"
response_file="${input_dir}/response_1718702682986.txt.txt"

# Skip the header line of the response.txt file and read the rest
tail -n +2 "$response_file" | while IFS=$'\t' read -r run name barcode sample_status action analysis_id path_r1 path_r2 file_status_r1 file_status_r2 md5_checksum_r1 md5_checksum_r2; do
    # Construct the current and new file names for R1 and R2
    current_r1="${input_dir}${path_r1##*/}"
    new_r1="${input_dir}${name}_S1_R1_001.fastq.gz"
    current_r2="${input_dir}${path_r2##*/}"
    new_r2="${input_dir}${name}_S1_R2_001.fastq.gz"

    # Check if the current files exist before renaming
    if [[ -f "$current_r1" && -f "$current_r2" ]]; then
        echo "Renaming $current_r1 to $new_r1"
        echo "Renaming $current_r2 to $new_r2"
        mv "$current_r1" "$new_r1"
        mv "$current_r2" "$new_r2"
    else
        echo "Files for $name not found, skipping..."
    fi
done