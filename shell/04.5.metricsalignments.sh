#!/bin/bash
#SBATCH --job-name=checkcov
#SBATCH --output=checkcov_%j.out
#SBATCH --error=checkcov_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Define the root directories for BAM
BAM_ROOT_DIRS=(
    "~/Zymoproj/UKsamples/BAM"
    "~/Zymoproj/CHsamples/BAM"
    "~/Zymoproj/USsamples/BAM"
)

# Iterate over each BAM root directory
for BAM_ROOT_DIR in "${BAM_ROOT_DIRS[@]}"; do
    # Define the input and output directories
    bam_dir="$BAM_ROOT_DIR/MD"
    output_dir="${BAM_ROOT_DIR%/BAM}/covstats"
    
    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"
    
    # Loop through each BAM file in the MD directory
    for bam_file in "$bam_dir"/*DuplMark.bam; do
        # Get the filename without the extension
        filename=$(basename "$bam_file" .bam)
        
        # Create the output file path
        output_file="$output_dir/$filename.coverage.txt"
        
        # Generate coverage per chromosome using samtools
        samtools depth -a "$bam_file" | awk -v filename="$filename" '{sum[$1]+=$3; count[$1]++} END {for (chr in sum) print filename, chr, sum[chr]/count[chr]}' > "$output_file"
    done
    
    # Get a single file with all the coverage information
    cat "$output_dir"/*.coverage.txt > "$output_dir/all_per_chr_coverage.txt"
    
    # Get coverage stats of the BAM files per sample
    for file in "$output_dir"/*.coverage.txt; do
        # Extract filename without path and extension
        filename_no_ext="${file##*/}"
        filename_base="${filename_no_ext%%.*}"
        
        # Get the coverage stats
        awk '{sum+=$2; sumsq+=$2*$2} END {print "mean:",sum/NR; print "stdev:",sqrt(sumsq/NR - (sum/NR)**2)}' "$file" > "$output_dir/$filename_base.stats"
    done
    
    # Get a single file with all the coverage stats information
    cat "$output_dir"/*.stats > "$output_dir/all_coverage_stats.txt"
    
    # Check that all the headers in the BAM files contain all the chromosomes in all the DuplMark.bam files inside BAM/
    for bamfile in "$bam_dir"/*.DuplMark.bam; do
        # Get the filename without the extension
        filename=$(basename "$bamfile" .bam)
        
        # Create the output file path
        output_file="$output_dir/$filename.perchrmap.txt"
        
        # Generate mapping rate per chromosome using samtools
        samtools idxstats "$bamfile" | awk -v filename="$filename" '{total[$1]+=$3+$4; mapped[$1]+=$3} END {print "Sample\tChromosome\tMapped Reads\tTotal Reads\tMapping Rate"; for(i in total) printf "%s\t%s\t%d\t%d\t%.2f%%\n", filename, i, mapped[i], total[i], (mapped[i]/total[i])*100;}' > "$output_file"
    done
    
    # Get a single file with all the mapping rate per chromosome information
    cat "$output_dir"/*.perchrmap.txt > "$output_dir/all_per_chr_mapping_rate.txt"
done