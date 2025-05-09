#!/bin/bash
#SBATCH --job-name=bwa_alignment_%j
#SBATCH --output=bwa_alignment_%j.out
#SBATCH --error=bwa_alignment_%j.err

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Define variables
reference="IPO323"  # Reference genome
input_dirs=("~/Zymoproj/CHsamples/fastq" "~/Zymoproj/USsamples/fastq")  # Input directories
output_dirs=("~/Zymoproj/CHsamples/BAM" "~/Zymoproj/USsamples/BAM")    # Output directories

# Loop through both input and output directories
for i in "${!input_dirs[@]}"; do
    input_dir="${input_dirs[$i]}"
    output_dir="${output_dirs[$i]}"

    # Process each pair of FASTQ files in the current input directory
    for forward_reads in "$input_dir"/*_1.fastq; do
        # Extract the sample name from the forward reads file
        sample_name=$(basename "$forward_reads" _1.fastq)
        # Define the reverse reads file name based on the sample name
        reverse_reads="$input_dir/${sample_name}_2.fastq"
        # Define the output SAM file name based on the sample name
        output_sam="$output_dir/${sample_name}.sam"

        # Submit SLURM job for each pair of reads
        sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=bwa_${sample_name}
#SBATCH --output=bwa_${sample_name}.out
#SBATCH --error=bwa_${sample_name}.err
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00

# Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# Activate conda environment
conda activate samtools

# Restate the directories
input_dir="$input_dir"
output_dir="$output_dir"
reference="$reference"
output_sam="$output_sam"

# Run BWA alignment
bwa mem -t 32 "$reference" "$forward_reads" "$reverse_reads" > "\$output_sam"

# Sort the SAM file
sorted_bam="\$output_dir/${sample_name}.bam"
samtools sort -@ 8 -o "\$sorted_bam" "\$output_sam"

# Remove the unsorted SAM file
rm "\$output_sam"
EOF
    done
done
