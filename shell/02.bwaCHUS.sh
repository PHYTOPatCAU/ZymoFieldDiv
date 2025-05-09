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

for i in "${!input_dirs[@]}"; do
    input_dir="${input_dirs[$i]}"
    output_dir="${output_dirs[$i]}"

    for forward_reads in "$input_dir"/*_1.fastq; do
        sample_name=$(basename "$forward_reads" _1.fastq)
        reverse_reads="$input_dir/${sample_name}_2.fastq"
        output_sam="$output_dir/${sample_name}.sam"
        sorted_bam="$output_dir/${sample_name}.bam"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=bwa_${sample_name}
#SBATCH --output=bwa_${sample_name}.out
#SBATCH --error=bwa_${sample_name}.err
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

bwa mem -t 32 $reference $forward_reads $reverse_reads > $output_sam
samtools sort -@ 8 -o $sorted_bam $output_sam
rm $output_sam
EOF

    done
done
