#!/bin/sh
#SBATCH --job-name=bwa_alignment2_%j
#SBATCH --time=30:00:00
#SBATCH --error=bwa_alignment2_%j.err
#SBATCH --output=bwa_alignment2_%j.out
#Load modules

module load gcc12-env/12.1.0
module load miniconda3/4.12.0

#activate conda env

conda activate samtools

#Previous steps for bwa index your reference genome:

#bwa index -p IPO323 IPO323.fa

#Set input and output directories
input_dir="~/Zymoproj/UKsamples/raw_reads/"
output_dir="~/Zymoproj/UKsamples/BAM_fixed"
reference="IPO323"
#Loop through all pairs of FASTQ files in the input directory

for forward_reads in "$input_dir"/*_R1.fastq.gz; do
# Extract the sample name from the forward reads file
sample_name=$(basename "$forward_reads" _R1.fastq.gz)
# Define the reverse reads file name based on the sample name
reverse_reads="$input_dir/${sample_name}_R2.fastq.gz"
# Define the output SAM file name based on the sample name
output_sam="$output_dir/${sample_name}.sam"
# Submit SLURM job for each pair of reads
sbatch <<EOF

#!/bin/sh
#SBATCH --job-name=bwa_${sample_name}
#SBATCH --output=bwa_${sample_name}.out
#SBATCH --error=bwa_${sample_name}.err
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

#Load modules
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
#load conda env
conda activate samtools
# restate the directories
input_dir="~/Zymoproj/UKsamples/raw_reads/"
output_dir="~/Zymoproj/UKsamples/BAM_fixed"
reference="IPO323"
output_sam="$output_sam"  # Pass the output_sam variable here so the SLURM knows what it is

#Run BWA alignment

bwa mem -t 32 "$reference" "$forward_reads" "$reverse_reads" > "\$output_sam"
#Sort the SAM file

sorted_bam="\$output_dir/${sample_name}.bam"
samtools sort -@ 8 -o "\$sorted_bam" "$output_sam"

#Remove the unsorted SAM file

#rm "\$output_sam"


EOF
done