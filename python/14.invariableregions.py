#this script uses as input the bed files generated in 14.invariableregionsbed.sh
import matplotlib.pyplot as plt
from Bio import SeqIO
from pybedtools import BedTool
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import matplotlib.patches as mpatches
import os

# Load reference genome
ref_genome = list(SeqIO.parse("Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa", "fasta"))

# Set the folders containing the BED files
folders = ["/work_beegfs/suaph296/Zymoproj/UKsamples/VCF", "/work_beegfs/suaph296/Zymoproj/USsamples/VCF", "/work_beegfs/suaph296/Zymoproj/CHsamples/VCF"]

# Create a dictionary to store lengths of each chromosome
chrom_lengths = {record.id: len(record) for record in ref_genome}

# Define the function to get telomere length
def get_telomere_length(chrom_id):
    if chrom_id == 'mt':  # mitochondrial genome
        return 0
    chrom_id = int(chrom_id)  # Convert chrom_id to integer
    if 1 <= chrom_id <= 13:  # core chromosomes
        return 134
    elif 14 <= chrom_id <= 21:  # accessory chromosomes
        return 118
    else:
        raise ValueError(f'Invalid chromosome ID: {chrom_id}')

# Create a dictionary to store the positions of the telomeres
telomeres = {}
for chrom in ref_genome:
    chrom_id = chrom.id
    if chrom_id == 'mt':  # Ignore mitochondrial chromosome
        continue
    telomere_length = get_telomere_length(chrom_id)
    telomeres[chrom_id] = (0, telomere_length), (chrom_lengths[chrom_id] - telomere_length, chrom_lengths[chrom_id])

# Calculate the non-telomeric regions
non_telomeric_regions = {}
for chrom_id, chrom_length in chrom_lengths.items():
    telomere_length = get_telomere_length(chrom_id)
    non_telomeric_regions[chrom_id] = (telomere_length, chrom_length - telomere_length)

# Load the table with centromeric regions
centromeric_regions = BedTool("centromeric_regions.bed")

# Define format function for x-axis
def format_func(value, tick_number):
    return f'{value / 1e6}Mbp'

# Iterate over each folder
for folder in folders:
    # Loop through each BED file in the folder
    for bed_file in os.listdir(folder):
        if bed_file.endswith("_invariable_regions.bed"):
            # Load invariable regions from BED file
            invariable_regions = BedTool(os.path.join(folder, bed_file))

            # Subtract the centromeric and telomeric regions from the invariable regions
            true_invariable_regions = invariable_regions.subtract(centromeric_regions)

            # Create a BedTool object for the telomeric regions
            telomeres_bed = BedTool([(chrom_id, start, end) for chrom_id, regions in telomeres.items() for start, end in regions])

            # Subtract the telomeric regions from the true invariable regions
            true_invariable_regions = true_invariable_regions.subtract(telomeres_bed)
            true_invariable_regions.saveas(os.path.join(folder, 'true_invariable_regions.bed'))

            # Filter the true invariable regions to keep only those that are within the non-telomeric regions
            true_invariable_regions_list = []
            for chrom_id, (start, end) in non_telomeric_regions.items():
                true_invariable_regions_filtered = true_invariable_regions.filter(lambda b: b.chrom == chrom_id and b.start > start and b.end < end)
                true_invariable_regions_list.extend(true_invariable_regions_filtered)

            # Create a BedTool object from the list of true invariable regions
            true_invariable_regions = BedTool(true_invariable_regions_list)

            # Plot the true invariable regions along with the centromeres and telomeres
            fig, ax = plt.subplots(figsize=(10, len(chrom_lengths) * 0.5))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))

            # Iterate over chromosomes
            for i, chrom in enumerate(ref_genome):
                chrom_id = chrom.id
                chrom_length = chrom_lengths[chrom_id]

                # Plot chromosome as a line
                ax.plot([0, chrom_length], [i, i], color='black', linewidth=2)

                # Plot true invariable regions
                for region in true_invariable_regions.filter(lambda x: x.chrom == chrom_id):
                    start, end = int(region.start), int(region.end)
                    ax.plot([start, end], [i, i], color='blue', linewidth=7)

                # Plot centromeres
                centromeric_region = [region for region in centromeric_regions if region.chrom == chrom_id]
                if centromeric_region:
                    start, end = centromeric_region[0].start, centromeric_region[0].end
                    ax.plot([start, end], [i, i], color='purple', linewidth=5)

                # Plot telomeres
                if chrom_id != 'mt' and chrom_id in telomeres:
                    for start, end in telomeres[chrom_id]:
                        ax.plot([start, end], [i, i], color='green', linewidth=5)

                # Label chromosome
                ax.text(chrom_length + 0.05 * chrom_length, i, chrom_id, verticalalignment='center')

            # Customize plot
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_func))
            ax.set_xlabel('Position')
            ax.set_yticks(range(len(chrom_lengths)))
            ax.set_yticklabels([chrom.id for chrom in ref_genome])
            ax.set_title(f'Invariable Regions per Chromosome ({bed_file})')
            ax.set_xlim(0, max(chrom_lengths.values()) * 1.1)
            ax.invert_yaxis()  # To have the first chromosome on top
            ax.spines['left'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            # Create a legend
            blue_patch = mpatches.Patch(color='blue', label='Invariable Regions')
            purple_patch = mpatches.Patch(color='purple', label='Centromeres')
            green_patch = mpatches.Patch(color='green', label='Telomeres')
            ax.legend(handles=[blue_patch, purple_patch, green_patch], loc='center right')

            # Save plot as SVG
            plt.tight_layout()
            plt.savefig(os.path.join(folder, f'true_invariable_regions_{bed_file}.svg'), format='svg')
            plt.close()