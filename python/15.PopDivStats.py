import os
import pandas as pd
import allel
import numpy as np
import glob

# Directory containing the sorted VCF files
directory = '/work_beegfs/suaph296/Zymoproj/merged'

# Initialize an empty DataFrame
df = pd.DataFrame(columns=['file', 'chromosome', 'window_start', 'window_end', 'Pi', 'Wattersons_Theta', 'Tajimas_D', 'filter'])

# Use glob to find all strict and relaxed VCF files
strict_files = glob.glob(os.path.join(directory, '*.max-m-80.biallelic-only.mac1.recode.vcf.gz'))
relaxed_files = glob.glob(os.path.join(directory, '*.SNP.qualityfilter.excl.rn1.vcf.gz'))

# Combine strict and relaxed files into one list
vcf_files = [(file, 'strict') for file in strict_files] + [(file, 'relaxed') for file in relaxed_files]

# Loop over all the VCF files
for file, filter_type in vcf_files:
    # Extract the filename from the full path
    filename = os.path.basename(file)
    print(f"Processing file: {filename} ({filter_type})")

    # Read the VCF file
    try:
        callset = allel.read_vcf(file, fields=['calldata/GT', 'variants/CHROM', 'variants/POS'])
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        continue

    # Get the genotype data
    if 'calldata/GT' not in callset:
        print(f"No genotype data found in file {filename}. Skipping.")
        continue

    gt = allel.GenotypeArray(callset['calldata/GT'])
    chroms = callset['variants/CHROM']
    pos = callset['variants/POS']

    # Process each chromosome separately
    unique_chroms = np.unique(chroms)
    for chrom in unique_chroms:
        # Filter data for the current chromosome
        chrom_mask = (chroms == chrom)
        chrom_pos = pos[chrom_mask]
        chrom_gt = gt.compress(chrom_mask, axis=0)

        # Skip if no data for this chromosome
        if chrom_gt.shape[0] == 0:
            print(f"No data for chromosome {chrom} in file {filename}. Skipping.")
            continue

        # Count alleles
        chrom_ac = chrom_gt.count_alleles()

        # Ensure positions and allele counts are sorted
        sorted_indices = np.argsort(chrom_pos)
        chrom_pos = chrom_pos[sorted_indices]
        chrom_ac = chrom_ac[sorted_indices]

        # Skip if allele count array is empty or invalid
        if chrom_ac.shape[0] == 0 or chrom_ac.ndim != 2:
            print(f"Invalid allele count data for chromosome {chrom} in file {filename}. Skipping.")
            continue

        # Define window size and step
        window_size = 100000  # 100 kb
        step_size = 100000    # 100 kb

        # Process data in windows
        for window_start in range(0, chrom_pos[-1], step_size):
            window_end = window_start + window_size
            window_mask = (chrom_pos >= window_start) & (chrom_pos < window_end)
            window_ac = chrom_ac[window_mask]

            # Skip if no data in the window
            if window_ac.shape[0] == 0:
                continue

            # Compute Pi, Watterson's Theta, and Tajima's D
            try:
                pi = allel.sequence_diversity(chrom_pos[window_mask], window_ac)
                wattersons_theta = allel.watterson_theta(chrom_pos[window_mask], window_ac)
                tajimas_d = allel.tajima_d(window_ac)
            except Exception as e:
                print(f"Error processing file {filename}, chromosome {chrom}, window {window_start}-{window_end}: {e}")
                continue

            # Add the results to the DataFrame
            new_row = pd.DataFrame({
                'file': [filename],
                'chromosome': [chrom],
                'window_start': [window_start],
                'window_end': [window_end],
                'Pi': [pi],
                'Wattersons_Theta': [wattersons_theta],
                'Tajimas_D': [tajimas_d],
                'filter': [filter_type]
            })
            df = pd.concat([df, new_row], ignore_index=True)

# Export to a CSV file
output_csv = f'{directory}/divstatsBigFilt_windows.csv'
df.to_csv(output_csv, index=False)
print(f"CSV file saved successfully at {output_csv}")
