#Using the output from 16.extractEffLoci.sh, this script calculates nucleotide diversity (Pi), Watterson's Theta, and Tajima's D for each effector locus in each region. The script reads the VCF files for each region, computes the statistics, and saves the results to a CSV file for each region and a combined CSV file for all regions. The script uses the scikit-allel library to read the VCF files and compute the statistics. The script also handles cases where the genotype data is not found or the callset is None for a VCF file. The script prints the nucleotide diversity, Watterson's Theta, and Tajima's D for each effector locus in each region. The script saves the results to CSV files in the same directory as the VCF files for each region and a combined CSV file for all regions. The script also prints a message indicating whether the CSV files were saved successfully or not.
import os
import pandas as pd
import allel
import numpy as np

# Base directory containing the region directories
base_directory = '/work_beegfs/suaph296/Zymoproj/merged/'

# Regions to process
regions = ['UK', 'US', 'CH']

# Initialize an empty DataFrame for all regions
all_regions_df = pd.DataFrame(columns=['file', 'Pi', 'Wattersons_Theta', 'Tajimas_D', 'directory', 'region'])

# Loop over each region
for region in regions:
    # Initialize an empty DataFrame for the current region
    region_df = pd.DataFrame(columns=['file', 'Pi', 'Wattersons_Theta', 'Tajimas_D', 'directory'])

    # Directory containing the VCF files for the current region
    region_directory = os.path.join(base_directory, f'{region}_Effectors')

    # Loop over all the VCF files in the region directory
    for filename in os.listdir(region_directory):
        if filename.endswith('relaxed.mac1.recode.vcf.gz'):
            # Full path to the VCF file
            file = os.path.join(region_directory, filename)
            # Read the VCF file
            callset = allel.read_vcf(file)
            # Check if callset is not None and genotype data is present
            if callset is not None and 'calldata/GT' in callset:
                # Get the genotype data
                gt = allel.GenotypeArray(callset['calldata/GT'])
                # Get the chromosome ID from the filename
                chrom_id = callset['variants/CHROM'][0]
                # Count alleles
                ac = gt.count_alleles()
                # Compute Pi
                pi = allel.sequence_diversity(range(len(ac)), ac)
                # Compute Watterson's Theta
                theta_w = allel.watterson_theta(range(len(ac)), ac)
                # Sum over the rows to get the total count for each allele
                ac_sum = ac.sum(axis=0)
                # Compute site frequency spectrum
                sfs = allel.sfs(ac_sum)
                # Compute Tajima's D
                tajimas_d = allel.tajima_d(ac)
                # Handle cases where Tajima's D is nan
                if np.isnan(tajimas_d):
                    tajimas_d = 'indet'
                print(f'Nucleotide diversity for {chrom_id} in {region}: {pi}')
                print(f'Watterson\'s Theta for {chrom_id} in {region}: {theta_w}')
                print(f"Tajima's D for {chrom_id} in {region}: {tajimas_d}")
                # Add the results to the DataFrame if Tajima's D is not 'indet'
                if tajimas_d != 'indet':
                    new_row = pd.DataFrame({'file': [filename], 'Pi': [pi], 'Wattersons_Theta': [theta_w], 'Tajimas_D': [tajimas_d], 'directory': [region_directory]})
                    region_df = pd.concat([region_df, new_row], ignore_index=True)
            else:
                print(f'Genotype data not found or callset is None for file: {filename}')
                # Add a row with 'indet' values
                new_row = pd.DataFrame({'file': [filename], 'Pi': ['indet'], 'Wattersons_Theta': ['indet'], 'Tajimas_D': ['indet'], 'directory': [region_directory]})
                region_df = pd.concat([region_df, new_row], ignore_index=True)

    # Save the region-specific DataFrame to a CSV file
    region_csv = os.path.join(base_directory, f'DivStats{region}.csv')
    region_df.to_csv(region_csv, index=False)

    # Add a column for the region and concatenate with the all regions DataFrame
    region_df['region'] = region
    all_regions_df = pd.concat([all_regions_df, region_df], ignore_index=True)

# Save the combined DataFrame for all regions to a CSV file
all_regions_csv = os.path.join(base_directory, 'DivStatsAll.csv')
all_regions_df.to_csv(all_regions_csv, index=False)

# Debugging: Confirm the files were saved
for region in regions:
    region_csv = os.path.join(base_directory, f'DivStats{region}.csv')
    if os.path.exists(region_csv):
        print(f"CSV file saved successfully for {region} at {region_csv}")
    else:
        print(f"Failed to save CSV file for {region} at {region_csv}")

if os.path.exists(all_regions_csv):
    print(f"Combined CSV file saved successfully at {all_regions_csv}")
else:
    print(f"Failed to save combined CSV file at {all_regions_csv}")