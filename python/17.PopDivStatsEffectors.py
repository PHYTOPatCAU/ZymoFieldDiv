# This script uses the input VCFfiles from 15.splitVCFtoChr.sh
import os
import pandas as pd
import allel

# Initialize an empty DataFrame
df = pd.DataFrame(columns=['file', 'region', 'chromosome', 'Pi', 'Tajimas_D', 'Wattersons_Theta'])

# Directories containing the VCF files
directories = {
    'UK': '~/Zymoproj/merged/UKsamples',
    'US': '~/Zymoproj/merged/USsamples',
    'CH': '~/Zymoproj/merged/CHsamples'
}

# Loop over all the directories
for region, directory in directories.items():
    # Loop over all the VCF files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.vcf.gz'):
            # Full path to the VCF file
            file = os.path.join(directory, filename)
            # Read the VCF file
            callset = allel.read_vcf(file)
            # Get the genotype data
            gt = allel.GenotypeArray(callset['calldata/GT'])
            # Get the chromosome ID from the filename
            chrom_id = callset['variants/CHROM'][0]  
            # Count alleles
            ac = gt.count_alleles()
            # Compute Pi
            pi = allel.sequence_diversity(range(len(ac)), ac)
            # Compute Watterson's theta
            theta_w = allel.watterson_theta(range(len(ac)), ac)
            # Sum over the rows to get the total count for each allele
            ac_sum = ac.sum(axis=0)
            # Compute site frequency spectrum
            sfs = allel.sfs(ac_sum)
            # Compute Tajima's D
            tajimas_d = allel.tajima_d(ac)
            print(f'Nucleotide diversity for {chrom_id} in {region}: {pi}')
            print(f'Site Frequency Spectrum for {chrom_id} in {region}: {sfs}')
            print(f"Tajima's D for {chrom_id} in {region}: {tajimas_d}")
            print(f"Watterson's Theta for {chrom_id} in {region}: {theta_w}")
            # Add the results to the DataFrame
            new_row = pd.DataFrame({
                'file': [filename],
                'region': [region],
                'chromosome': [chrom_id],
                'Pi': [pi],
                'Tajimas_D': [tajimas_d],
                'Wattersons_Theta': [theta_w]
            })
            df = pd.concat([df, new_row], ignore_index=True)

# Print the DataFrame
print(df)

# Export to a csv file
df.to_csv('~/Zymoproj/merged/divstats_all_eff_regions.csv', index=False)