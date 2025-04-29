#!/bin/bash

# Step 1: Merge the VCF files into a single file
bcftools merge \
    Zymoproj/CHsamples/VCF/Zt_CH2.SNP.qualityfilter.excl.rn1.vcf.gz \
    Zymoproj/USsamples/VCF/Zt_US.SNP.qualityfilter.excl.rn1.vcf.gz \
    Zymoproj/UKsamples/VCF/Zt_UK.SNP.qualityfilter.excl.rn1.vcf.gz \
    -Oz -o Zymoproj/merged/WW_corrected_relaxed.vcf.gz

# Step 2: Index the merged VCF file
bcftools index Zymoproj/merged/WW_corrected_relaxed.vcf.gz

# Step 3: Run the annotation 
java -mx200G -jar snpEff/snpEff.jar -v Zymoseptoria_tritici Zymoproj/merged/Zt_WW_relaxed_mac1.recode.vcf.gz > snpEff/SNP_eff_Zt_WW_relaxed.vcf

# Filter the VCF to keep only SNPs with MODIFIER, HIGH, or MODERATE impact
java -Xmx200G -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" snpEff/SNP_eff_Zt_WW_relaxed.vcf > snpEff/SNP_eff_Zt_WW_relaxed_impact.vcf

# Extract the regions of interest from the filtered VCF
#cat snpEff/SNP_eff_Zt_WW_strict_impact.vcf | perl snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx50G -jar snpEff/SnpSift.jar extractFields - CHROM POS REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" > snpEff/SNP_filtered_Zt_WW_strict.eff.txt
# Fiter to keep only the variants with AF < 0.05 for the recovered SNPs
#java -Xmx200G -jar snpEff/SnpSift.jar filter "AF < 0.05" snpEff/SNP_eff_Zt_WW_relaxed_impact.vcf > snpEff/SNP_eff_Zt_WW_relaxed_rare_events.vcf

# Define the VCF file
vcf_file="snpEff/SNP_eff_Zt_WW_relaxed_impact.vcf"

# Define the variant types to count
variant_types=(
    "missense_variant"
    "inframe_insertion"
    "inframe_deletion"
    "stop_gained"
    "frameshift_variant"
    "start_lost"
    "splice_acceptor_variant"
    "splice_donor_variant"
)

# Loop through each variant type and count occurrences
for variant in "${variant_types[@]}"; do
    count=$(grep -o "$variant" "$vcf_file" | wc -l)
    echo "$variant: $count"
done
#for effectors
# Step 3: Run the annotation 
java -mx200G -jar snpEff/snpEff.jar -v Zymoseptoria_tritici /work_beegfs/suaph296/Zymoproj/merged/GOI_WW_strict/Effectors_all_merged.vcf.gz > snpEff/SNP_eff_Zt_WW_strict.vcf

#extract the regions of interest
#cat snpEff/SNP_eff_Zt_WW_strict.vcf | perl snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx50G -jar snpEff/SnpSift.jar extractFields - CHROM POS REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" > snpEff/SNP_filtered_Zt_WW_strict.eff.txt

# Filter the VCF to keep only SNPs with MODIFIER, HIGH, or MODERATE impact
java -Xmx200G -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE')" snpEff/SNP_eff_Zt_WW_strict.vcf > snpEff/SNP_eff_Zt_WW_strict_impact.vcf

# Extract the regions of interest from the filtered VCF
#cat snpEff/SNP_eff_Zt_WW_strict_impact.vcf | perl snpEff/scripts/vcfEffOnePerLine.pl | java -Xmx50G -jar snpEff/SnpSift.jar extractFields - CHROM POS REF ALT AF "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" > snpEff/SNP_filtered_Zt_WW_strict.eff.txt
# Fiter to keep only the variants with AF < 0.05 for the recovered SNPs
java -Xmx200G -jar snpEff/SnpSift.jar filter "AF < 0.05" snpEff/SNP_eff_Zt_WW_strict_impact.vcf > snpEff/SNP_eff_Zt_WW_strict_rare_events.vcf

# Define the VCF file
vcf_file="snpEff/SNP_eff_Zt_WW_strict_rare_events.vcf"

# Define the variant types to count
variant_types=(
    "missense_variant"
    "inframe_insertion"
    "inframe_deletion"
    "stop_gained"
    "frameshift_variant"
    "start_lost"
    "splice_acceptor_variant"
    "splice_donor_variant"
)

# Loop through each variant type and count occurrences
for variant in "${variant_types[@]}"; do
    count=$(grep -o "$variant" "$vcf_file" | wc -l)
    echo "$variant: $count"
done
