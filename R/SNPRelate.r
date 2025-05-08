library(tidyverse)
library(ggrepel)
library(SNPRelate)
library(ggplot2)
# UK
# Define the list of VCF files and corresponding output names
vcf_files <- list(
  relaxed = "~/Zymoproj/merged/Zt_UK.relaxed.mac1.recode.vcf.gz",
  strict = "~/Zymoproj/merged/Zt_UK.max-m-80.biallelic-only.mac1.recode.vcf.gz"
)

# Loop through each VCF file
for (dataset in names(vcf_files)) {
  vcf_file <- vcf_files[[dataset]]
  gds_file <- paste0("Zt_UK_", dataset, ".gds")
  
  # Convert VCF to GDS
  snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only") 
  genofile <- openfn.gds(gds_file)
  
  # Perform PCA
  pca <- snpgdsPCA(genofile, autosome.only = FALSE)
  
  # Read the population information
  info <- read.csv2(file = "site_UK_all_SNP_call.csv", header = TRUE)
  
  # Create a data frame with PCA results
  tab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],
                    PC2 = pca$eigenvect[,2],
                    stringsAsFactors = FALSE)
  tab_full <- tab %>% 
    left_join(info, by = "sample.id")
  
  # Calculate the percentage of variance explained by each PC
  pc.percent <- pca$varprop * 100 
  label1 <- paste("PC", 1, " ", format(pc.percent[1], digits = 2), "%", sep = "")
  label2 <- paste("PC", 2, " ", format(pc.percent[2], digits = 2), "%", sep = "")
  
  # Create the labels column in the tab_full dataframe
  tab_full$labels <- tab_full$sample.id
  
  # Create the PCA plot
  p1 <- ggplot(tab_full, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light() +
    geom_text_repel(aes(label = labels), size = 3.5, max.overlaps = 20, box.padding = 0.5) +
    theme()
  
  # Save the plot to a file
  svg(file = paste0("PCA_UK_", dataset, "_2.svg")) 
  print(p1)
  dev.off()
  
  # Remove outliers
  pc1_95 <- quantile(tab_full$PC1, 0.95)
  pc2_95 <- quantile(tab_full$PC2, 0.95)
   tab_full_NO <- tab_full %>% 
  filter(PC1 <= pc1_95 & PC2 <= pc2_95)
  
  # Create the PCA plot without specific rows
  p2 <- ggplot(tab_full_NO, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light() +
    geom_text_repel(aes(label = labels), size = 3.5, max.overlaps = 20, box.padding = 0.5) +
    theme()
  
  # Save the plot to a file
  svg(file = paste0("PCA_UK_", dataset, "2NO.svg")) 
  print(p2)
  dev.off()
  
  # Create the PCA plot without labels
  p3 <- ggplot(tab_full_NO, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light()
  
  # Save the plot to a file
  svg(file = paste0("PCA_UK_", dataset, "NO2_nolbl.svg")) 
  print(p3)
  dev.off()
}
###############################################################################################################
#CH
# Define the list of VCF files and corresponding output names
vcf_files <- list(
  relaxed = "~/Zymoproj/merged/Zt_CH.relaxed.mac1.recode.vcf.gz",
  strict = "~/Zymoproj/merged/Zt_CH.max-m-80.biallelic-only.mac1.recode.vcf.gz"
)

# Loop through each VCF file
for (dataset in names(vcf_files)) {
  vcf_file <- vcf_files[[dataset]]
  gds_file <- paste0("Zt_CH_", dataset, ".gds")
  
  # Convert VCF to GDS
  snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only") 
  genofile <- openfn.gds(gds_file)
  
  # Perform PCA
  pca <- snpgdsPCA(genofile, autosome.only = FALSE)
  
  # Read the population information
  info <- read.csv2(file = "site_CH_SNP_call2.csv", header = TRUE)
  
  # Create a data frame with PCA results
  tab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],
                    PC2 = pca$eigenvect[,2],
                    stringsAsFactors = FALSE)
  tab_full <- tab %>% 
    left_join(info, by = "sample.id")
  
  # Calculate the percentage of variance explained by each PC
  pc.percent <- pca$varprop * 100 
  label1 <- paste("PC", 1, " ", format(pc.percent[1], digits = 2), "%", sep = "")
  label2 <- paste("PC", 2, " ", format(pc.percent[2], digits = 2), "%", sep = "")
  
  # Create the labels column in the tab_full dataframe
  tab_full$labels <- tab_full$sample.id
  
  # Create the PCA plot
  p1 <- ggplot(tab_full, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3.5) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light() +
    geom_text_repel(aes(label = sample.id), size = 3.5, max.overlaps = 20, box.padding = 0.5, colour = "black") +
    theme()
  
  # Save the plot to a file
  svg(file = paste0("PCA_CH_", dataset, "_2.svg")) 
  print(p1)
  dev.off()
  
  # Remove outliers
  pc1_95 <- quantile(tab_full$PC1, 0.95)
  pc2_95 <- quantile(tab_full$PC2, 0.95)
   tab_full_NO <- tab_full %>% 
  filter(PC1 <= pc1_95 & PC2 <= pc2_95)
  
  # Create the PCA plot without specific rows
  p2 <- ggplot(tab_full_NO, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3.5) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light()
  
  # Save the plot to a file
  svg(file = paste0("PCA_CH_", dataset, "2_NO.svg")) 
  print(p2)
  dev.off()
}
###############################################################################################################
#US
# Define the list of VCF files and corresponding output names
vcf_files <- list(
  relaxed = "~/Zymoproj/merged/Zt_US.relaxed.mac1.recode.vcf.gz",
  strict = "~/Zymoproj/merged/Zt_US.max-m-80.biallelic-only.mac1.recode.vcf.gz"
)

# Loop through each VCF file
for (dataset in names(vcf_files)) {
  vcf_file <- vcf_files[[dataset]]
  gds_file <- paste0("Zt_US_", dataset, ".gds")
  
  # Convert VCF to GDS
  snpgdsVCF2GDS(vcf_file, gds_file, method="biallelic.only") 
  genofile <- openfn.gds(gds_file)
  
  # Perform PCA
  pca <- snpgdsPCA(genofile, autosome.only = FALSE)
  
  # Read the population information
  info <- read.csv2(file = "site_US_SNP_call.csv", header = TRUE)
  
  # Create a data frame with PCA results
  tab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],
                    PC2 = pca$eigenvect[,2],
                    stringsAsFactors = FALSE)
  tab_full <- tab %>% 
    left_join(info, by = "sample.id")
  
  # Calculate the percentage of variance explained by each PC
  pc.percent <- pca$varprop * 100 
  label1 <- paste("PC", 1, " ", format(pc.percent[1], digits = 2), "%", sep = "")
  label2 <- paste("PC", 2, " ", format(pc.percent[2], digits = 2), "%", sep = "")
  
  # Create the labels column in the tab_full dataframe
  tab_full$labels <- tab_full$sample.id
  
  # Create the PCA plot
  p1 <- ggplot(tab_full, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3.5) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light() +
    geom_text_repel(aes(label = sample.id), size = 3.5, max.overlaps = 20, box.padding = 0.5, colour = "black") +
    theme()
  
  # Save the plot to a file
  svg(file = paste0("PCA_US_", dataset, "_2.svg")) 
  print(p1)
  dev.off()
  
  # Remove outliers
  pc1_95 <- quantile(tab_full$PC1, 0.95)
  pc2_95 <- quantile(tab_full$PC2, 0.95)
   tab_full_NO <- tab_full %>% 
  filter(PC1 <= pc1_95 & PC2 <= pc2_95)
  
  # Create the PCA plot without specific rows
  p2 <- ggplot(tab_full_NO, aes(x = PC1, y = PC2, color = pop)) +
    geom_point(size = 3.5) +
    xlab(label1) +
    ylab(label2) +
    labs(color = 'Field') +
    theme_light()
  
  # Save the plot to a file
  svg(file = paste0("PCA_US_", dataset, "2_NO.svg")) 
  print(p2)
  dev.off()
}
###############################################################################################################
#WW
# Define the VCF file
vcf_file <- "Zt_WW_relaxed_core_v2.vcf.gz"

# Convert VCF to GDS
snpgdsVCF2GDS(vcf_file, "Zt_WW_v2.gds", method="biallelic.only") 
genofile <- openfn.gds("Zt_WW_v2.gds")

# Perform PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE)

# Read the population information
info <- read.csv(file = "site_WW_SNP_call_3_justpop.csv", header = TRUE)

# Create a data frame with PCA results
tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],
                  PC2 = pca$eigenvect[,2],
                  stringsAsFactors = FALSE)
tab_full <- tab %>% 
  left_join(info, by = "sample.id")

# Calculate the percentage of variance explained by each PC
pc.percent <- pca$varprop * 100 
label1 <- paste("PC", 1, " ", format(pc.percent[1], digits = 2), "%", sep = "")
label2 <- paste("PC", 2, " ", format(pc.percent[2], digits = 2), "%", sep = "")

# Create the PCA plot
p1 <- tab_full %>% 
  ggplot(aes(x = PC1, y = PC2, color = pop)) + 
  xlab(label1) +
  ylab(label2) + 
  geom_point(size = 3.5) + 
  labs(color = 'Field') +
  theme_light() + 
  geom_text_repel(aes(label = sample.id), size = 3.5, max.overlaps = 20, box.padding = 0.5) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_text(face = "bold"))

# Save the plot to a file
svg(file = "PCA_WW_all.svg") 
print(p1)
dev.off()

# Remove outliers
  pc1_95 <- quantile(tab_full$PC1, 0.95)
  pc2_95 <- quantile(tab_full$PC2, 0.95)
  tab_full_NO <- tab_full %>% 
  filter(PC1 <= pc1_95 & PC2 <= pc2_95)

# Create the PCA plot without specific rows
p2 <- tab_full_NO %>% 
  ggplot(aes(x = PC1, y = PC2, color = pop)) + 
  xlab(label1) +
  ylab(label2) + 
  geom_point(size = 3.5) + 
  labs(color = 'Field') +
  theme_light()

# Save the plot to a file
svg(file = "PCA_WW_all_NO.svg") 
print(p2)
dev.off()

# Create the PCA plot without labels
p3 <- tab_full %>% 
  ggplot(aes(x = PC1, y = PC2, color = pop)) + 
  xlab(label1) +
  ylab(label2) + 
  geom_point(size = 3.5) + 
  labs(color = 'Field') +
  theme_light()

# Save the plot to a file
svg(file = "PCA_WW_all_nolbl.svg") 
print(p3)
dev.off()
