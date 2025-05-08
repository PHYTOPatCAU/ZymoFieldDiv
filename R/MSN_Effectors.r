library(vcfR)
library(ade4)
library(ape)
library(adegenet)
library(poppr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(igraph)
#set wd
setwd("~/Zymoproj/merged/")
# List of VCF file names (with extension)for all three fields functionally characterized effectors
vcf_files <- c("sorted_ZtIPO323_060700.vcf.gz",  # AvrStb6
               "sorted_ZtIPO323_119610.vcf.gz",  # Zt10
               "sorted_ZtIPO323_122170.vcf.gz",  # Zt2
               "sorted_ZtIPO323_017790.vcf.gz",  # Zt13
               "sorted_ZtIPO323_117500.vcf.gz",  # Zt11
               "sorted_ZtIPO323_106990.vcf.gz",  # Zt5
               "sorted_ZtIPO323_019470.vcf.gz",  # Zt12
               "sorted_ZtIPO323_085670.vcf.gz",  # Avr3D1
               "sorted_ZtIPO323_007790.vcf.gz")  # AvrStb9
#effectors snps not found:  "Zt14.vcf.gz" or too few snps: ZtIPO323_012280.vcf.gz
for (vcf_file in vcf_files) {
  # Load the VCF file
  vcf_data <- read.vcfR(vcf_file)
  
  # Filter out loci with more than two alleles
  biallelic_vcf <- vcf_data[is.biallelic(vcf_data), ]
  
  # Change into a table poppr can use
  ww_gl <- vcfR2genlight(biallelic_vcf)
  
  # Read the population data
  pop <- read.csv2("site_WW_SNP_call_2.csv", header = TRUE)
  
  # Put the population data into your table
  pop(ww_gl) <- pop$pop
  ploidy(ww_gl) <- 1 # We have to tell the computer the data is haploid
  
  # Change table type so it is accepted by the distance calculator
  ww_sc <- as.snpclone(ww_gl)
  
  # Calculate genetic distance
  dist <- poppr::bitwise.dist(ww_sc)
  
  # Extract the effector name from the file name
  effector <- sub("\\.vcf\\.gz$", "", vcf_file)
  
  # Extract sample IDs
  sample_ids <- ww_sc@ind.names
  
  # Create a light blue color palette with transparency
  set1_palette <-adjustcolor(c("red", "green", "blue"), alpha.f = 0.5)
  
  # Calculate the minimum spanning network
  msn <- poppr.msn(ww_sc, dist, mlg.compute = "original")
  
  # Make the minimum spanning network plot
svg(file = paste0("MSN_", effector, "_WW.svg"), width = 20, height = 15)
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # Adjust layout to provide space for legends
par(mar = c(5, 5, 2, 2))  # Adjust margins: bottom, left, top, right
plot_poppr_msn(
  x = ww_sc,
  poppr_msn = msn,
  gscale = TRUE,
  gadj = 2,
  mlg.compute = "original",
  glim = c(0, 0.8),
  gweight = 2,
  wscale = FALSE,
  nodescale = 5,
  nodebase = NULL,
  nodelab = 2,
  inds = "n",  # Label all individuals
  mlg = TRUE,  # Label nodes by multilocus genotype
  quantiles = TRUE,
  cutoff = NULL,
  palette = set1_palette,  
  layfun = function(graph) layout_with_drl(graph, options = list(simmer.attraction = 0.1)),  # reduce overlap
  beforecut = FALSE,
  pop.leg = TRUE,  # Ensure population legend appears and does not overlap
  size.leg = FALSE,  # Remove the size legend
  scale.leg = FALSE  # Remove the scale of the distance
)
dev.off()
}
#For within-field effector distribution is only in the UK field. (Per plant GPS position)
setwd("~/Zymoproj/merged/UK_Effectors/")
# List of VCF files to process for the next part
vcf_files <- c(
  "ZtIPO323_117500.vcf.gz",  # Zt11
  "ZtIPO323_119610.vcf.gz",  # Zt10
  "ZtIPO323_060700.vcf.gz",  # AvrStb6
  "ZtIPO323_085670.vcf.gz",  # Avr3D1
  "ZtIPO323_007790.vcf.gz"   # AvrStb9
)

#effectors snps not found: "ZtIPO323_085670.vcf.gz" 
for (vcf_file in vcf_files) {
  # Load the VCF file
  vcf_data <- read.vcfR(vcf_file)
  
  # Filter out loci with more than two alleles
  biallelic_vcf <- vcf_data[is.biallelic(vcf_data), ]
  
  # Change into a table poppr can use
  uk_gl <- vcfR2genlight(biallelic_vcf)
  
  # Read the population data
  pop <- read.table("site_UK_all_Rose_SNP_call_with_tramline.txt", header = TRUE)
  
  # Put the population data into your table
  pop(uk_gl) <- pop$pop
  ploidy(uk_gl) <- 1 # We have to tell the computer the data is haploid
  
  # Change table type so it is accepted by the distance calculator
  uk_sc <- as.snpclone(uk_gl)
  
  # Calculate genetic distance
  dist <- poppr::bitwise.dist(uk_sc)
  
  # Extract the effector name from the file name
  effector <- sub("\\.vcf\\.gz$", "", vcf_file)
  
  # Extract sample IDs
  sample_ids <- uk_sc@ind.names
  # Create a light color palette with transparency
  set1_palette <- adjustcolor(c("purple", "cyan", "magenta"), alpha = 0.5)
  create_mlg_colors <- function(mlg_vector, palette) {
  mlg_sizes <- table(mlg_vector)
  mlg_colors <- rep("gray", length(mlg_sizes))
  large_mlgs <- names(mlg_sizes[mlg_sizes > 5])
  mlg_colors[large_mlgs] <- palette[1:length(large_mlgs)]
  return(setNames(mlg_colors, names(mlg_sizes)))
}
  # Calculate the minimum spanning network
  msn <- poppr.msn(uk_sc, dist, mlg.compute = "original")
  
  # Make the minimum spanning network plot
  svg(file = paste0("MSN_", effector, "_UK.svg"), width = 20, height = 15)
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1))  # Adjust layout to provide space for legends
  par(mar = c(5, 5, 2, 2))  # Adjust margins: bottom, left, top, right
  plot_poppr_msn(
    x = uk_sc,
    poppr_msn = msn,
    gscale = TRUE,
    gadj = 2,
    mlg.compute = "original",
    glim = c(0, 0.8),
    gweight = 2,
    wscale = FALSE,
    nodescale = 5,
    nodebase = NULL,
    nodelab = 2,
    inds = "ALL",  # Label all individuals
    mlg = TRUE,  # Label nodes by multilocus genotype
    quantiles = TRUE,
    cutoff = NULL,
    palette = set1_palette,  
    layfun = function(graph) layout_with_drl(graph, options = list(simmer.attraction = 0.1)),  # reduce overlap
    beforecut = FALSE,
    pop.leg = TRUE,  # Ensure population legend appears and does not overlap
    size.leg = FALSE,  # Remove the size legend
    scale.leg = FALSE  # Remove the scale of the distance
  )
  dev.off()

}

# List of VCF files to process for the next part
vcf_files_of_interest <- c(
  "ZtIPO323_117500.vcf.gz",  # Zt11
  "ZtIPO323_119610.vcf.gz",  # Zt10
  "ZtIPO323_060700.vcf.gz",  # AvrStb6
  "ZtIPO323_085670.vcf.gz",  # Avr3D1
  "ZtIPO323_007790.vcf.gz"   # AvrStb9
)

# Define the MLG colors for each VCF file
mlg_colors_list <- list(
  "ZtIPO323_117500.vcf.gz" = c("150" = "purple", "32" = "cyan", "36" = "magenta", "default" = "grey"),
  "ZtIPO323_119610.vcf.gz" = c("112" = "purple", "94" = "cyan", "default" = "grey"),
  "ZtIPO323_060700.vcf.gz" = c("1" = "purple", "2" = "cyan", "default" = "grey"),
  "ZtIPO323_085670.vcf.gz" = c("141" = "purple", "133" = "cyan", "152" = "magenta", "157" = "green", "default" = "grey"),
  "ZtIPO323_007790.vcf.gz" = c("148" = "purple", "default" = "grey")
)

for (vcf_file in vcf_files_of_interest) {
  # Load the VCF file
  vcf_data <- read.vcfR(vcf_file)
  
  # Filter out loci with more than two alleles
  biallelic_vcf <- vcf_data[is.biallelic(vcf_data), ]
  
  # Change into a table poppr can use
  uk_gl <- vcfR2genlight(biallelic_vcf)
  
  # Read the population data
  pop <- read.table("site_UK_all_Rose_SNP_call_with_tramline.txt", header = TRUE)
  
  # Put the population data into your table
  pop(uk_gl) <- pop$pop
  ploidy(uk_gl) <- 1 # We have to tell the computer the data is haploid
  
  # Change table type so it is accepted by the distance calculator
  uk_sc <- as.snpclone(uk_gl)
  
  # Calculate genetic distance
  dist <- poppr::bitwise.dist(uk_sc)
  
  # Extract the effector name from the file name
  effector <- sub("\\.vcf\\.gz$", "", vcf_file)
  
  # Extract MLG information
  mlg_info <- data.frame(sample.id = indNames(uk_sc), MLG = as.factor(mlg.vector(uk_sc)))
  
  # Merge MLG information with geographic data
  pop_data <- merge(pop, mlg_info, by = "sample.id")
  
  # Remove rows with NA in tramline
  pop_data <- pop_data[!is.na(pop_data$tramline), ]
  
  # Convert lat and long to numeric
  pop_data$lat <- as.numeric(gsub(",", ".", pop_data$lat))
  pop_data$long <- as.numeric(gsub(",", ".", pop_data$long))
  
  # Remove rows where tramline is '1'
  pop_data <- subset(pop_data, tramline != '1')
  
  # Get the MLG colors for the current VCF file
  mlg_colors <- mlg_colors_list[[vcf_file]]
  
  # Create a function to map MLGs to colors, defaulting to grey if not specified
  color_mapping <- function(mlg) {
    if (mlg %in% names(mlg_colors)) {
      return(mlg_colors[mlg])
    } else {
      return(mlg_colors["default"])
    }
  }
  
  # Apply the color mapping function to the MLG column
  pop_data$color <- sapply(pop_data$MLG, color_mapping)
  
  
  # Plot the geographic distribution using ggplot2
  svg(filename = paste0("Geographic_Distribution_", effector, "_UK.svg"), width = 10, height = 8)
  ggplot(pop_data, aes(x = long, y = lat, color = as.factor(MLG))) +
    geom_point(size = 3) +
    geom_jitter(size = 3, width = 0.0001) + 
    scale_color_manual(values = mlg_colors) +
    labs(x = "Longitude", y = "Latitude", color = "MLG") +
    theme_minimal()
  
  dev.off()
}

# Extract plant and leaf information from sample names
sample_info <- strsplit(pop_data$sample.id, "-")
plant_info <- sapply(sample_info, `[`, 2)
leaf_info <- sapply(sample_info, `[`, 3)

# Extract the letter part of the leaf information
leaf_letter <- substr(leaf_info, 1, 1)

# Step 1: Create a data frame combining sample.id, plant_info, leaf_info, and leaf_letter
combined_df <- data.frame(sample.id = pop_data$sample.id, plant = plant_info, leaf = leaf_info, leaf_letter = leaf_letter)

# Step 2: Merge the combined data frame with pop_data
merged_df <- merge(pop_data, combined_df, by = "sample.id")

# Step 3: Calculate same_plant and same_leaf for each row
merged_df <- merged_df %>%
  group_by(plant) %>%
  mutate(same_plant = n() > 1) %>%
  ungroup() %>%
  group_by(plant, leaf_letter) %>%
  mutate(same_leaf = n() > 1) %>%
  ungroup()

# Step 4: Create a mapping from sample names to categories
category_mapping <- merged_df %>%
  mutate(category = ifelse(same_plant, "same_plant",
                           ifelse(same_leaf, "same_leaf", "same_field"))) %>%
  select(sample.id, category)

# Step 5: Add a new column to pop_data that indicates the category of each sample
MLGdistrib <- pop_data %>%
  left_join(category_mapping, by = "sample.id")

# Print the MLGdistrib data frame
print(MLGdistrib)
# Subset the data
same_plant_data <- MLGdistrib %>% filter(category == "same_plant")
same_leaf_data <- MLGdistrib %>% filter(category == "same_leaf")
##barplot of distribution
# Subset the data to include only rows where same_leaf is TRUE
same_leaf_df <- merged_df %>% filter(same_leaf == TRUE)
barplotmlgperplant <- ggplot(same_leaf_df, aes(x = plant, fill = MLG)) +
  geom_bar() +
  labs(title = "MLGs per Plant for Same Leaf Samples", x = "Plant", y = "Count of MLGs") +
  theme_minimal()
svg(file = paste0("barplotmlgperplant_", effector, "_UK.svg"))
print(barplotmlgperplant)
dev.off()
