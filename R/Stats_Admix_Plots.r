#This is a collection of miniscripts used to plot the data from various scripts: 04.5,13,15,17
library(ggplot2)
library(ggpattern)
library(tidyverse)
library(reshape2)
library(scales)
library(data.table)
#set wd
setwd("/work_beegfs/suaph296/Zymoproj/merged/")

#Quality check for depth per chromosome
# Define the directories
dirs <- list(
  UK = "/work_beegfs/suaph296/Zymoproj/UKsamples",
  US = "/work_beegfs/suaph296/Zymoproj/USsamples",
  CH = "/work_beegfs/suaph296/Zymoproj/CHsamples"
)

# Define the output directory
output_dir <- "/work_beegfs/suaph296/Zymoproj/merged"

# Define the accessory contigs
accessory_contigs <- c("chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21")

# Define the colors and the values at which the colors should change for raw heatmap
colors_raw <- c("black", "red", "yellow", "white")
values_raw <- c(1, 10, 2000, 30000)

# Define the colors and the values at which the colors should change for accessory contigs heatmap
colors_acc <- c("white", "yellow", "red", "black")
values_acc <- c(2, 1.5, 0.625 / 1.5, 0.25)

# Iterate through each directory
for (region in names(dirs)) {
  dir <- dirs[[region]]
  
  # Load the data
  depthchr <- read.csv(file.path(dir, "covstats", "all_per_chr_coverage.txt"), sep = "\t")
  depthchr <- depthchr %>% group_by(sample)
  depthperhrnm <- depthchr %>% mutate(depth = as.numeric(depth))
  
  # Rename and filter chromosomes
  chr <- depthchr %>% 
    mutate(chr = if_else(chr %in% as.character(1:21), paste0("chr", chr), chr)) %>% 
    filter(chr %in% c(paste0("chr", 1:21), "mt"))
  
  # Filter the data for chromosomes from chr1 to chr21
  only_chr <- chr %>% filter(chr %in% c(paste0("chr", 1:21), "mt"))
  
  # Calculate mean depth for each sample
  mean_depth_per_sample <- only_chr %>%
    group_by(sample) %>%
    summarize(mean_depth = mean(depth, na.rm = TRUE))
  
  # Calculate the overall mean coverage across all samples
  overall_mean_depth <- mean(mean_depth_per_sample$mean_depth, na.rm = TRUE)
  print(paste("Overall mean depth for", region, ":", overall_mean_depth))
  
  # Calculate the total coverage for each sample
  total_coverage <- chr %>% group_by(sample) %>% summarise(total = sum(depth))
  
  # Order the samples by total coverage
  chr$sample <- factor(chr$sample, levels = total_coverage$sample[order(total_coverage$total)])
  
  # Order the chromosomes
  chr$chr <- factor(chr$chr, levels = c(paste0("chr", 1:21), "mt"))
  
  # Recreate the heatmap with the custom color gradient for raw heatmap
  heatmapraw <- ggplot(chr, aes(x = chr, y = sample, fill = depth)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors_raw, values = scales::rescale(values_raw)) +
    theme_minimal() +
    labs(x = "Chromosome", y = "Sample", fill = "Coverage")
  
  # Save the heatmap to a file
  svg(file.path(output_dir, paste0("heatmap_", region, ".svg")), width = 20, height = 30)
  print(heatmapraw)
  dev.off()
  
  # Subset only the accessory contigs
  acc <- only_chr %>% filter(chr %in% accessory_contigs)
  
  # Divide each chr depth by the mean core coverage already calculated in the mean_per_sample
  normacc <- acc %>% 
    left_join(mean_depth_per_sample, by = "sample") %>% 
    mutate(depth = depth / mean_depth)
  
  # Filter out specific samples
  filtered_normacc <- normacc %>% filter(!sample %in% c("FU_008", "FU_110"))
  
  # Recreate the heatmap with the custom color gradient for accessory contigs
  heatmapacc <- ggplot(filtered_normacc, aes(x = chr, y = sample, fill = depth)) +
    geom_tile() +
    scale_fill_gradientn(colors = colors_acc, values = scales::rescale(values_acc), breaks = c(1, 2, 3, 4, 5, 8)) +
    theme_minimal() +
    labs(x = "Chromosome", y = "Sample", fill = "Coverage")
  
  # Save the heatmap to a file
  svg(file.path(output_dir, paste0("heatmap_acc_", region, ".svg")), width = 20, height = 30)
  print(heatmapacc)
  dev.off()
}
#Genotyping rate per field file
# Define the regions and types
regions <- c("UK", "US", "CH")
types <- c("relaxed", "strict")

# Iterate through each region and type
for (region in regions) {
  for (type in types) {
    # Define the file names
    imiss_file <- paste0("Zt_", region, ".", type, ".imiss")
    frq_file <- paste0("Zt_", region, ".", type, ".frq")
    
    # Load the PLINK .imiss file into a data frame
    data <- read.delim(imiss_file, sep="", header=TRUE)
    
    # Calculate the genotyping rate for each sample & add the column
    data <- data %>%
      mutate(GenotypingRate = N_GENO / (N_GENO + N_MISS) * 100)
    
    # Create the bar plot for genotyping rate
    svg(file=paste0(region, "_GT_", type, ".svg"))
    ggplot(data, aes(x=reorder(1:nrow(data), -GenotypingRate), y=GenotypingRate)) +
      geom_bar(stat="identity") +
      labs(x="", y="Genotyping Rate (%)") +
      theme_bw() +
      geom_hline(yintercept = mean(data$GenotypingRate), color = "skyblue", linetype = "dashed")
    dev.off()
    
    ## MAF calc and plot
    # Load maf data
    maf_data <- fread(frq_file)
    
    # Plot raw MAF
    mafplot <- ggplot(maf_data, aes(x = MAF)) +
      geom_histogram(bins = 20) +
      labs(x = "MAF", y = "SNP count") +
      ylim(0, 1e6) +
      theme_bw()
    svg(file=paste0(region, "_MAF_", type, "_20.svg"))
    print(mafplot)
    dev.off()
    
    # Make the histogram normalized by sample size
    maf_data <- maf_data %>%
      mutate(NormCount = NCHROBS / sum(NCHROBS))
    mafplot <- ggplot(maf_data, aes(x = MAF)) +
      geom_histogram(bins = 20) +
      labs(title = "Minor Allele Frequency (MAF) Distribution",
           x = "MAF", y = "SNP count") +
      ylim(0, 1e6) +
      theme_bw() +
      scale_y_continuous(labels = function(x) paste0(x / 1e6, "M SNPs"), limits = c(0, 1e6))
    
    # Save the normalized MAF plot
    svg(file=paste0(region, "_MAF_norm_", type, "_20.svg"))
    print(mafplot)
    dev.off()
  }
}
#SNP per chr counts
# Read the data from the file
data <- read.table("~/Desktop/Work/Uni-Kiel/snp_counts_all_corrected_mac1.txt", header = TRUE)
contiglength <- read.table("~/Desktop/Work/Uni-Kiel/contiglengthrefZt.txt", h=T)
# Merge the data frames
merged_data <- merge(data, contiglength, by = "chr")
# Create the new column
norm <- merged_data %>%
  mutate(length = length / 1e3) %>%
  rename(length_in_Kb = length)
norm <- norm %>%
  mutate(SNPs_per_Kb = snp_count / length_in_Kb)
# Change the chromosome names
norm$chr <- paste0("chr", norm$chr)

# Convert the chromosome column to a factor and specify the levels
norm$chr <- factor(norm$chr, levels = paste0("chr", 1:21))

# Create the bar plot
# Create a new variable for alpha levels
norm$alpha <- ifelse(norm$filters %in% c("Strict", "Relaxed"), norm$filters, NA)

bar_plot <- ggplot() +
  geom_bar(data = subset(norm, filters == "Strict"), aes(x = chr, y = SNPs_per_Kb, fill = country, alpha = alpha), stat = "identity", position = "dodge") +
  geom_bar(data = subset(norm, filters == "Relaxed"), aes(x = chr, y = SNPs_per_Kb, fill = country, alpha = alpha), stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") +
  scale_alpha_manual(name = "Filters", values = c("Strict" = 0.5, "Relaxed" = 0.75), guide = "legend") +
  theme_minimal() +
  labs(title="SNP density corrected by Chr length",x = "Chromosome", y = "SNPs per Kb")
bar_plot
# Print the plot
print(bar_plot)
svg("SNPperKb_all_filters_new.svg",  width = 12 ,  height =6.5)
bar_plot
dev.off()

# Load data for Pi and Tajima_D at Chr level
divstatsall <- read.csv("/work_beegfs/suaph296/Zymoproj/merged/DivStatsAll.csv")

# Rename columns to chr, Pi, Tajima_D and country
divstatsall <- divstatsall %>%
  rename(chr = file, Pi = Pi, Tajima_D = Tajimas_D, country = directory)

# Remove the '.vcf.gz' from the chr column
divstatsall$chr <- gsub(".vcf.gz", "", divstatsall$chr)

# Add chr to the beginning of the chr column except for mt
divstatsall$chr <- ifelse(divstatsall$chr == "mt", "mt", paste0("chr", divstatsall$chr))

# Remove the 'samples' from the end of the country column
divstatsall$country <- gsub("samples", "", divstatsall$country)

# Define the chromosome levels
chr_levels <- c(paste0("chr", 1:21), "mt")

# Convert chr to a factor with levels in the correct order
divstatsall$chr <- factor(divstatsall$chr, levels = chr_levels)

# Barplot for the Pi values with fill=country
bar_plot <- ggplot(divstatsall, aes(x = chr, y = Pi, fill = country)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Chromosome", y = "Pi") +
  scale_fill_brewer(palette = "Set1")
bar_plot

# Save the corrected fig
svg("Pi_all_corrected_mac1.svg", width = 9, height = 6.5)
bar_plot
dev.off()

# Barplot for the Tajima_D values with fill=country
bar_plot <- ggplot(divstatsall, aes(x = chr, y = Tajima_D, fill = country)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Chromosome", y = "Tajima_D") +
  scale_fill_brewer(palette = "Set1")
bar_plot

# Save the corrected fig
svg("Tajima_D_all_corrected_mac1.svg", width = 9, height = 6.5)
bar_plot + ylim(-2.5, 1)
dev.off()

#load data dor Pi and Tajima_D per Effector
# Load data
divstatsmac1 <- read.csv("/work_beegfs/suaph296/Zymoproj/merged/divstats_all_eff_regions.csv")

# Rename columns
divstatsmac1 <- divstatsmac1 %>%
  rename(Effector = Effector, Pi = Pi, Wattersons_Theta = Wattersons_Theta, Tajima_D = Tajimas_D, region = country)

# Remove 'vcf.gz' from Effector names
divstatsmac1 <- divstatsmac1 %>%
  mutate(Effector = str_remove(Effector, ".vcf.gz$"),
         country = str_remove(country, "^Zt_"))

# Remove rows where Pi contains 'indet'
divstatsmac1 <- divstatsmac1 %>%
  filter(!str_detect(Pi, "indet"))

# Convert Pi to numeric
divstatsmac1 <- divstatsmac1 %>%
  mutate(Pi = as.numeric(Pi))

# Reorder Effector based on Pi values
divstatsmac1 <- divstatsmac1 %>%
  arrange(Pi) %>%
  mutate(Effector = factor(Effector, levels = unique(Effector)))

# Histogram for the Pi values with fill=country
pi_histogram <- ggplot(divstatsmac1, aes(x = Effector, y = Pi, fill = country)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "Effector", y = "Pi") 
pi_histogram

# Save the Pi histogram
svg("Pi_histogram_eff2.svg", width = 9, height = 6.5)
pi_histogram
dev.off()
# Remove rows where Pi contains 'indet'
divstatsmac1 <- divstatsmac1 %>%
  filter(!str_detect(Tajima_D, "indet"))
# Convert Tajima_D to numeric

divstatsmac1 <- divstatsmac1 %>%
  mutate(Tajima_D = as.numeric(Tajima_D))

# Reorder Effector based on Tajima_D values
divstatsmac1 <- divstatsmac1 %>%
  arrange(Tajima_D) %>%
  mutate(Effector = factor(Effector, levels = unique(Effector)))

# Histogram for the Tajima_D values with fill=country
tajima_histogram <- ggplot(divstatsmac1, aes(x = Effector, y = Tajima_D, fill = country)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  ylim(-9, 9) +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) +
  labs(x = "Effector", y = "Tajima_D")
tajima_histogram

# Save the Tajima_D histogram
svg("Tajima_D_histogram2.svg", width = 9, height = 6.5)
print(tajima_histogram)
dev.off()

#Admixture plots

CV_error <- read.delim("cv_errors.txt")
colnames(CV_error) <- c("K", "cv_error")
CV_error$K <- gsub("\\(K=", "", CV_error$K)
CV_error$K <- gsub("\\):", "", CV_error$K)
CV_error$K <- factor(CV_error$K, levels = as.character(1:10))
#save as svg
svg(filename = "CV_errors.svg", width = 8, height = 8)
ggplot(CV_error, aes(x = K, y = cv_error)) +
  geom_boxplot() +
  labs(x = "K", y = "CV Error") +
  theme_minimal()
dev.off()
# Add  labels using "site_WW_SNP_call_2.csv"
labels <- read.csv2("site_WW_SNP_call_justpop2.csv",h=T)

# Add a column with population indices to order the barplots
labels$n <- factor(labels$pop, levels = unique(labels$pop))
levels(labels$n) <- c(1:length(levels(labels$n)))
labels$n <- as.integer(as.character(labels$n))

# read in the different admixture output files
minK = 2
maxK = 2
prefix = "maf05pruned_data"
tbl <- lapply(minK:maxK, function(x) read.table(paste0(prefix, ".", x, ".Q")))

# Prepare spaces to separate the populations/species
rep <- as.vector(table(labels$n))
spaces <- 0
for(i in 1:length(rep)) {
  spaces <- c(spaces, rep(0, rep[i] - 1), 0.5)
}
spaces <- spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for K=2
svg(filename = "admixture_plot_K2_corrected.svg", width = 10, height = 7)
par(mar = c(8, 4, 4, 2) + 0.1)
# Plot K=2
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), 
              col = rainbow(n = 2), 
              xaxt = "n", 
              border = NA, 
              ylab = "K=2", 
              yaxt = "n", 
              space = spaces)

# Add the individual names below the barplot with smaller font size
axis(1, at = bp, labels = labels$sample.id[order(labels$n)], las = 2, tick = FALSE, cex = 0.1)

# Add the population names below the individual names
pop.names <- unique(labels$pop[order(labels$n)])
pop.positions <- sapply(pop.names, function(x) mean(bp[labels$pop[order(labels$n)] == x]))
axis(1, at = pop.positions, labels = pop.names, las = 1, tick = FALSE, cex = 0.6, line = +3)
par(mar = c(5, 4, 4, 2) + 0.1)
dev.off()
#BETTER VERSION
# Plot K=2
# Increase the height of the plot (too crowded with ind names)
svg(filename = "admixture_plot_K2_corrected.svg", width = 15, height = 7)

# Plot K=2
bp <- barplot(t(as.matrix(tbl[[1]][order(labels$n),])), 
              col = rainbow(n = 2), 
              xaxt = "n", 
              border = NA, 
              ylab = "K=2", 
              yaxt = "n", 
              space = spaces)

# Add the population names below the barplot
pop.names <- unique(labels$pop[order(labels$n)])
pop.positions <- sapply(pop.names, function(x) mean(bp[labels$pop[order(labels$n)] == x]))
axis(1, at = pop.positions, labels = pop.names, las = 2, tick = FALSE, cex = 0.6, srt = 90)

dev.off()