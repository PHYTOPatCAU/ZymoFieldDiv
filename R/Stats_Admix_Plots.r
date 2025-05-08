#This is a collection of miniscripts used to plot the data from various scripts: 04.5,13,15,17,18,20
library(ggplot2)
library(ggpattern)
library(tidyverse)
library(reshape2)
library(scales)
library(data.table)
library(multcompView)
#set wd
setwd("~/Zymoproj/merged/")

#Quality check for depth per chromosome
# Define the directories
dirs <- list(
  UK = "~/Zymoproj/UKsamples",
  US = "~/Zymoproj/USsamples",
  CH = "~/Zymoproj/CHsamples"
)

# Define the output directory
output_dir <- "~/Zymoproj/merged"

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
# Define the regions (countries) and types (strict and relaxed)
regions <- c("UK", "US", "CH")
types <- c("relaxed", "strict")

# Function to process .imiss files and create a genotyping rate plot
process_imiss_file <- function(imiss_file, region, type) {
  # Load the .imiss file
  data <- read.delim(imiss_file, sep = "", header = TRUE)
  
  # Calculate the genotyping rate
  data <- data %>%
    mutate(GenotypingRate = N_GENO / (N_GENO + N_MISS) * 100)
  
  # Create the genotyping rate plot
  plot <- ggplot(data, aes(x = reorder(1:nrow(data), -GenotypingRate), y = GenotypingRate)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Genotyping Rate -", region, type),
         x = "Samples",
         y = "Genotyping Rate (%)") +
    theme_bw() +
    geom_hline(yintercept = mean(data$GenotypingRate), color = "skyblue", linetype = "dashed")
  
  # Save the plot
  svg(file = paste0(region, "_GT_", type, ".svg"))
  print(plot)
  dev.off()
}

# Function to process .frq files and create MAF plots
process_frq_file <- function(frq_file, region, type) {
  # Load the .frq file
  maf_data <- fread(frq_file)
  
  # Create the raw MAF plot
  raw_maf_plot <- ggplot(maf_data, aes(x = MAF)) +
    geom_histogram(bins = 20, fill = "gray", color = "black") +
    labs(title = paste("MAF Distribution -", region, type),
         x = "MAF",
         y = "SNP Count") +
    theme_bw()
  
  # Save the raw MAF plot
  svg(file = paste0(region, "_MAF_", type, "_raw.svg"))
  print(raw_maf_plot)
  dev.off()
  
  # Normalize the MAF data by sample size
  maf_data <- maf_data %>%
    mutate(NormCount = NCHROBS / sum(NCHROBS))
  
  # Create the normalized MAF plot
  norm_maf_plot <- ggplot(maf_data, aes(x = MAF, y = NormCount)) +
    geom_histogram(stat = "identity", bins = 20, fill = "darkgray", color = "black") +
    labs(title = paste("Normalized MAF Distribution -", region, type),
         x = "MAF",
         y = "Normalized SNP Count") +
    theme_bw()
  
  # Save the normalized MAF plot
  svg(file = paste0(region, "_MAF_", type, "_norm.svg"))
  print(norm_maf_plot)
  dev.off()
}

# Function to process .lmiss files and create a missingness per site plot
process_lmiss_file <- function(lmiss_file, region, type) {
  # Load the .lmiss file
  lmiss_data <- read.delim(lmiss_file, sep = "", header = TRUE)
  
  # Calculate the missingness rate for each site
  lmiss_data <- lmiss_data %>%
    mutate(MissingnessRate = F_MISS * 100)
  
  # Create a histogram for missingness rate
  plot <- ggplot(lmiss_data, aes(x = MissingnessRate)) +
    geom_histogram(bins = 30, fill = ifelse(type == "relaxed", "lightgrey", "darkgrey"), color = "black") +
    labs(title = paste("Missingness Rate Distribution -", region, type),
         x = "Missingness Rate (%)",
         y = "Site Count") +
    theme_bw()
  
  # Save the plot
  svg(file = paste0(region, "_miss_per_site_", type, ".svg"))
  print(plot)
  dev.off()
}

# Loop through each region and type
for (region in regions) {
  for (type in types) {
    # Define the file names
    imiss_file <- paste0("Zt_", region, ".", type, ".imiss")
    frq_file <- paste0("Zt_", region, ".", type, ".frq")
    lmiss_file <- paste0("Zt_", region, ".", type, ".lmiss")
    
    # Process the .imiss file
    process_imiss_file(imiss_file, region, type)
    
    # Process the .frq file
    process_frq_file(frq_file, region, type)
    
    # Process the .lmiss file
    process_lmiss_file(lmiss_file, region, type)
  }
}
# Script for processing diversity statistics per field and plotting for Pi and Tajimas_D at Chr level
divstatsall <- read.csv("~/Zymoproj/merged/DivStatsAll.csv")

# Preprocess the data
divstatscorrectedall <- divstatscorrectedall %>%
  mutate(
    file = gsub(".SNP.qualityfilter.excl.rn1.sorted.vcf.gz", "", file),
    chr = ifelse(chr == "mt", "mt", paste0("chr", chr)),
    country = case_when(
      grepl("Zt_US", file) ~ "US",
      grepl("Zt_UK", file) ~ "UK",
      grepl("Zt_CH2", file) ~ "CH",
      TRUE ~ NA_character_
    ),
    Tajimas_D = as.numeric(Tajimas_D),
    Pi = as.numeric(Pi),
    Wattersons_Theta = as.numeric(Wattersons_Theta)
  ) %>%
  filter(
    !is.na(Tajimas_D) & !is.na(Pi) & !is.na(Wattersons_Theta),
    Tajimas_D != "indet" & Pi != "indet" & Wattersons_Theta != "indet"
  )

# Define core and accessory chromosomes
core_chromosomes <- paste0("chr", 1:13)
accessory_chromosomes <- paste0("chr", 14:21)

# Filter data for core and accessory chromosomes
core_data <- divstatscorrectedall %>% filter(chr %in% core_chromosomes)
accessory_data <- divstatscorrectedall %>% filter(chr %in% accessory_chromosomes)

# Function to perform Kruskal-Wallis and Wilcoxon rank-sum tests
perform_tests <- function(data, metric) {
  kruskal_result <- kruskal.test(as.formula(paste(metric, "~ country")), data = data)
  pairwise_result <- pairwise.wilcox.test(data[[metric]], data$country, p.adjust.method = "bonferroni")
  list(kruskal = kruskal_result, pairwise = pairwise_result)
}

# Perform tests for core and accessory chromosomes
core_tests <- list(
  Tajimas_D = perform_tests(core_data, "Tajimas_D"),
  Pi = perform_tests(core_data, "Pi"),
  Wattersons_Theta = perform_tests(core_data, "Wattersons_Theta")
)

accessory_tests <- list(
  Tajimas_D = perform_tests(accessory_data, "Tajimas_D"),
  Pi = perform_tests(accessory_data, "Pi"),
  Wattersons_Theta = perform_tests(accessory_data, "Wattersons_Theta")
)

# Function to perform permutation tests
permutation_test_all_pairs <- function(data, metric, group_col, n_permutations = 1000) {
  groups <- unique(data[[group_col]])
  p_values <- matrix(NA, nrow = length(groups), ncol = length(groups), dimnames = list(groups, groups))
  
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group1 <- groups[i]
      group2 <- groups[j]
      observed_diff <- abs(mean(data[[metric]][data[[group_col]] == group1], na.rm = TRUE) - 
                           mean(data[[metric]][data[[group_col]] == group2], na.rm = TRUE))
      perm_diffs <- replicate(n_permutations, {
        permuted_data <- data
        permuted_data[[group_col]] <- sample(permuted_data[[group_col]])
        abs(mean(permuted_data[[metric]][permuted_data[[group_col]] == group1], na.rm = TRUE) - 
            mean(permuted_data[[metric]][permuted_data[[group_col]] == group2], na.rm = TRUE))
      })
      p_value <- mean(perm_diffs >= observed_diff)
      p_values[group1, group2] <- p_value
      p_values[group2, group1] <- p_value
    }
  }
  return(p_values)
}

# Perform permutation tests for core and accessory chromosomes
core_permutation_tests <- list(
  Tajimas_D = permutation_test_all_pairs(core_data, "Tajimas_D", "country"),
  Pi = permutation_test_all_pairs(core_data, "Pi", "country"),
  Wattersons_Theta = permutation_test_all_pairs(core_data, "Wattersons_Theta", "country")
)

accessory_permutation_tests <- list(
  Tajimas_D = permutation_test_all_pairs(accessory_data, "Tajimas_D", "country"),
  Pi = permutation_test_all_pairs(accessory_data, "Pi", "country"),
  Wattersons_Theta = permutation_test_all_pairs(accessory_data, "Wattersons_Theta", "country")
)

# Function to create boxplots with significance stars
plot_metric <- function(data, metric, title, p_values) {
  get_stars <- function(p_value) {
    if (p_value < 0.001) {
      return("***")
    } else if (p_value < 0.01) {
      return("**")
    } else if (p_value < 0.05) {
      return("*")
    } else {
      return("")
    }
  }
  
  plot <- ggplot(data, aes(x = country, y = !!sym(metric), fill = country)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = title, x = "Country", y = metric) +
    theme(text = element_text(size = 18),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 16),
          strip.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16))
  
  comparisons <- combn(unique(data$country), 2, simplify = FALSE)
  for (comparison in comparisons) {
    group1 <- comparison[1]
    group2 <- comparison[2]
    p_value <- p_values[group1, group2]
    if (!is.na(p_value) && p_value < 0.05) {
      plot <- plot + geom_signif(comparisons = list(comparison),
                                 annotations = get_stars(p_value),
                                 y_position = max(data[[metric]], na.rm = TRUE) + 0.1,
                                 tip_length = 0.01, size = 1.5)
    }
  }
  return(plot)
}

# Create and save plots for core and accessory chromosomes
metrics <- c("Tajimas_D", "Pi", "Wattersons_Theta")
titles <- c("Tajima's D", "Pi", "Watterson's Theta")

for (i in seq_along(metrics)) {
  metric <- metrics[i]
  title <- titles[i]
  
  core_plot <- plot_metric(core_data, metric, paste("Core Chromosomes -", title), core_permutation_tests[[metric]])
  accessory_plot <- plot_metric(accessory_data, metric, paste("Accessory Chromosomes -", title), accessory_permutation_tests[[metric]])
  
  svg(paste0(metric, "_core_accessory.svg"), width = 16, height = 10)
  grid.arrange(core_plot, accessory_plot, ncol = 1, nrow = 2)
  dev.off()
}

# Export Kruskal-Wallis and permutation test results
kruskal_results <- data.frame(
  Metric = rep(metrics, each = 2),
  Chromosome_Type = rep(c("Core", "Accessory"), times = length(metrics)),
  P_Value = c(
    sapply(core_tests, function(x) x$kruskal$p.value),
    sapply(accessory_tests, function(x) x$kruskal$p.value)
  )
)
write.table(kruskal_results, "kruskal_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

summary_table <- data.frame(
  Metric = rep(metrics, each = 3),
  Comparison = rep(c("CH vs UK", "CH vs US", "UK vs US"), times = length(metrics)),
  Core_P_Value = unlist(lapply(core_permutation_tests, function(x) c(x["CH", "UK"], x["CH", "US"], x["UK", "US"]))),
  Accessory_P_Value = unlist(lapply(accessory_permutation_tests, function(x) c(x["CH", "UK"], x["CH", "US"], x["UK", "US"])))
)
write.table(summary_table, "permutation_test_summary.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#load data for Pi and Tajimas_D per Effector
# Load data
divstatsmac1 <- read.csv("~/Zymoproj/merged/divstats_all_eff_regions.csv")

# Preprocess the data
divstatsmac1 <- divstatsmac1 %>%
  rename(Effector = Effector, Pi = Pi, Wattersons_Theta = Wattersons_Theta, Tajima_D = Tajimas_D, region = country) %>%
  mutate(
    Effector = str_remove(Effector, ".vcf.gz$"),
    country = str_remove(country, "^Zt_"),
    Pi = as.numeric(Pi),
    Tajima_D = as.numeric(Tajima_D),
    Wattersons_Theta = as.numeric(Wattersons_Theta)
  ) %>%
  filter(!is.na(Pi) & !is.na(Tajima_D) & !is.na(Wattersons_Theta))

# Reorder Effector based on Pi values
divstatsmac1 <- divstatsmac1 %>%
  arrange(Pi) %>%
  mutate(Effector = factor(Effector, levels = unique(Effector)))

# Perform ANOVA and Tukey's HSD test for Pi, Tajima's D, and Watterson's Theta
perform_anova_tukey <- function(data, metric) {
  anova_result <- aov(as.formula(paste(metric, "~ country")), data = data)
  tukey_result <- TukeyHSD(anova_result)
  list(anova = summary(anova_result), tukey = tukey_result)
}

anova_tukey_results <- list(
  Pi = perform_anova_tukey(divstatsmac1, "Pi"),
  Tajima_D = perform_anova_tukey(divstatsmac1, "Tajima_D"),
  Wattersons_Theta = perform_anova_tukey(divstatsmac1, "Wattersons_Theta")
)

# Print ANOVA and Tukey's HSD results
for (metric in names(anova_tukey_results)) {
  cat("\n### ANOVA Results for", metric, "###\n")
  print(anova_tukey_results[[metric]]$anova)
  
  cat("\n### Tukey's HSD Results for", metric, "###\n")
  print(anova_tukey_results[[metric]]$tukey)
}

# Save Tukey's HSD results to text files
for (metric in names(anova_tukey_results)) {
  write.table(
    as.data.frame(anova_tukey_results[[metric]]$tukey),
    file = paste0(metric, "_tukey_results.txt"),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE
  )
}

# Calculate the 5% extreme values for Tajima's D
tajima_d_limits <- quantile(divstatsmac1$Tajima_D, probs = c(0.025, 0.975), na.rm = TRUE)
lower_limit <- tajima_d_limits[1]
upper_limit <- tajima_d_limits[2]

cat("Lower 2.5% limit for Tajima's D:", lower_limit, "\n")
cat("Upper 2.5% limit for Tajima's D:", upper_limit, "\n")

# Plotting Pi
pi_histogram <- ggplot(divstatsmac1, aes(x = Effector, y = Pi, fill = country)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  ylim(0, 1) +
  labs(x = "Effector", y = "Pi")

svg("Pi_histogram_eff2.svg", width = 9, height = 6.5)
print(pi_histogram)
dev.off()

# Plotting Tajima's D with significance limits
tajima_histogram <- ggplot(divstatsmac1, aes(x = Effector, y = Tajima_D, fill = country)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  scale_y_continuous(breaks = seq(-9, 9, by = 3)) +
  geom_hline(yintercept = lower_limit, color = "blue", linetype = "dashed", size = 1) +
  geom_hline(yintercept = upper_limit, color = "blue", linetype = "dashed", size = 1) +
  labs(x = "Effector", y = "Tajima's D", title = "Tajima's D Distribution with Significance Limits")

svg("Tajima_D_histogram_eff2.svg", width = 9, height = 6.5)
print(tajima_histogram)
dev.off()

# Boxplot for Tajima's D with significance limits
box_plot_tajima_d <- ggplot(divstatsmac1, aes(x = country, y = Tajima_D, fill = country)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Country", y = "Tajima's D", fill = "Field") +
  theme(text = element_text(size = 18),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  geom_hline(yintercept = lower_limit, color = "blue", linetype = "dashed", size = 1) +
  geom_hline(yintercept = upper_limit, color = "blue", linetype = "dashed", size = 1)

svg("Tajima_D_boxplot_eff2.svg", width = 10, height = 8)
print(box_plot_tajima_d)
dev.off()

# Generate compact letter display (CLD) for Tukey's HSD results
plot_tukey_cld <- function(data, metric, tukey_result, title) {
  cld <- multcompLetters4(aov(as.formula(paste(metric, "~ country")), data = data), tukey_result[[metric]]$tukey)
  cld_df <- data.frame(
    country = names(cld$Letters),
    Letters = cld$Letters
  )
  
  plot <- ggplot(data, aes(x = country, y = !!sym(metric), fill = country)) +
    geom_boxplot() +
    geom_text(data = cld_df, aes(x = country, y = max(data[[metric]], na.rm = TRUE) + 0.1, label = Letters), size = 6) +
    labs(title = title, x = "Country", y = metric) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "none",
      text = element_text(size = 16),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
  
  return(plot)
}

# Generate and save plots for each metric with CLD
for (metric in names(anova_tukey_results)) {
  plot <- plot_tukey_cld(divstatsmac1, metric, anova_tukey_results, paste("Tukey's HSD for", metric))
  
  svg(paste0(metric, "_tukey_cld_plot.svg"), width = 10, height = 8)
  print(plot)
  dev.off()
}
# Define a permutation function
perm_test <- function(data, metric, n_perm = 1000) {
  observed_diff <- mean(data[[metric]][data$Genome == "Core"], na.rm = TRUE) - 
    mean(data[[metric]][data$Genome == "Accessory"], na.rm = TRUE)
  
  perm_diffs <- replicate(n_perm, {
    permuted_genome <- sample(data$Genome)
    mean(data[[metric]][permuted_genome == "Core"], na.rm = TRUE) - 
      mean(data[[metric]][permuted_genome == "Accessory"], na.rm = TRUE)
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  list(observed_diff = observed_diff, p_value = p_value)
}

# Run permutation tests per field
fields <- unique(divstatsmac1$Field)
results_list <- list()

for (field in fields) {
  field_data <- divstatsmac1 %>% filter(Field == field)
  perm_test_Pi <- perm_test(field_data, "Pi")
  perm_test_TajimaD <- perm_test(field_data, "Tajima_D")
  perm_test_WattersonTheta <- perm_test(field_data, "Wattersons_Theta")
  
  results_list[[field]] <- data.frame(
    Field = field,
    Metric = c("Pi", "Tajima's D", "Watterson's Theta"),
    Observed_Difference = c(perm_test_Pi$observed_diff, perm_test_TajimaD$observed_diff, perm_test_WattersonTheta$observed_diff),
    p_value = c(perm_test_Pi$p_value, perm_test_TajimaD$p_value, perm_test_WattersonTheta$p_value)
  )
}

results_table <- do.call(rbind, results_list)
print(results_table)

# Define colors
colors <- brewer.pal(3, "Set1")
names(colors) <- c("CH", "UK", "US")

# Visualise the data with improved plots using facet_wrap
Pi_subgenomeplot <- ggplot(divstatsmac1, aes(x = Genome, y = Pi, fill = Field, alpha = Genome)) +
  geom_boxplot() +
  labs(title = "Pi Comparison between Core and Accessory Genes", x = "Genome", y = "Pi") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = colors) +
  scale_alpha_manual(values = c("Core" = 1, "Accessory" = 0.5)) +
  facet_wrap(~ Field) +
  theme(legend.position = "none", 
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  geom_text(data = results_table %>% filter(Metric == "Pi"), aes(x = 1.5, y = max(divstatsmac1$Pi, na.rm = TRUE), label = paste("p-value:", round(p_value, 3))), size = 6, hjust = 0.5, inherit.aes = FALSE)

TD_subgenomeplot <- ggplot(divstatsmac1, aes(x = Genome, y = Tajima_D, fill = Field, alpha = Genome)) +
  geom_boxplot() +
  labs(title = "Tajima's D Comparison between Core and Accessory Genes", x = "Genome", y = "Tajima's D") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = colors) +
  scale_alpha_manual(values = c("Core" = 1, "Accessory" = 0.5)) +
  facet_wrap(~ Field) +
  theme(legend.position = "none", 
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  geom_text(data = results_table %>% filter(Metric == "Tajima's D"), aes(x = 1.5, y = 4, label = paste("p-value:", round(p_value, 3))), size = 6, hjust = 0.5, inherit.aes = FALSE)

WT_subgenomeplot <- ggplot(divstatsmac1_filtered, aes(x = Genome, y = Wattersons_Theta, fill = Field, alpha = Genome)) +
  geom_boxplot() +
  labs(title = "Watterson's Theta Comparison between Core and Accessory Genes", x = "Genome", y = "Watterson's Theta") +
  theme_minimal(base_size = 16) +
  scale_fill_manual(values = colors) +
  scale_alpha_manual(values = c("Core" = 1, "Accessory" = 0.5)) +
  facet_wrap(~ Field) +
  theme(legend.position = "none", 
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)) +
  geom_text(data = results_table %>% filter(Metric == "Watterson's Theta"), aes(x = 1.5, y = 0.3, label = paste("p-value:", round(p_value, 3))), size = 6, hjust = 0.5, inherit.aes = FALSE)


## There is an unaesthetic outlier I'm removing for ease of visualisation (but this does not affect the P value)
divstatsmac1_filtered <- divstatsmac1 %>% filter(Wattersons_Theta <= 0.5)

# Save the plots
svg("Pi_subgenomeplot_facet.svg", width = 10, height = 8)
print(Pi_subgenomeplot)
dev.off()

svg("TD_subgenomeplot_facet.svg", width = 10, height = 8)
print(TD_subgenomeplot)
dev.off()

svg("WT_subgenomeplot_facet.svg", width = 10, height = 8)
print(WT_subgenomeplot)
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

# Add the individual names below the barplot with a smaller font size
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
