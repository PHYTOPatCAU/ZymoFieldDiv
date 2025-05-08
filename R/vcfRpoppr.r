library(vcfR)
library(ade4)
library(ape)
library(adegenet)
library(poppr)
library(seqinr)
library(pegas)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(geosphere)
library(vegan)
# All are output variant call files from 10.5.downtreamfilters.sh
#Beacuse the metadata was slightly different for the fields we accomodated the category mapping for each one
# Define the list of VCF files
vcf_files <- c("~/Zymoproj/merged/Zt_UK.relaxed.mac1.recode.vcf.gz",
               "~/Zymoproj/merged/Zt_UK.max-m-80.biallelic-only.mac1.recode.vcf.gz")

# Loop through each VCF file
for (vcf_file_path in vcf_files) {
  # Load the VCF file
  vcf_file <- read.vcfR(vcf_file_path)
  
  # Transform the vcfR object into a genind/genlight
  glUK <- vcfR2genlight(vcf_file)
  ploidy(glUK) <- 1
  pop <- read.csv2(file = "site_UK_all_SNP_call.csv", header = TRUE)
  pop(glUK) <- pop$pop
  UK_sc <- as.snpclone(glUK)
  
  # Gen dist
  dist <- poppr::bitwise.dist(UK_sc)
  
  # Hierarchical clustering
  hc <- hclust(as.dist(dist))
  
  # Clone definition Singh, Karisto, Croll 2021
  clusters <- cutree(hc, h = 0.01)
  
  # Count the number of samples in each cluster
  cluster_counts <- table(clusters)
  
  # Count the overall samples below the threshold
  below_threshold_count <- sum(cluster_counts[cluster_counts >= 2])
  
  # Extract plant and leaf information from sample names
  sample_info <- strsplit(names(clusters), "-")
  plant_info <- sapply(sample_info, `[`, 2)
  leaf_info <- sapply(sample_info, `[`, 3)
  
  # Data frame with sample names, cluster assignments, and plant and leaf info
  data_clusters <- data.frame(
    sample = names(clusters),
    cluster = clusters,
    plant = plant_info,
    leaf = leaf_info
  )
  
  # Sort the data by plant and leaf
  data_clusters <- data_clusters %>%
    arrange(plant, leaf)
  
  # Add cluster counts to the data_clusters data frame
  data_clusters <- data_clusters %>%
    group_by(cluster) %>%
    mutate(cluster_size = n())
  
  # Filter clusters that have more than one sample
  clusters_with_multiple_samples <- data_clusters %>%
    filter(cluster_size > 1)
  
  # Print the number of clusters with more than one sample
  num_clusters_with_multiple_samples <- clusters_with_multiple_samples %>%
    distinct(cluster) %>%
    nrow()
  print(num_clusters_with_multiple_samples)
  
  # Export the clusters with more than one sample to a CSV file
  write.csv(clusters_with_multiple_samples, file = paste0("cluster_assignments_clonal_groups_", basename(vcf_file_path), ".csv"), row.names = FALSE)
  
  # Export the data_clusters data frame to a CSV file
  write.csv(data_clusters, file = paste0("cluster_assignments_", basename(vcf_file_path), ".csv"), row.names = FALSE)
  
  # Df for same plant/same leaf
  cluster_info <- data_clusters %>%
    group_by(cluster) %>%
    summarize(same_plant = length(unique(plant)) == 1,
              same_leaf = length(unique(leaf)) == 1,
              count_same_plant = sum(duplicated(plant) | duplicated(plant, fromLast = TRUE)),
              count_same_leaf = sum(duplicated(paste(plant, leaf)) | duplicated(paste(plant, leaf), fromLast = TRUE)))
  
  # Convert the dist matrix into a long format data frame
  dist_df <- melt(as.matrix(dist))
  
  # Create a mapping from sample names to categories
  category_mapping <- data_clusters %>%
    left_join(cluster_info, by = "cluster") %>%
    select(sample, same_plant, same_leaf) %>%
    mutate(category = ifelse(same_plant, "same_plant",
                             ifelse(same_leaf, "same_leaf", "same_field")))
  
  # Add a new column to dist_df that indicates the category of each pairwise comparison
  dist_df <- dist_df %>%
    left_join(category_mapping, by = c("Var1" = "sample")) %>%
    rename(category1 = category) %>%
    left_join(category_mapping, by = c("Var2" = "sample")) %>%
    rename(category2 = category) %>%
    mutate(category = ifelse(category1 == category2, category1, "same_field"))
  
  # Adjust the category assignment so that same_field only includes comparisons that are not in same_plant or same_leaf
  dist_df$category[dist_df$category == "same_field" & (dist_df$category1 == "same_plant" | dist_df$category2 == "same_plant")] <- "same_field"
  dist_df$category[dist_df$category == "same_field" & (dist_df$category1 == "same_leaf" | dist_df$category2 == "same_leaf")] <- "same_field"
  
  # Save the violin plot
  svg(file = paste0("Distance_violinplot_", basename(vcf_file_path), ".svg"))
  ggplot(dist_df, aes(x = category, y = value, fill = category)) +
    geom_violin() +
    theme_bw() +
    labs(x = "Category", y = "Pairwise Genetic Distance") +
    scale_fill_brewer(palette = "Set1")
  dev.off()
  
  # NTJ
  # Identify clusters with 2 or more samples
  cluster_sizes <- table(clusters)
  large_clusters <- names(cluster_sizes[cluster_sizes >= 2])
  
  # Convert to phylogenetic tree
  phylo_tree <- as.phylo(hc)
  
  svg(file = paste0("UK_NJTclust01_HC_circular_", basename(vcf_file_path), ".svg"), width = 20, height = 16)
  plot(phylo_tree, type = "fan", show.tip.label = FALSE)
  
  # Highlight large clusters
  for (cluster in large_clusters) {
    tips_in_cluster <- which(clusters == cluster)
    mrca_node <- getMRCA(phylo_tree, tips_in_cluster)
    nodelabels(node = mrca_node, pch = 21, col = "red", bg = "yellow", cex = 1.5)
  }
  
  # Max-min-mean values
  genetic_dist <- mean(dist)
  print(genetic_dist)
  
  min_gen_dist <- min(dist)
  print(min_gen_dist)
  
  max_gen_dist <- max(dist)
  print(max_gen_dist)
  
  # Histogram for gen dist
  upper_tri <- dist[upper.tri(dist)]
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100) +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.3))
  svg(file = paste0("Hist_genetic_dist_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
  
  # Histogram normalized by number of pairwise comparisons
  distance_matrix <- as.matrix(dist)
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100, stat = "density") +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.15)) +
    scale_y_continuous(limits = c(0, 170)) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "skyblue")
  svg(file = paste0("Hist_genetic_dist_norm_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
  
  # Mantel test
  geo <- data.frame(pop$long, pop$lat)
  d.geo <- distm(geo, fun = distHaversine)
  dist.geo <- as.dist(d.geo)
  class(dist.geo)
  genvsdist <- mantel(dist, dist.geo, method = "spearman", permutations = 1000, na.rm = TRUE)
  print(genvsdist)
  pdf(file = paste0("Rose_mantel_", basename(vcf_file_path), ".pdf"))
  plot((dist.geo), (dist), pch = 20, cex = 1, xlab = "Geographic Distance in m", ylab = "Genetic Distance", col = "lightblue")
  dev.off()
  
  # Convert genetic distance matrix to a matrix
  dist_matrix <- as.matrix(dist)
  
  # Ensure sample names are consistent
  rownames(dist.geo) <- colnames(dist.geo) <- rownames(dist_matrix) <- colnames(dist_matrix) <- names(clusters)
  
  # Extract the upper triangular part of the genetic distance matrix
  upper_tri_genetic <- as.data.frame(as.table(dist_matrix))
  upper_tri_genetic <- upper_tri_genetic[upper_tri_genetic$Var1 != upper_tri_genetic$Var2, ]
  upper_tri_genetic <- upper_tri_genetic %>%
    rename(Sample1 = Var1, Sample2 = Var2, GeneticDistance = Freq)
  
  # Extract the upper triangular part of the geographic distance matrix
  upper_tri_geo <- as.data.frame(as.table(dist.geo))
  upper_tri_geo <- upper_tri_geo[upper_tri_geo$Var1 != upper_tri_geo$Var2, ]
  upper_tri_geo <- upper_tri_geo %>%
    rename(Sample1 = Var1, Sample2 = Var2, GeographicDistance = Freq)
  
  # Merge the two data frames on the sample pairs
  merged_distances <- merge(upper_tri_genetic, upper_tri_geo, by = c("Sample1", "Sample2"))
  
  # Filter the merged data frame to include only clusters with more than one sample
  merged_distances <- merged_distances %>%
    filter(Sample1 %in% clusters_with_multiple_samples$sample & Sample2 %in% clusters_with_multiple_samples$sample)
  
  # Print the merged data frame
  head(merged_distances)
  
  # Ensure sample names are consistent in clusters_with_multiple_samples
  clusters_with_multiple_samples <- clusters_with_multiple_samples %>%
    mutate(sample = as.character(sample))
  
  # Generate all possible pairs within each cluster
  sample_pairs <- clusters_with_multiple_samples %>%
    group_by(cluster) %>%
    summarise(pairs = list(combn(sample, 2, simplify = FALSE))) %>%
    unnest(pairs) %>%
    mutate(Sample1 = map_chr(pairs, 1),
           Sample2 = map_chr(pairs, 2)) %>%
    select(Sample1, Sample2)
  
  # Use inner join to filter merged_distances to include only these pairs
  filtered_sample_pairs <- merged_distances %>%
    inner_join(sample_pairs, by = c("Sample1", "Sample2"))
  
  # Print the filtered sample pairs with geographic distances
  print(filtered_sample_pairs)
  
  # Calculate the median of the GeographicDistance
  median_geographic_distance <- median(filtered_sample_pairs$GeographicDistance, na.rm = TRUE)
  
  # Calculate the mean of the GeographicDistance
  mean_geographic_distance <- mean(filtered_sample_pairs$GeographicDistance, na.rm = TRUE)
  
  # Print the results
  print(paste("Median Geographic Distance:", median_geographic_distance))
  print(paste("Mean Geographic Distance:", mean_geographic_distance))
  
  # Export the filtered sample pairs with geographic distances to a CSV file
  write.csv(filtered_sample_pairs, file = paste0("cluster_pairs_with_geographic_distances_", basename(vcf_file_path), ".csv"), row.names = FALSE)
  
  # Export the merged data frame to a CSV file
  write.csv(merged_distances, file = paste0("merged_distances_", basename(vcf_file_path), ".csv"), row.names = FALSE)
  
  # Read the CSV file into a data frame
  clonal_data <- read.csv("clonal_dist_gps_UK.csv")
  
  # Initialize an empty list to store distances
  distances <- list()
  
  # Loop through the data frame in steps of 2 to get consecutive pairs
  for (i in seq(1, nrow(clonal_data), by = 2)) {
    # Extract the coordinates of each pair
    coord1 <- c(clonal_data$long[i], clonal_data$lat[i])
    coord2 <- c(clonal_data$long[i + 1], clonal_data$lat[i + 1])
    
    # Calculate the distance using distHaversine
    distance <- distHaversine(coord1, coord2)
    
    # Store the result in the list
    distances[[i]] <- data.frame(
      sample1 = clonal_data$sample.id[i],
      sample2 = clonal_data$sample.id[i + 1],
      distance = distance
    )
  }
  
  # Convert the list to a data frame
  distances_df <- do.call(rbind, distances)
  
  # Print the results
  print(distances_df)
  
  # Export the table to a text file
  write.table(distances_df, file = paste0("distances_clonal_pairs_diff_gps_", basename(vcf_file_path), ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
}
###############################################################################################
#CH
# Define the list of VCF files
vcf_files <- c("~/Zymoproj/merged/Zt_CH.relaxed.mac1.recode.vcf.gz",
               "~/Zymoproj/merged/Zt_CH.max-m-80.biallelic-only.mac1.recode.vcf.gz")

# Loop through each VCF file
for (vcf_file_path in vcf_files) {
  # Load the VCF file
  vcf_file <- read.vcfR(vcf_file_path)
  
  # Transform the vcfR object into a genind/genlight
  glCH <- vcfR2genlight(vcf_file)
  ploidy(glCH) <- 1
  pop <- read.csv2(file = "site_CH_SNP_call2.csv", header = TRUE)
  pop(glCH) <- pop$pop
  CH_sc <- as.snpclone(glCH)
  
  # Gen dist
  dist <- poppr::bitwise.dist(CH_sc)
  
  # Hierarchical clustering
  hc <- hclust(as.dist(dist))
  
  # Clone definition Singh, Karisto, Croll 2021
  clusters <- cutree(hc, h = 0.01)
  
  # Count the number of samples in each cluster
  cluster_counts <- table(clusters)
  
  # Count the overall samples below the threshold
  below_threshold_count <- sum(cluster_counts[cluster_counts >= 2])
  
  # Extract plot information from sample names
  plot_info <- sapply(strsplit(names(clusters), "_"), function(x) substr(x[2], 2, 2))
  
  # Data frame with sample names, cluster assignments, and plot info
  data_clusters <- data.frame(
    sample = names(clusters),
    cluster = clusters,
    plot = plot_info
  )
  
  # Sort the data by plot
  data_clusters <- data_clusters %>%
    arrange(plot)
  
  # Df for same plot
  cluster_info <- data_clusters %>%
    group_by(cluster) %>%
    summarize(same_plot = length(unique(plot)) == 1,
              count_same_plot = sum(duplicated(plot) | duplicated(plot, fromLast = TRUE)))
  
  # Create a mapping from sample names to categories
  category_mapping <- data_clusters %>%
    left_join(cluster_info, by = "cluster") %>%
    select(sample, same_plot) %>%
    mutate(category = ifelse(same_plot, "same_plot", "same_field"))
  
  # Convert the dist matrix into a long format data frame
  dist_df <- melt(as.matrix(dist))
  
  # Add a new column to dist_df that indicates the category of each pairwise comparison
  dist_df <- dist_df %>%
    left_join(category_mapping, by = c("Var1" = "sample")) %>%
    rename(category1 = category) %>%
    left_join(category_mapping, by = c("Var2" = "sample")) %>%
    rename(category2 = category) %>%
    mutate(category = ifelse(category1 == category2, category1, "same_field"))
  
  # Adjust the category assignment so that same_field only includes comparisons that are not in same_plot
  dist_df$category[dist_df$category == "same_field" & dist_df$category1 == "same_plot"] <- "same_field"
  
  # Save the violin plot
  svg(file = paste0("Distance_violinplot_", basename(vcf_file_path), ".svg"))
  ggplot(dist_df, aes(x = category, y = value, fill = category)) +
    geom_violin() +
    theme_bw() +
    labs(x = "Category", y = "Pairwise Genetic Distance") +
    scale_fill_brewer(palette = "Set1")
  dev.off()
  
  # NTJ
  # Identify clusters with 2 or more samples
  cluster_sizes <- table(clusters)
  large_clusters <- names(cluster_sizes[cluster_sizes >= 2])
  
  # Convert to phylogenetic tree
  phylo_tree <- as.phylo(hc)
  
  svg(file = paste0("CH_NJTclust01_HC_circular_", basename(vcf_file_path), ".svg"), width = 20, height = 16)
  plot(phylo_tree, type = "fan", show.tip.label = FALSE)
  
  # Highlight large clusters
  for (cluster in large_clusters) {
    tips_in_cluster <- which(clusters == cluster)
    mrca_node <- getMRCA(phylo_tree, tips_in_cluster)
    nodelabels(node = mrca_node, pch = 21, col = "red", bg = "yellow", cex = 1.5)
  }
  
  # Max-min-mean values
  genetic_dist <- mean(dist)
  print(genetic_dist)
  
  min_gen_dist <- min(dist)
  print(min_gen_dist)
  
  max_gen_dist <- max(dist)
  print(max_gen_dist)
  
  # Histogram for gen dist
  upper_tri <- dist[upper.tri(dist)]
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100) +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.3))
  svg(file = paste0("Hist_genetic_dist_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
  
  # Histogram normalized by number of pairwise comparisons
  distance_matrix <- as.matrix(dist)
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100, stat = "density") +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.15)) +
    scale_y_continuous(limits = c(0, 170)) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "skyblue")
  svg(file = paste0("Hist_genetic_dist_norm_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
}

###############################################################################################
#US
# Define the list of VCF files
vcf_files <- c("~/Zymoproj/merged/Zt_US.relaxed.mac1.recode.vcf.gz",
               "~/Zymoproj/merged/Zt_US.max-m-80.biallelic-only.mac1.recode.vcf.gz")

# Loop through each VCF file
for (vcf_file_path in vcf_files) {
  # Load the VCF file
  vcf_file <- read.vcfR(vcf_file_path)
  
  # Transform the vcfR object into a genind/genlight
  glUS <- vcfR2genlight(vcf_file)
  ploidy(glUS) <- 1
  pop <- read.csv2(file = "site_US_SNP_call.csv", header = TRUE)
  pop(glUS) <- pop$pop
  US_sc <- as.snpclone(glUS)
  
  # Gen dist
  dist <- poppr::bitwise.dist(US_sc)
  
  # Hierarchical clustering
  hc <- hclust(as.dist(dist))
  
  # Clone definition Singh, Karisto, Croll 2021
  clusters <- cutree(hc, h = 0.01)
  
  # Count the number of samples in each cluster
  cluster_counts <- table(clusters)
  
  # Count the overall samples below the threshold
  below_threshold_count <- sum(cluster_counts[cluster_counts >= 2])
  
  # Extract plot information from sample names
  plot_info <- sapply(strsplit(names(clusters), "_"), function(x) gsub(".*_([A-Z][0-9])$", "\\1", x[3]))
  
  # Data frame with sample names, cluster assignments, and plot info
  data_clusters <- data.frame(
    sample = names(clusters),
    cluster = clusters,
    plot = plot_info
  )
  
  # Sort the data by plot
  data_clusters <- data_clusters %>%
    arrange(plot)
  
  # Df for same plot
  cluster_info <- data_clusters %>%
    group_by(cluster) %>%
    summarize(same_plot = length(unique(plot)) == 1,
              count_same_plot = sum(duplicated(plot) | duplicated(plot, fromLast = TRUE)))
  
  # Create a mapping from sample names to categories
  category_mapping <- data_clusters %>%
    left_join(cluster_info, by = "cluster") %>%
    select(sample, same_plot) %>%
    mutate(category = ifelse(same_plot, "same_plot", "same_field"))
  
  # Convert the dist matrix into a long format data frame
  dist_df <- melt(as.matrix(dist))
  
  # Add a new column to dist_df that indicates the category of each pairwise comparison
  dist_df <- dist_df %>%
    left_join(category_mapping, by = c("Var1" = "sample")) %>%
    rename(category1 = category) %>%
    left_join(category_mapping, by = c("Var2" = "sample")) %>%
    rename(category2 = category) %>%
    mutate(category = ifelse(category1 == category2, category1, "same_field"))
  
  # Adjust the category assignment so that same_field only includes comparisons that are not in same_plot
  dist_df$category[dist_df$category == "same_field" & dist_df$category1 == "same_plot"] <- "same_field"
  
  # Save the violin plot
  svg(file = paste0("Distance_violinplot_", basename(vcf_file_path), ".svg"))
  ggplot(dist_df, aes(x = category, y = value, fill = category)) +
    geom_violin() +
    theme_bw() +
    labs(x = "Category", y = "Pairwise Genetic Distance") +
    scale_fill_brewer(palette = "Set1")
  dev.off()
  
  # NTJ
  # Identify clusters with 2 or more samples
  cluster_sizes <- table(clusters)
  large_clusters <- names(cluster_sizes[cluster_sizes >= 2])
  
  # Convert to phylogenetic tree
  phylo_tree <- as.phylo(hc)
  
  svg(file = paste0("US_NJTclust01_HC_circular_", basename(vcf_file_path), ".svg"), width = 20, height = 16)
  plot(phylo_tree, type = "fan", show.tip.label = FALSE)
  
  # Highlight large clusters
  for (cluster in large_clusters) {
    tips_in_cluster <- which(clusters == cluster)
    mrca_node <- getMRCA(phylo_tree, tips_in_cluster)
    nodelabels(node = mrca_node, pch = 21, col = "red", bg = "yellow", cex = 1.5)
  }
  
  # Max-min-mean values
  genetic_dist <- mean(dist)
  print(genetic_dist)
  
  min_gen_dist <- min(dist)
  print(min_gen_dist)
  
  max_gen_dist <- max(dist)
  print(max_gen_dist)
  
  # Histogram for gen dist
  upper_tri <- dist[upper.tri(dist)]
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100) +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.3))
  svg(file = paste0("Hist_genetic_dist_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
  
  # Histogram normalized by number of pairwise comparisons
  distance_matrix <- as.matrix(dist)
  distribution <- ggplot(data.frame(distances = upper_tri), aes(x = distances)) +
    geom_histogram(bins = 100, stat = "density") +
    labs(x = "Distance", y = "Frequency") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 0.15)) +
    scale_y_continuous(limits = c(0, 170)) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "skyblue")
  svg(file = paste0("Hist_genetic_dist_norm_", basename(vcf_file_path), ".svg"))
  distribution
  dev.off()
}





