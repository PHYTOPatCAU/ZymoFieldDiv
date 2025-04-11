library(pak)
library(Icens)
library(r4pde)
library(knitr)
library(jtools)
library(minpack.lm)
library(ggplot2)
library(tidyverse)
setwd("C:/Users/suaph296/Desktop/Zymo_project/Data/NGS")

# Load the data
df <- read.delim("C:/Users/suaph296/Desktop/Zymo_project/Data/NGS/rarefactioncurveall_corrected.txt")

# Filter the data to include only CH, UK, and US fields
df <- df %>% filter(country %in% c("CH", "UK", "US"))
df$country <- factor(df$country, levels = c("CH", "UK", "US"))
df$country_field <- factor(paste(df$country, df$field, sep = "_"), levels = unique(paste(df$country, df$field, sep = "_")))

# Define custom colors for the selected fields
custom_colors <- c("US.Oregon" = "green", "CH.Eschikon" = "red", "UK.Rosemaund" = "blue")

# Create the plot with custom colors
svg(file = "Rarefaction_field_CH_UK_US.svg")
ggplot(df, aes(x = subset, y = snps / 1e6, color = interaction(country, field), shape = country)) +
  geom_point(stat = "identity") +
  scale_y_continuous(name = "Number of SNPs (millions)") +
  scale_color_manual(values = custom_colors) +
  theme_bw()
dev.off()

# Round the 'subset' column to the nearest double digit
round_to_double_digit <- function(x) {
  return(round(x / 10) * 10)
}
df$subset <- round_to_double_digit(df$subset)

# Create the boxplot with custom colors
svg(file = "Boxplot_field_call.svg")
ggplot(df, aes(x = as.factor(subset), y = snps / 1e6, fill = interaction(country, field))) +
  geom_boxplot() +
  scale_y_continuous(name = "Number of SNPs (millions)") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  xlab("Subset of isolates")
dev.off()

# Prepare the data for modeling
df_filtered <- df
df_filtered$field <- factor(paste(df_filtered$country, df_filtered$field, sep = "_"), levels = unique(paste(df_filtered$country, df_filtered$field, sep = "_")))

# Function to fit the power model and predict SNP counts
fit_and_predict_power_model <- function(data) {
  model <- lm(log(Y) ~ log(x), data = data)
  predict_snp <- function(x) {
    exp(predict(model, newdata = data.frame(x = x)))
  }
  return(list(model = model, predict_snp = predict_snp))
}

# Prepare the data for the test
prepare_data <- function(data) {
  data %>%
    mutate(snp_count_millions = snps) %>%
    select(subset, snp_count_millions, field) %>%
    rename(x = subset, Y = snp_count_millions)
}
grad_data_relaxed <- prepare_data(df_filtered)

# Fit the model for each field and store in a named list
models <- grad_data_relaxed %>%
  group_by(field) %>%
  group_map(~ fit_and_predict_power_model(.x), .keep = TRUE)

# Convert the list to a named list
models <- setNames(models, as.character(unique(grad_data_relaxed$field)))

# Function to find the x value where the increase in SNP count is less than 5% of the maximum SNP count
find_x_for_5_percent <- function(model, max_snp_count) {
  x_values <- seq(10, 10000, by = 10)
  predicted_snp_counts <- model$predict_snp(x_values)
  threshold <- 0.01 * max_snp_count
  x_for_1_percent <- NA
  snp_count_for_1_percent <- NA
  
  for (i in 2:length(predicted_snp_counts)) {
    if ((predicted_snp_counts[i] - predicted_snp_counts[i - 1]) < threshold) {
      x_for_1_percent <- x_values[i]
      snp_count_for_1_percent <- predicted_snp_counts[i]
      break
    }
  }
  
  return(list(
    x_values = x_values,
    predicted_snp_counts = predicted_snp_counts,
    x_for_1_percent = x_for_1_percent,
    snp_count_for_1_percent = snp_count_for_1_percent
  ))
}

# Initialize an empty list to store results
results_list <- list()
predictions_list <- list()

# Loop through each field and calculate the x value for 5% increase
for (field in unique(grad_data_relaxed$field)) {
  max_snp_count <- max(grad_data_relaxed$Y[grad_data_relaxed$field == field])
  result <- find_x_for_1_percent(models[[field]], max_snp_count)
  
  # Debugging statements
  if (is.na(result$x_for_1_percent)) {
    cat("Field:", field, "- No x value found for 1% increase\n")
  } else {
    cat("Field:", field, "- x value for 1% increase:", result$x_for_1_percent, "SNP count:", result$snp_count_for_1_percent, "\n")
  }
  
  results_list[[field]] <- data.frame(
    field = field,
    max_snp_count = max_snp_count,
    x_for_1_percent = result$x_for_1_percent,
    snp_count = result$snp_count_for_1_percent
  )
  
  predictions_list[[field]] <- data.frame(
    field = field,
    x_values = result$x_values,
    predicted_snp_counts = result$predicted_snp_counts
  )
}

# Combine the results into a single data frame
results1 <- do.call(rbind, results_list)
predictions <- do.call(rbind, predictions_list)

print(results5)
# Export the results to a CSV file
write.csv(results1, "samplesizefor1%limit.csv", row.names = FALSE)

# Create the plot with custom colors
svg(file = "Rarefaction__power_model_1percent_cutoff.svg")
modelplot <- ggplot(df_filtered, aes(x = subset, y = snps / 1e6, color = field)) +
  geom_point() +
  geom_line(data = predictions, aes(x = x_values, y = predicted_snp_counts / 1e6, color = field), alpha = 0.2) +
  geom_point(data = results1, aes(x = x_for_1_percent, y = snp_count / 1e6, color = field), shape = 8, size = 3) +
  geom_vline(data = results1, aes(xintercept = x_for_5_percent, color = field), linetype = "dotted") +
  scale_x_continuous(limits = c(0, 850)) +
  scale_y_continuous(limits = c(0, 7), name = "Number of SNPs (millions)") +
  scale_color_manual(values = custom_colors) +
  theme_bw() +
  labs(title = "Rarefaction Curves with Power Model Fits and 1% Increase Points", x = "n samples")
modelplot
dev.off()
