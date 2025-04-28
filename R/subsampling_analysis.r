#install.packages("pak")
#pak::pkg_install("Icens")
#pak::pkg_install("emdelponte/r4pde")
library(pak)
library(Icens)
library(r4pde)
library(knitr)
library(jtools)
library(minpack.lm)
library(ggplot2)
library(tidyverse)
setwd("/work_beegfs/suaph296/Zymoproj/merged")

# Load the data
df <- read.table("rarefaction_mac1_all.txt", header = TRUE)

# Plot the data
svg(file = "Rarefaction_all_mac1_2021.svg")
ggplot(df, aes(x = subset, y = snp_count / 1e6, color = country, shape = filter)) +
  geom_point() +
  scale_y_continuous(name = "Number of SNPs (millions)") +
  theme_bw()
dev.off()

# Separate the data into strict and relaxed data frames
strict_data <- df %>% filter(filter == "Strict")
relaxed_data <- df %>% filter(filter == "Relaxed")

# Prepare the data for the test
prepare_data <- function(data) {
  data %>%
    mutate(snp_count_millions = snp_count / 1e6) %>%
    select(subset, snp_count_millions, country) %>%
    rename(x = subset, Y = snp_count_millions)
}

grad_data_strict <- prepare_data(strict_data)
grad_data_relaxed <- prepare_data(relaxed_data)

# Logistic curve model function
logistic_growth <- function(x, K, r, x0) {
  K / (1 + exp(-r * (x - x0)))
}

# Extract statistics function
extract_stats <- function(model, data) {
  if (inherits(model, "lm")) {
    summary_model <- summary(model)
    coefficients <- summary_model$coefficients
    r_squared <- summary_model$r.squared
    
    data.frame(
      a_intercept = coefficients[1, 1],
      se_a = coefficients[1, 2],
      b_x = coefficients[2, 1],
      se_b = coefficients[2, 2],
      a_back_intercept = exp(coefficients[1, 1]),
      R2 = r_squared
    )
  } else if (inherits(model, "nls")) {
    coefficients <- summary(model)$coefficients
    r_squared <- 1 - sum(residuals(model)^2) / sum((data$Y - mean(data$Y))^2)
    
    data.frame(
      a_intercept = coefficients[1, 1],
      se_a = coefficients[1, 2],
      b_x = coefficients[2, 1],
      se_b = coefficients[2, 2],
      a_back_intercept = exp(coefficients[1, 1]),
      R2 = r_squared
    )
  }
}

# Function to fit models and extract stats for each country
fit_models_per_country <- function(data) {
  reg_exp <- lm(log(Y) ~ x, data = data)
  reg_p <- lm(log(Y) ~ log(x), data = data)
  reg_pm <- lm(log(Y) ~ log(x + 0.4), data = data)
  
  # Adjusted starting values and bounds for the logistic model
  start_vals <- list(K = max(data$Y), r = 0.1, x0 = mean(data$x))
  lower_bounds <- c(K = 0, r = 0, x0 = min(data$x))
  upper_bounds <- c(K = max(data$Y) * 10, r = 10, x0 = max(data$x))
  
  reg_log <- tryCatch(
    nls(Y ~ logistic_growth(x, K, r, x0), data = data,
        start = start_vals,
        algorithm = "port",
        lower = lower_bounds,
        upper = upper_bounds),
    error = function(e) NULL
  )
  
  reg_lin <- lm(Y ~ x, data = data)
  
  results <- bind_rows(
    cbind(Model = "Exponential", extract_stats(reg_exp, data)),
    cbind(Model = "Power", extract_stats(reg_p, data)),
    cbind(Model = "Modified_Power", extract_stats(reg_pm, data)),
    if (!is.null(reg_log)) cbind(Model = "Logistic", extract_stats(reg_log, data)) else NULL,
    cbind(Model = "Linear", extract_stats(reg_lin, data))
  )
  
  return(results)
}

# Apply the fitting process per country for relaxed data
results_per_country_relaxed <- grad_data_relaxed %>%
  group_by(country) %>%
  do(fit_models_per_country(.))

# Print and save the results
print(knitr::kable(results_per_country_relaxed, caption = "Model Fit Results Per Country (Relaxed)"))
write.table(results_per_country_relaxed, file = "results_table_per_country_relaxed.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Apply the fitting process per country for strict data
results_per_country_strict <- grad_data_strict %>%
  group_by(country) %>%
  do(fit_models_per_country(.))

# Print and save the results
print(knitr::kable(results_per_country_strict, caption = "Model Fit Results Per Country (Strict)"))
write.table(results_per_country_strict, file = "results_table_per_country_strict.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#for all of them the best fitted is power!!!
#so assuming power as the growth behaviour:
fit_power_model <- function(data) {
  reg_p <- lm(log(Y) ~ log(x), data = data)
  return(reg_p)
}
#estimate the flattening of the curve
estimate_flattening_point <- function(model, data) {
  # Predict SNP counts using the power model
  predict_snp <- function(x) {
    exp(predict(model, newdata = data.frame(x = x)))
  }
  
  # Calculate the derivative of the power model
  derivative <- function(x) {
    coef(model)[2] * exp(coef(model)[1]) * x^(coef(model)[2] - 1)
  }
  
  # Extrapolate beyond the observed data
  max_x <- max(data$x)
  x_values <- seq(max_x, max_x * 10, length.out = 1000)  
  derivatives <- sapply(x_values, derivative)
  
  # Find the point where the derivative is sufficiently small
  flattening_threshold <- 1e-6  
  flattening_point <- x_values[which.min(abs(derivatives - flattening_threshold))]
  
  # Predict the number of SNPs at the flattening point
  snps_at_flattening <- predict_snp(flattening_point)
  
  return(list(flattening_point = flattening_point, snps_at_flattening = snps_at_flattening))
}
#apply to both strict and relaxed data
# Apply the fitting process per country for strict data
results_per_country_strict <- grad_data_strict %>%
  group_by(country) %>%
  do({
    model <- fit_power_model(.)
    flattening_info <- estimate_flattening_point(model, .)
    data.frame(
      country = unique(.$country),
      flattening_point = flattening_info$flattening_point,
      snps_at_flattening = flattening_info$snps_at_flattening
    )
  })

# Print and save the results
print(knitr::kable(results_per_country_strict, caption = "Flattening Points and SNP Counts Per Country (Strict)"))
write.table(results_per_country_strict, file = "flattening_points_per_country_strict.txt", sep = "\t", row.names = FALSE, quote = FALSE)
#apply to both strict and relaxed data
# Apply the fitting process per country for relaxed data
results_per_country_relaxed <- grad_data_relaxed %>%
  group_by(country) %>%
  do({
    model <- fit_power_model(.)
    flattening_info <- estimate_flattening_point(model, .)
    data.frame(
      country = unique(.$country),
      flattening_point = flattening_info$flattening_point,
      snps_at_flattening = flattening_info$snps_at_flattening
    )
  })

# Print and save the results
print(knitr::kable(results_per_country_relaxed, caption = "Flattening Points and SNP Counts Per Country (Strict)"))
write.table(results_per_country_strict, file = "flattening_points_per_country_strict.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#both projections together escaled up 10 times
combined_results <- bind_rows(
  results_per_country_strict %>% mutate(filter = "Strict"),
  results_per_country_relaxed %>% mutate(filter = "Relaxed")
)
combined_results <- combined_results %>%
  arrange(country)
combined_results
#plot the results

# Fit the power model and predict SNP counts for each country
fit_and_predict_power_model <- function(data) {
  model <- lm(log(Y) ~ log(x), data = data)
  predict_snp <- function(x) {
    exp(predict(model, newdata = data.frame(x = x)))
  }
  return(list(model = model, predict_snp = predict_snp))
}

# Generate predictions for the strict dataset
predictions_strict <- grad_data_strict %>%
  group_by(country) %>%
  do({
    fit <- fit_and_predict_power_model(.)
    data.frame(
      x = seq(min(.$x), max(.$x) * 10, length.out = 1000),
      y = fit$predict_snp(seq(min(.$x), max(.$x) * 10, length.out = 1000)),
      country = unique(.$country),
      filter = "Strict"
    )
  })

# Generate predictions for the relaxed dataset
predictions_relaxed <- grad_data_relaxed %>%
  group_by(country) %>%
  do({
    fit <- fit_and_predict_power_model(.)
    data.frame(
      x = seq(min(.$x), max(.$x) * 10, length.out = 1000),
      y = fit$predict_snp(seq(min(.$x), max(.$x) * 10, length.out = 1000)),
      country = unique(.$country),
      filter = "Relaxed"
    )
  })

# Combine predictions from both datasets
predictions <- bind_rows(predictions_strict, predictions_relaxed)
predictions
#because the curves don't really flatten near the amount of samples we could realistically obtain in a survey, we assume that any less than an increase of 1% in the observed snp count is not meaningful for a survey, so we follow up the analysis with this:
# Function to fit the power model and predict SNP counts
fit_and_predict_power_model <- function(data) {
  model <- lm(log(Y) ~ log(x), data = data)
  predict_snp <- function(x) {
    exp(predict(model, newdata = data.frame(x = x)))
  }
  return(list(model = model, predict_snp = predict_snp))
}

# For strict_data
max_strict_per_country <- strict_data %>%
  group_by(country) %>%
  summarize(max_value = max(snp_count, na.rm = TRUE)) %>%
  mutate(five_percent = max_value * 0.01)

# For relaxed_data
max_relaxed_per_country <- relaxed_data %>%
  group_by(country) %>%
  summarize(max_value = max(snp_count, na.rm = TRUE)) %>%
  mutate(five_percent = max_value * 0.01)
# Print the results
print(max_strict_per_country)
print(max_relaxed_per_country)
# Function to find the sample size at which the increase is less than 1% of the max SNP count
find_sample_size_for_1_percent_increase <- function(data, max_snp_count) {
  fit <- fit_and_predict_power_model(data)
  x_values <- seq(min(data$x), max(data$x) * 10, length.out = 1000)
  y_values <- fit$predict_snp(x_values)
  
  # Calculate the threshold for 1% increase
  threshold <- max_snp_count * 1.05
  
  # Debugging: Check the lengths of y_values and threshold
  print(paste("Length of x_values:", length(x_values)))
  print(paste("Length of y_values:", length(y_values)))
  print(paste("Threshold value:", threshold))
  
  # Ensure the threshold is a single numeric value
  if (length(threshold) != 1) {
    stop("Threshold is not a single numeric value.")
  }
  
  # Find the sample size where the predicted SNP count is less than the threshold
  sample_size <- x_values[which(y_values >= threshold)[1]]
  return(sample_size)
}

# Calculate the maximum SNP count for each country from the observed data
max_snp_counts_strict <- grad_data_strict %>%
  group_by(country) %>%
  summarize(max_snp_count = max(Y))

max_snp_counts_relaxed <- grad_data_relaxed %>%
  group_by(country) %>%
  summarize(max_snp_count = max(Y))

# Initialise an empty data frame to store the results
sample_sizes_strict <- data.frame()

# Find the sample size for each country in the strict dataset
for (country in unique(grad_data_strict$country)) {
  country_data <- grad_data_strict %>% filter(country == !!country)
  max_snp_count <- max_snp_counts_strict %>% filter(country == !!country) %>% pull(max_snp_count)
  sample_size <- find_sample_size_for_1_percent_increase(country_data, max_snp_count)
  
  sample_sizes_strict <- rbind(sample_sizes_strict, data.frame(country = country, sample_size = sample_size, type = "Strict"))
}
print(sample_sizes_strict)
# Initialise an empty data frame to store the results
sample_sizes_relaxed <- data.frame()
# Find the sample size for each country in the relaxed dataset
for (country in unique(grad_data_relaxed$country)) {
  country_data <- grad_data_relaxed %>% filter(country == !!country)
  max_snp_count <- max_snp_counts_relaxed %>% filter(country == !!country) %>% pull(max_snp_count)
  sample_size <- find_sample_size_for_1_percent_increase(country_data, max_snp_count)
  
  sample_sizes_relaxed <- rbind(sample_sizes_relaxed, data.frame(country = country, sample_size = sample_size, type = "Relaxed"))
}
print(sample_sizes_relaxed)

# Combine the strict and relaxed data frames into a single data frame
predictions_at_1 <- rbind(
  sample_sizes_strict %>% mutate(threshold_value = c(1509969.3, 1253291.55, 601444.2)),
  sample_sizes_relaxed %>% mutate(threshold_value = c(2128122.15, 1743038.85, 842622.9))
)

# Rename the 'type' column to 'filter'
predictions_at_1 <- predictions_at_1 %>% rename(filter = type)

# Print the combined data frame
print(predictions_at_1)

# Plot the data points, fitted power model, and annotations
svg(file = "Rarefaction_all_mac1_2021_with_power_model.svg")
modelplot <- ggplot(df, aes(x = subset, y = snp_count / 1e6, color = country, shape = filter)) +
  geom_point() +
  geom_line(data = predictions, aes(x = x, y = y / 1e6, color = country, linetype = filter), alpha = 0.2) +
  geom_point(data = combined_results, aes(x = flattening_point, y = snps_at_flattening / 1e6, color = country, shape = filter), size = 3) +
  geom_point(data = predictions_at_1, aes(x = sample_size, y = threshold_value / 1e6), shape = 8, color = "grey", size = 3) +
  geom_vline(data = predictions_at_1, aes(xintercept = sample_size), linetype = "dotted", color = "grey") +
  theme_bw() +
  labs(title = "Rarefaction Curves with Power Model Fits and Flattening Points")
modelplot
dev.off()

#now with the reasonable numbers predicted for meaningful saturation of diversity (using just the predictions that have an increase of more than 1%)
svg(file = "Rarefaction_all_mac1_2021_with_power_model_1cutoff.svg")
modelplot <- ggplot(df, aes(x = subset, y = snp_count / 1e6, color = country, shape = filter)) +
  geom_point() +
  geom_line(data = predictions, aes(x = x, y = y / 1e6, color = country, linetype = filter), alpha = 0.2) +
  geom_point(data = combined_results, aes(x = flattening_point, y = snps_at_flattening / 1e6, color = country, shape = filter), size = 3) +
  geom_point(data = predictions_at_1, aes(x = sample_size, y = threshold_value / 1e6), shape = 8, color = "grey", size = 3) +
  geom_vline(data = predictions_at_1, aes(xintercept = sample_size), linetype = "dotted", color = "grey") +
  scale_x_continuous(limits = c(0, 200)) +
  scale_y_continuous(limits = c(0, 3), name = "Number of SNPs (millions)") +
  theme_bw() +
  labs(title = "Rarefaction Curves with Power Model Fits and Flattening Points")

modelplot
dev.off()
###just for relaxed
df_relaxed <- df %>% filter(filter == "Relaxed")
# Filter the predictions, combined_results, and predictions_at_1 data frames
predictions_relaxed <- predictions %>% filter(filter == "Relaxed")
combined_results_relaxed <- combined_results %>% filter(filter == "Relaxed")
predictions_at_1_relaxed <- predictions_at_1 %>% filter(filter == "Relaxed")

# Plot the data points and the predictions for the relaxed dataset
svg(file = "Rarefaction_relaxed_mac1_2021_with_power_model_1cutoff_relx.svg")
modelplot <- ggplot(df_relaxed, aes(x = subset, y = snp_count / 1e6, color = country)) +
  geom_point() +
  geom_line(data = predictions_relaxed, aes(x = x, y = y / 1e6, color = country), alpha = 0.2) +
  geom_point(data = combined_results_relaxed, aes(x = flattening_point, y = snps_at_flattening / 1e6, color = country), size = 3) +
  geom_point(data = predictions_at_1_relaxed, aes(x = sample_size, y = threshold_value / 1e6), shape = 8, color = "grey", size = 3) +
  geom_vline(data = predictions_at_1_relaxed, aes(xintercept = sample_size), linetype = "dotted", color = "grey") +
  scale_x_continuous(limits = c(0, 200)) +
  scale_y_continuous(limits = c(0, 3), name = "Number of SNPs (millions)") +
  theme_bw() 

modelplot
dev.off()
