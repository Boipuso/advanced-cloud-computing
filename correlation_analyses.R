# script to conduct correlation analyses of the Indices time series values of single beach trees 
# for the advanced GEE course - BEECHDECLINE project 
# last updated: 2024-07-07
# author: Henning Riecken

# steps for correlation analysis 
# 1. Data pre-processing
# 2.1 Pearson correlation heatmap
# 2.2 Spearman multi-correlation
# 3.1 combined Scatter plots with stats
# 3.2 single Scatter plots with stats

##### 1. pre-process the data structure #####
setwd("C:/Users/henni/Documents/Uni/Eagle/Semester2/GEE2/Project/output/")

data <- read.csv("Indices_with_FID_cleaned.csv", sep = ';', header = T)

library(tidyverse)

# Reshape the data frame 
long_df <- data %>%
  pivot_longer(
    cols = -fid,  # Select all columns except 'id'
    names_to = c("Index", "Year"),  # Split column names into 'Index' and 'Year'
    names_pattern = "(.+)_(\\d+)"  # Regular expression to match 'Index' and 'Year'
  ) %>%
  pivot_wider(
    names_from = Index,  # Use 'Index' to create new columns
    values_from = value  # The values for the new columns come from the 'value' column
  )

# sort the table by year 
long_df <- long_df %>%
  mutate(Year = as.numeric(Year)) %>%  # Ensure Year is numeric
  arrange(Year)  # Sort by Year

# filter out the yars 2018 and 2020 since ground truth data is missing
filtered_df <- long_df %>%
  filter(!(Year %in% c(2018, 2020)))

# divide the "LL" column by 100 
filtered_df$LL <- filtered_df$LL / 100


# export filtered data frame 
# write.csv(filtered_df, "filtered_data.csv", row.names = FALSE)



#### 2.1 Pearson correlation - heatmap#### 
library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

# exclude 'fid' and 'Year' for the correlation analysis
cor_data <- filtered_df %>%
  select(-fid, -Year)

# Compute correlations between 'LL' and other columns including p-values and adjusted R-squared values
results <- cor_data %>%
  reframe(across(-LL, ~{
    model <- lm(LL ~ .)  # apply linear model
    summary_model <- summary(model)  # Get the summary of the model
    cor_test <- cor.test(LL, ., method = "pearson")
    cor_value <- cor_test$estimate # Extract correlation coefficient
    p_value <- cor_test$p.value # Extract p-value
    adjusted_r_squared <- summary_model$adj.r.squared  # Extract adjusted R-squared
    list(correlation = cor_value, p_value = p_value, adjusted_r_squared = adjusted_r_squared)
  }))


# correlation coefficiant i in a strange format
# Function to extract numeric value from the r values
extract_number <- function(item) {
  matches <- regmatches(item, regexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", item))
  if (length(matches) > 0) {
    return(as.numeric(matches))
  } else {
    return(NA)
  }
}

# Apply the function to each item in the list
numbers <- sapply(results, extract_number) # output is a matrix
print(numbers)

# Assign row names to the matrix
rownames(numbers) <- c("Correlation", "P-value", "adj.R-squared")
print(numbers)

# function to assure that the p-values are formatted correctly (e.g. 1e-12) if p-value < 0.01
format_values <- function(matrix, row, col) {
  if (row == 2 & "P-value" < 0.01)  {
    return(format(matrix[row, col], scientific = TRUE, digits = 3))  # Format p-values
  } else {
    return(format(matrix[row, col], nsmall = 2, digits = 1))  # Format other values with 2 digits after decimal
  }
}

# Apply the formatting to the entire matrix
formatted_matrix <- matrix(nrow = nrow(numbers), ncol = ncol(numbers), dimnames = dimnames(numbers))
for (i in 1:nrow(numbers)) {
  for (j in 1:ncol(numbers)) {
    formatted_matrix[i, j] <- format_values(numbers, i, j)
  }
}

# create a heatmap out of the results
pheatmap(numbers, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = formatted_matrix, main = "Heatmap of Indices")



#### 2.2 Spearman multi-correlation #####
# not really helping but just for interest
filtered_df$LL <- as.numeric(filtered_df$LL)
filtered_df$NDVI <- as.numeric(filtered_df$NDVI) 

spearman_cor <- cor(filtered_df, method = "spearman")

print(spearman_cor)



#### 3.1 combined Scatter plots with stats ####
library(gridExtra)

# Create the output folder if it doesn't exist
if (!dir.exists("scatter_plots")) {
  dir.create("scatter_plots")
}

# List of indices
indices <- c('DI', 'EVI', 'G', 'MSI', 'NDII', 'NDVI', 'NDWI', 'NIR', 'SWIR1', 'SWIR2')

# Function to create regression plots without title 
create_regression_plot <- function(data, response, predictor) {
  formula <- as.formula(paste(response, "~", predictor))
  model <- lm(formula, data = data)
  summary_model <- summary(model)
  cor_value <- cor(data[[response]], data[[predictor]])  # Extract correlation coefficient
  p_value <- summary_model$coefficients[2, 4]
  adj_r_squared <- summary_model$adj.r.squared  # Extract adjusted R-squared
  
  # Format p-value
  formatted_p_value <- ifelse(p_value < 0.001, format.pval(p_value, digits = 3, eps = 0.001), round(p_value, 3))
  
  # position the geom_label with the stats depending on the correlation value, for better visuality
  vjust_value <- ifelse(cor_value > 0, 4.0, 1.2)
  
  plot <- ggplot(data, aes_string(x = predictor, y = response)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, col = "blue") +
    geom_label(aes(x = Inf, y = Inf, label = paste("r =", round(cor_value, 2), "\np =", formatted_p_value, "\nadj.R² =", round(adj_r_squared, 2))), 
               hjust = 1.0, vjust = vjust_value, size = 3, color = "black", fill = "lightgrey") +
    theme_minimal()
  
  return(plot)
}

# Create regression plots for each index
plots <- lapply(indices, function(index) create_regression_plot(filtered_df, "LL", index))

# Combine all plots with a single title 
combined_plots <- grid.arrange(grobs = plots, ncol = 4, top = "Correlation of true Leaf Loss (LL) with predictor Indices", gp = gpar(fontface = "bold", cex = 3)) # gp doesnt work for some reason

# safe the plot
ggsave("scatter_plots/combined_plots.png", combined_plots, width = 40, height = 20, units = "cm")


#### 3.2 single Scatter plots with stats ####
# Function to create individual regression plots with title 
create_regression_plot_with_title <- function(data, response, predictor) {
  formula <- as.formula(paste(response, "~", predictor))
  model <- lm(formula, data = data)
  summary_model <- summary(model)
  cor_value <- cor(data[[response]], data[[predictor]])  # Extract correlation coefficient
  p_value <- summary_model$coefficients[2, 4] # Extract p-value
  adj_r_squared <- summary_model$adj.r.squared  # Extract adjusted R-squared
  
  # Format p-value
  formatted_p_value <- ifelse(p_value < 0.001, format.pval(p_value, digits = 3, eps = 0.001), round(p_value, 3))
  
  # position the geom_label with the stats depending on the correlation value, for better visuality
  vjust_value <- ifelse(cor_value > 0, 6.0, 1.2)
  
  plot <- ggplot(data, aes_string(x = predictor, y = response)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, col = "blue") +
    geom_label(aes(x = Inf, y = Inf, label = paste("r =", round(cor_value, 2), "\np =", formatted_p_value, "\nadj.R² =", round(adj_r_squared, 2))), 
               hjust = 1.1, vjust = vjust_value, size = 5, color = "black", fill = "lightgrey") +
    labs(title = paste('Correlation of', indices[i], 'with Leaf Loss'), x = indices[i], y = 'Leaf Loss in %') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  
  return(plot)
}

# safe the plot
for (i in 1:length(indices)) {
  plot <- create_regression_plot_with_title(filtered_df, "LL", indices[i])
  ggsave(filename = paste0("scatter_plots/", indices[i], ".png"), plot = plot, width = 8, height = 6)
}









