# Load libraries
library(ggplot2)  
library(gridExtra) 
library(MASS)  # For statistical functions like Box-Cox transformation

# Load lipid phenotype data from a specified file
lipid_phenotype = read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/T2D_depr_BMItransformed_RG_final_lipids.sample", header = TRUE)

# Function to create a histogram for a given column in the data
create_histogram <- function(data, column_name) {
  ggplot(data, aes_string(x = column_name)) + 
    geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +  
    labs(title = paste("Histogram of", column_name), x = column_name, y = "Frequency") +  # Set title and axis labels
    theme_minimal()  # Apply a minimal theme
}

# Create histograms for specific columns in the dataset
hdlc_histogram <- create_histogram(lipid_phenotype, "hdlc")
chol_histogram <- create_histogram(lipid_phenotype, "chol")
ldlc_histogram <- create_histogram(lipid_phenotype, "ldlc")
triglycerides_histogram <- create_histogram(lipid_phenotype, "triglycerides")

# Arrange the histograms in a 2x2 layout
grid.arrange(hdlc_histogram, chol_histogram, ldlc_histogram, triglycerides_histogram, nrow = 2, ncol = 2)

# Apply log transformation to specified columns to normalize the data
lipid_phenotype$hdlc_log <- log(lipid_phenotype$hdlc + 0.01)
lipid_phenotype$chol_log <- log(lipid_phenotype$chol + 0.01)
lipid_phenotype$ldlc_log <- log(lipid_phenotype$ldlc + 0.01)
lipid_phenotype$triglycerides_log <- log(lipid_phenotype$triglycerides + 0.01)

# Function to perform Box-Cox transformation, handling missing values (NA)
boxcox_transform <- function(column) {
  original_length <- length(column)  # Store the original length of the column
  na_indices <- is.na(column)  # Identify NA indices
  column_no_na <- column[!na_indices]  # Remove NA values
  
  # Perform Box-Cox transformation and find the optimal lambda
  bc <- boxcox(column_no_na ~ 1, lambda = seq(-5, 5, 0.1))
  lambda <- bc$x[which.max(bc$y)]
  
  # Apply the transformation using the optimal lambda
  transformed_no_na <- if(lambda == 0) {
    log(column_no_na)
  } else {
    (column_no_na^lambda - 1) / lambda
  }
  
  transformed <- rep(NA, original_length)  # Create a vector with the original length filled with NA
  transformed[!na_indices] <- transformed_no_na  # Replace non-NA values with transformed values
  
  return(transformed)  # Return the transformed column
}

# Apply Box-Cox transformation to specified columns
lipid_phenotype$hdlc_bc <- boxcox_transform(lipid_phenotype$hdlc)
lipid_phenotype$chol_bc <- boxcox_transform(lipid_phenotype$chol)
lipid_phenotype$ldlc_bc <- boxcox_transform(lipid_phenotype$ldlc)
lipid_phenotype$triglycerides_bc <- boxcox_transform(lipid_phenotype$triglycerides)

# Create histograms for the Box-Cox transformed data
hdlc_bc <- create_histogram(lipid_phenotype, "hdlc_bc")
chol_bc <- create_histogram(lipid_phenotype, "chol_bc")
ldlc_bc <- create_histogram(lipid_phenotype, "ldlc_bc")
triglycerides_bc <- create_histogram(lipid_phenotype, "triglycerides_bc")

# Arrange the histograms for the transformed data in a 2x2 layout
grid.arrange(hdlc_bc, chol_bc, ldlc_bc, triglycerides_bc, nrow = 2, ncol = 2)

# Rename 'eid' to 'FID' and 'iid' to 'IID' to standardize column names
names(lipid_phenotype)[names(lipid_phenotype) == 'eid'] <- 'FID'
names(lipid_phenotype)[names(lipid_phenotype) == 'iid'] <- 'IID'

# Remove any double quotes from all column names
names(lipid_phenotype) <- gsub('"', '', names(lipid_phenotype))

# Rearrange columns to put FID and IID at the front of the dataframe
cols_to_move <- c('FID', 'IID')
other_cols <- setdiff(names(lipid_phenotype), cols_to_move)
lipid_phenotype <- lipid_phenotype[, c(cols_to_move, other_cols)]

# Save the modified lipid phenotype data to a file
write.table(lipid_phenotype, file = "/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/lipid_analysis/lipid_phenotype.sample", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Check the unique values in the T2D column to understand its distribution
unique(lipid_phenotype$T2D)

# Convert 1 to 0 and 2 to 1 in the T2D column to modify its coding scheme
lipid_phenotype$T2D <- ifelse(lipid_phenotype$T2D == 1, 0, ifelse(lipid_phenotype$T2D == 2, 1, lipid_phenotype$T2D))

# Verify the conversion of T2D values
unique(lipid_phenotype$T2D)

# Save the modified data with converted T2D values
write.table(lipid_phenotype, "lipid_analysis/lipid_phenotype_converted.sample", row.names = FALSE, quote = FALSE, sep = "\t")

# Reload the modified data
lipid_phenotype = read.table(file = "/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/lipid_analysis/lipid_phenotype.sample", header = TRUE)


# Plotting original histograms without titles for different columns
p1 <- ggplot(lipid_phenotype, aes(x = RG)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(x = "RG", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p2 <- ggplot(lipid_phenotype, aes(x = BMI)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "lightgreen", color = "black", alpha = 0.7) +
  labs(x = "BMI", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(),  
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p3 <- ggplot(lipid_phenotype, aes(x = ldlc)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, fill = "lightcoral", color = "black", alpha = 0.7) +
  labs(x = "LDLC", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p4 <- ggplot(lipid_phenotype, aes(x = hdlc)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, fill = "lightpink", color = "black", alpha = 0.7) +
  labs(x = "HDLC", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p5 <- ggplot(lipid_phenotype, aes(x = chol)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, fill = "lightyellow", color = "black", alpha = 0.7) +
  labs(x = "Total Cholesterol", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

# Plotting normalized histograms without titles for different columns
p1_r <- ggplot(lipid_phenotype, aes(x = RG_ln)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(x = "Normalized RG", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p2_r <- ggplot(lipid_phenotype, aes(x = BMI_transformed)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, fill = "lightgreen", color = "black", alpha = 0.7) +
  labs(x = "Normalized BMI", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p3_r <- ggplot(lipid_phenotype, aes(x = ldlc_bc)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, fill = "lightcoral", color = "black", alpha = 0.7) +
  labs(x = "Normalized LDLC", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p4_r <- ggplot(lipid_phenotype, aes(x = hdlc_bc)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, fill = "lightpink", color = "black", alpha = 0.7) +
  labs(x = "Normalized HDLC", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

p5_r <- ggplot(lipid_phenotype, aes(x = chol_bc)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.3, fill = "lightyellow", color = "black", alpha = 0.7) +
  labs(x = "Normalized Total Cholesterol", y = "Density") +
  theme_minimal() +
  theme(plot.title = element_blank(), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5))

# Arrange all the plots in a grid: original histograms in the first row and normalized histograms in the second row
grid.arrange(p1, p2, p3, p5, p1_r, p2_r, p3_r, p5_r, ncol = 4)
