library(dplyr)
library(viridis)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
# Load the .bim file
custom_palette <- colorRampPalette(c("#313695", "#4575B4", "#ABD9E9", "#FFFFBF", "#FEE090", "#F46D43", "#A50026"))(100)

t2d_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/T2D/t2d_c4_group.txt", header=TRUE)
BMI_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/BMI/bmi_c4_group.txt", header=TRUE)
LDLC_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/LDLC/ldlc_c4_group.txt", header=TRUE)
HDLC_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/HDLC/hdlc_c4_group.txt", header=TRUE)
Trig_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/trig/trig_c4_group.txt", header=TRUE)
RG_group <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/RG/rq_c4_group.txt", header=TRUE)
# Function to rename columns to avoid duplicates during merge
convert_to_zscore <- function(df,df_name = "default", pvalue_col = "Pvalue", pvalue_skat_col = "Pvalue_SKAT", beta_col = "BETA_Burden", se_col = "SE_Burden") {
  # Convert P-value columns to z-scores and add as columns in the dataframe
  df[[paste0(df_name, "_z_skato")]] <- qnorm(1 - df[[pvalue_col]] / 2)
  df[[paste0(df_name, "_z_skat")]] <- qnorm(1 - df[[pvalue_skat_col]] / 2)
  
  # Calculate z-scores for Burden and add as a column in the dataframe
  df[[paste0(df_name, "_z_burden")]] <- df[[beta_col]] / df[[se_col]]
  
  selected_columns <- c("Region", "Group", "max_MAF", 
                        paste0(df_name, "_z_skato"), 
                        paste0(df_name, "_z_skat"), 
                        paste0(df_name, "_z_burden"))
  
  return(df[selected_columns])
}

t2d_group <- convert_to_zscore(t2d_group,df_name = "t2d")
BMI_group <- convert_to_zscore(BMI_group, df_name = "bmi")
LDLC_group <- convert_to_zscore(LDLC_group, df_name = "ldlc")
#HDLC_group <- convert_to_zscore(HDLC_group, df_name = "hdlc")
#Trig_group <- convert_to_zscore(Trig_group, df_name = "trig")
RG_group <- convert_to_zscore(RG_group, df_name = "rg")



# Assuming your dataframes are named df1, df2, df3, df4, df5

# List of dataframes
dfs <- list(t2d_group, BMI_group, LDLC_group, RG_group)

# Merge all dataframes by Region, Group, and max_MAF
merged_df <- Reduce(function(x, y) merge(x, y, by = c("Region", "Group", "max_MAF")), dfs)

# Select only columns with "_z_skato" in their names
z_skato_columns <- grep("_z_skato", names(merged_df), value = TRUE)
# Add the Region, Group, and max_MAF columns to the selection
skato_columns <- c("Region", "Group", "max_MAF", z_skato_columns)

# Select only columns with "_z_skato" in their names
z_skat_columns <- grep("_z_skat$", names(merged_df), value = TRUE)
# Add the Region, Group, and max_MAF columns to the selection
skat_columns <- c("Region", "Group", "max_MAF", z_skat_columns)

# Select only columns with "_z_skato" in their names
z_burden_columns <- grep("_z_burden", names(merged_df), value = TRUE)
# Add the Region, Group, and max_MAF columns to the selection
burden_columns <- c("Region", "Group", "max_MAF", z_burden_columns)


# Create the final dataframe with the selected columns
z_skato_df <- merged_df[, skato_columns]
z_skat_df <- merged_df[, skat_columns]
z_burden_df <- merged_df[, burden_columns]
# Function to generate a heatmap with Inf handling and annotations
generate_heatmap_with_inf <- function(df, z_score_columns, group_col = "Group", maf_col = "max_MAF", annotation_digits = 2, plot_title = "Heatmap") {
  df <- na.omit(df)
  
  max_finite_value <- max(unlist(df[z_score_columns])[is.finite(unlist(df[z_score_columns]))], na.rm = TRUE)
  placeholder_value <- max_finite_value * 1.5
  
  df_for_heatmap <- df
  df_for_heatmap[z_score_columns] <- lapply(df_for_heatmap[z_score_columns], function(x) {
    x[is.infinite(x)] <- placeholder_value
    return(x)
  })
  
  z_score_matrix <- as.matrix(df_for_heatmap[, z_score_columns])
  combined_row_names <- paste0(df[[group_col]], " (MAF: ", round(df[[maf_col]], 5), ")")
  rownames(z_score_matrix) <- combined_row_names
  
  annotation_matrix <- apply(z_score_matrix, c(1, 2), function(x) {
    if (x == placeholder_value) return("Inf")
    return(formatC(x, format = "f", digits = annotation_digits))
  })
  
  # Simplify x-axis labels
  simple_labels <- gsub("_z_(skato|skat|burden)", "", z_score_columns)
  pheatmap(z_score_matrix, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = annotation_matrix,  
           fontsize_number = 10,
           color = colorRampPalette(c("lightyellow", "orange", "darkred"))(100),  
           breaks = seq(0, max_finite_value, length.out = 100),
           labels_row = combined_row_names,
           labels_col = gsub("_z_(skato|skat|burden)", "", z_score_columns),
           main = plot_title,
           number_color = "black",  
           angle_col = 45)
}


# Generate heatmaps with the updated function and simplified labels
generate_heatmap_with_inf(z_skat_df, z_score_columns = c("t2d_z_skat", "bmi_z_skat", "ldlc_z_skat", "rg_z_skat"), plot_title = "Z-Scores from SKAT Test")

generate_heatmap_with_inf(z_skato_df, z_score_columns = c("t2d_z_skato", "bmi_z_skato", "ldlc_z_skato", "rg_z_skato"), plot_title = "Z-Scores from SKAT-O Test")

generate_heatmap_with_inf(z_burden_df, z_score_columns = c("t2d_z_burden", "bmi_z_burden", "ldlc_z_burden", "rg_z_burden"), plot_title = "Z-Scores from Burden Test")





t2d <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/T2D/t2d_c4.txt", header=TRUE)
BMI <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/BMI/bmi_c4.txt", header=TRUE)
LDLC <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/LDLC/ldlc_c4.txt", header=TRUE)
HDLC <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/HDLC/hdlc_c4.txt", header=TRUE)
Trig <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/trig/trig_c4.txt", header=TRUE)
RG <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/SAIGE_analysis/output/RG/rg_c4.txt", header=TRUE)

# Extract Burden p-values
burden_t2d <- data.frame(Group = t2d$Group, T2D_BURDEN_p = t2d$Pvalue_Burden)
burden_ldlc <- data.frame(Group = LDLC$Group, LDLC_BURDEN_p = LDLC$Pvalue_Burden)
burden_rg <- data.frame(Group = RG$Group, RG_BURDEN_p = RG$Pvalue_Burden)

# Merge burden p-values into one data frame
burden_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE),
                  list(burden_t2d, burden_ldlc, burden_rg))

# Extract SKAT p-values
skat_t2d <- data.frame(Group = t2d$Group, T2D_SKAT_p = t2d$Pvalue_SKAT)
skat_bmi <- data.frame(Group = BMI$Group, BMI_SKAT_p = BMI$Pvalue_SKAT)
skat_ldlc <- data.frame(Group = LDLC$Group, LDLC_SKAT_p = LDLC$Pvalue_SKAT)
skat_rg <- data.frame(Group = RG$Group, RG_SKAT_p = RG$Pvalue_SKAT)

# Merge SKAT p-values into one data frame
skat_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE),
                  list(skat_t2d, skat_bmi, skat_ldlc, skat_rg))

# Extract SKAT-O p-values
skato_t2d <- data.frame(Group = t2d$Group, T2D_SKATO_p = t2d$Pvalue)
skato_bmi <- data.frame(Group = BMI$Group, BMI_SKATO_p = BMI$Pvalue)
skato_ldlc <- data.frame(Group = LDLC$Group, LDLC_SKATO_p = LDLC$Pvalue)
skato_rg <- data.frame(Group = RG$Group, RG_SKATO_p = RG$Pvalue)

# Merge SKAT-O p-values into one data frame
skat_o_df <- Reduce(function(x, y) merge(x, y, by = "Group", all = TRUE),
                    list(skato_t2d, skato_bmi, skato_ldlc, skato_rg))












t2d <- convert_to_zscore(t2d,df_name = "t2d")
#BMI <- convert_to_zscore(BMI, df_name = "bmi")
LDLC <- convert_to_zscore(LDLC, df_name = "ldlc")
#HDLC <- convert_to_zscore(HDLC, df_name = "hdlc")
#Trig <- convert_to_zscore(Trig, df_name = "trig")
RG <- convert_to_zscore(RG, df_name = "rg")


new_dfs <- list(t2d, LDLC, RG)

# Merge all dataframes by Region and Group
merged_new_df <- Reduce(function(x, y) merge(x, y, by = c("Region", "Group")), new_dfs)

# Select only columns with "_z_burden" in their names
z_burden <- grep("_z_burden", names(merged_new_df), value = TRUE)
# Add the Region and Group columns to the selection
burden_columns <- c("Region", "Group", z_burden)
z_burden_df <- merged_new_df[, burden_columns]
z_burden_df <- na.omit(z_burden_df)

z_score_columns <- z_burden

# Define max and min finite values
max_finite_value <- max(unlist(z_burden_df[z_score_columns])[is.finite(unlist(z_burden_df[z_score_columns]))], na.rm = TRUE)
min_finite_value <- min(unlist(z_burden_df[z_score_columns])[is.finite(unlist(z_burden_df[z_score_columns]))], na.rm = TRUE)
positive_placeholder_value <- max_finite_value * 1.5
negative_placeholder_value <- min_finite_value * 1.5

# Replace infinite values in the z_score columns
z_burden_df[z_score_columns] <- lapply(z_burden_df[z_score_columns], function(x) {
  x[is.infinite(x) & x > 0] <- positive_placeholder_value  # Replace positive Inf with a large positive placeholder
  x[is.infinite(x) & x < 0] <- negative_placeholder_value  # Replace negative Inf with a large negative placeholder
  return(x)
})

# Convert dataframe to matrix for the heatmap
z_score_matrix <- as.matrix(z_burden_df[, z_score_columns])
rownames(z_score_matrix) <- z_burden_df[["Group"]]
significant_threshold <- 0.05/4



# Prepare the annotation matrix with bold formatting for significant values
annotation_matrix <- apply(z_score_matrix, c(1, 2), function(x) {
  # Find the row (group) and column (p-value type) corresponding to the current z-score value
  coords <- which(z_score_matrix == x, arr.ind = TRUE)
  group_name <- rownames(z_score_matrix)[coords[1, 1]]
  column_name <- colnames(z_score_matrix)[coords[1, 2]]
  
  # Determine the correct p-value column based on the z-score column name
  if (column_name == "t2d_z_burden") {
    p_column <- "T2D_BURDEN_p"
  } else if (column_name == "ldlc_z_burden") {
    p_column <- "LDLC_BURDEN_p"
  } else if (column_name == "rg_z_burden") {
    p_column <- "RG_BURDEN_p"
  } else {
    p_column <- NA  # In case the column is not recognized, set to NA
  }

  if (!is.na(p_column) && !is.na(group_name)) {
    p_value <- burden_df[burden_df$Group == group_name, p_column]
    
    # Check if p_value is present and less than the significant threshold
    formatted_value <- formatC(x, format = "f", digits = 2)
    if (!is.na(p_value) && length(p_value) > 0 && p_value < significant_threshold) {
      return(paste0(formatted_value, "*"))  # Append an asterisk for significance
    } else {
      return(formatted_value)
    }
  } else {
    return(formatC(x, format = "f", digits = 2))  # Return non-asterisk format if p-value or group is NA
  }
})



# Simplify x-axis labels
simple_labels <- gsub("_z_burden", "", z_score_columns)

# Adjust the breaks to ensure they are centered around zero
max_abs_value <- max(abs(c(min_finite_value, max_finite_value)))
breaks <- seq(-max_abs_value, max_abs_value, length.out = 100)
palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
modified_labels_row <- sapply(rownames(z_score_matrix), function(name) {
  if (name == "cluster1") {
    return("cluster1 (green)")
  } else if (name == "cluster2") {
    return("cluster2 (pink)")
  } else if (name == "cluster3") {
    return("cluster3 (cyan)")
  } else if (name == "cluster4") {
    return("cluster4 (purple)")
  } else {
    return(name)
  }
})
# Generate the heatmap
pheatmap(z_score_matrix, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = annotation_matrix,  
         fontsize_number = 14,     # Size of the numbers inside the heatmap
         fontsize_col = 14,        # Increase font size for column labels (ldlc, rg, t2d)
         fontsize_row = 14,
         color = palette,
         breaks = breaks,
         labels_row = modified_labels_row,
         labels_col = simple_labels,
         main = "Z-Scores from Burden Test\n(* indicates statistically significant values)",
         number_color = "black",  
         angle_col = 45)

