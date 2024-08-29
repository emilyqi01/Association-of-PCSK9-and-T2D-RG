library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

# Load the data
bmi_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/BMI_cluster.csv", header = TRUE)
bmi_table$z_bmi <- bmi_table$BETA / bmi_table$SE
bmi_table <- bmi_table %>% dplyr::select("SNP", "z_bmi")
ldlc_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/ldlc_cluster.csv", header = TRUE)
ldlc_table$z_ldlc <- ldlc_table$BETA / ldlc_table$SE
ldlc_table <- ldlc_table %>% dplyr::select("SNP", "z_ldlc")
hdlc_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/hdlc_cluster.csv", header = TRUE)
hdlc_table$z_hdlc <- hdlc_table$BETA / hdlc_table$SE
hdlc_table <- hdlc_table %>% dplyr::select("SNP", "z_hdlc")
RG_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/RG_cluster.csv", header = TRUE)
RG_table$z_rg <- RG_table$BETA / RG_table$SE
RG_table <- RG_table %>% dplyr::select("SNP", "z_rg")
t2d_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/t2d_cluster.csv", header = TRUE)
t2d_table$z_t2d <- t2d_table$BETA / t2d_table$SE
t2d_table <- t2d_table %>% dplyr::select("SNP", "z_t2d")
trig_table <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/triglycerides_cluster.csv", header = TRUE)
trig_table$z_trig <- trig_table$BETA / trig_table$SE
trig_table <- trig_table %>% dplyr::select("SNP", "z_trig")

# List of dataframes
dfs <- list(bmi_table, ldlc_table, hdlc_table, RG_table, t2d_table, trig_table)

# Use Reduce to merge all dataframes
merged_df <- Reduce(function(x, y) merge(x, y, by = "SNP"), dfs)

merged_df <- merged_df %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# Set the output name
out <- "output"
merged_df = merged_df  %>% dplyr::select("SNP", "z_t2d", "z_ldlc", "z_rg")
df_no_labels <- merged_df
df_no_labels_no_na <- df_no_labels[, 2:ncol(df_no_labels)]

## TRUNCATION
# Specify the phenotype (column) we want to truncate
target_phenotype <- "z_ldlc"

if (target_phenotype %in% names(df_no_labels_no_na)) {
  
  j <- target_phenotype
  
  # Replace Inf values with the maximum finite value
  if (Inf %in% df_no_labels_no_na[, j]) {
    ind <- which(df_no_labels_no_na[, j] == Inf)
    df_no_labels_no_na[ind, j] <- max(df_no_labels_no_na[-ind, j], na.rm = TRUE)
  }
  
  # Calculate the truncation threshold
  p_threshold <- 0.05 / 226
  
  # Convert the Bonferroni-corrected p-value to a z-score
  sd2 <- qnorm(p_threshold / 2, lower.tail = FALSE)  # two-tailed test
  # gwas significance
  #p_value <- 5*10^-8
  #sd2 <- abs(qnorm(p_value / 2))
  # Truncate values
  for (i in 1:nrow(df_no_labels_no_na)) {
    if (df_no_labels_no_na[i, j] > sd2) {
      df_no_labels_no_na[i, j] <- sd2
    } else if (df_no_labels_no_na[i, j] < -sd2) {
      df_no_labels_no_na[i, j] <- -sd2
    }
  }
} else {
  warning("The specified phenotype does not exist in the dataframe.")
}

# Reput in df_no_labels the imputation
df_no_labels[, 2:ncol(df_no_labels)] <- df_no_labels_no_na

# Set the number of clusters
k <- 4

# Convert the data frame to a matrix
df_matrix <- as.matrix(df_no_labels_no_na, rownames.force = TRUE)
second_column <- df_matrix[, 2]
# Determine the color breaks based on the data range
range_vals <- range(second_column, na.rm = TRUE)
color_breaks <- seq(range_vals[1], range_vals[2], length.out = 100)

# Use a more subtle color palette
palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

# Perform hierarchical clustering with pheatmap
cluster <- pheatmap(df_matrix, 
                    clustering_method = "ward.D2", 
                    breaks = color_breaks, 
                    color = palette, 
                    show_rownames = FALSE, 
                    width = 6, 
                    height = 12, 
                    cutree_rows = k)

# Extract cluster assignments
cluster_assignments <- cutree(cluster$tree_row, k)

# Add cluster assignments to the data frame
df_no_labels$Cluster <- factor(cluster_assignments)

# Create a data frame for row annotations
row_annotations <- data.frame(Cluster = factor(cluster_assignments))
rownames(row_annotations) <- rownames(df_matrix)
# Plot the heatmap with cluster annotations
cluster = pheatmap(df_matrix, 
         clustering_method = "ward.D2", 
         breaks = color_breaks, 
         color = palette, 
         show_rownames = FALSE, 
         width = 6, 
         height = 12, 
         cutree_rows = k, 
         annotation_row = row_annotations,
         filename = paste0("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/clustering_result/", out, "_hierarchical_clustering_new", k, "clusters_dendrogram.png"))



threshold <- max(cluster$tree_row$height) / 2 #typically half of the distance here
df_no_labels$cluster <- cutree(cluster$tree_row, h = threshold)
k <- length(unique(df_no_labels$cluster))


write.csv(df_no_labels, paste0("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/clustering_result/",out,"_hierarchical_clustering_new",k,"clusters_table.csv"), row.names = F)

#boxplot of phenos per cluster
for (cluster in unique(df_no_labels$cluster)) {
  
  temp_df <- df_no_labels_no_na[which(df_no_labels$cluster == cluster),]
  
  reformatted_df_no_labels_no_na <- data.frame()
  
  for (col in names(temp_df)) {
    
    temp_reformatted_df <- data.frame(z = temp_df[,col], 
                                      pheno = rep(col, nrow(temp_df)), 
                                      mean_value = rep(mean(temp_df[,col]),nrow(temp_df)))
    reformatted_df_no_labels_no_na <- rbind(reformatted_df_no_labels_no_na, temp_reformatted_df)
    
  }
  
  ggplot(reformatted_df_no_labels_no_na, aes(x = reorder(pheno, mean_value, mean), y = z, fill = mean_value)) +
    geom_boxplot() + 
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
    theme_minimal() + 
    labs(title = paste("Boxplot of cluster",cluster),y = "Z-score",x = "Phenotype") +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste0("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/clustering_result/",out,"_hierarchical_clustering_new", k, "cluster_",cluster,"_boxplot.png"), width = 10, height = 10)
  
}

# Set weighting factors for the specific phenotypes
k = 4
weight_ldlc <- 1.2
weight_t2d <- 1.2  
weight_rg <- 1.2  

# Apply the weights to the corresponding columns in the matrix
df_matrix_weighted <- df_matrix
df_matrix_weighted[, "z_ldlc"] <- df_matrix_weighted[, "z_ldlc"] * weight_ldlc
df_matrix_weighted[, "z_t2d"] <- df_matrix_weighted[, "z_t2d"] * weight_t2d
df_matrix_weighted[, "z_rg"] <- df_matrix_weighted[, "z_rg"] * weight_rg

# Determine the color breaks based on the range of the second column
second_column <- df_matrix_weighted[, "z_ldlc"]  # Assuming z_ldlc is the second column
range_vals <- range(second_column, na.rm = TRUE)
color_breaks <- seq(range_vals[1], range_vals[2], length.out = 100)

# Perform hierarchical clustering with pheatmap using the weighted matrix
cluster <- pheatmap(df_matrix_weighted, 
                    clustering_method = "ward.D2", 
                    breaks = color_breaks, 
                    color = palette, 
                    show_rownames = FALSE, 
                    width = 6, 
                    height = 12, 
                    cutree_rows = k,
                    filename = paste0("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/clustering_result/", out, "weighted_", k, "clusters_dendrogram.png"))


threshold <- max(cluster$tree_row$height) / 2 #typically half of the distance here
df_no_labels$cluster <- cutree(cluster$tree_row, k = k)
k <- length(unique(df_no_labels$cluster))


write.csv(df_no_labels, paste0("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/clustering_result/",out,"weighted",k,"clusters_table.csv"), row.names = F)
