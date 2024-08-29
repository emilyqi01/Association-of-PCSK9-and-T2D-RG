library(dplyr)
RG_linear <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/RG_cluster.csv", header=TRUE)
RG_linear <- RG_linear[!is.na(RG_linear$P), ]
RG_linear <- RG_linear[order(RG_linear$P), ]
adjust_a = 0.05 / 226
gwas = 5 * 10^-8
# Filter rows where the p-value is less than the adjusted threshold
significant_rows <- RG_linear[RG_linear$P < adjust_a, ]
# Print out the filtered rows
print(significant_rows)


t2d_linear <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/t2d_cluster.csv", header=TRUE)
t2d_linear <- t2d_linear[!is.na(t2d_linear$P), ]
t2d_linear <- t2d_linear[order(t2d_linear$P), ]
adjust_a = 0.05 / 226
# Filter rows where the p-value is less than the adjusted threshold
significant_rows <- t2d_linear[t2d_linear$P < adjust_a, ]
# Print out the filtered rows
print(significant_rows)

# Read the data from the CSV file
chol_linear <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/clustering_analysis/chol_cluster.csv", header=TRUE)
# Remove rows where the p-value is NA
chol_linear <- chol_linear[!is.na(chol_linear$P), ]
# Remove rows where the p-value is exactly 0
chol_linear <- chol_linear[chol_linear$P != 0, ]
# Sort the dataframe by the p-value column in ascending order
chol_linear <- chol_linear[order(chol_linear$P), ]
# Set the adjusted significance threshold
adjust_a = 0.05 / 226
# Filter rows where the p-value is less than the adjusted threshold
significant_rows <- chol_linear[chol_linear$P < gwas, ]
# Print out the filtered rows
print(significant_rows)
