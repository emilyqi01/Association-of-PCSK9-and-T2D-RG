# Load necessary libraries
library(dplyr)  # For data manipulation
library(tidyr)  # For data tidying functions
library(ggplot2)  # For data visualization

# Load input files from specified file paths

# Read the original genetic data file
pcsk9_extraction <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/pcsk9_extraction.bim", header = FALSE)
colnames(pcsk9_extraction) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")  # Assign column names

# Read the allele counts file
frqx <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/allele_counts.frq.counts", header = TRUE)

# Read the variant frequency file
frq_MAF <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/variant_frequencies.frq", header = TRUE)

# Read the Hardy-Weinberg equilibrium results file
hardy <- read.table("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/hwe_results.hwe", header = TRUE)

# Split the 'GENO' column into separate columns for observed homozygous, heterozygous, and homozygous alternate
hardy <- hardy %>%
  separate(GENO, into = c("O_HOM1", "O_HET", "O_HOM2"), sep = "/")

# Convert the new columns to numeric type
hardy <- hardy %>%
  mutate(across(starts_with("O_"), as.numeric))

# Merge data with allele frequencies
extraction <- left_join(pcsk9_extraction, frqx, by = c("CHR", "SNP", "A1", "A2"))

# Merge data with MAF (Minor Allele Frequency)
extraction_636 <- left_join(extraction, frq_MAF, by = c("CHR", "SNP", "A1", "A2"))

# Merge data with Hardy-Weinberg equilibrium results
extraction_hardy <- left_join(extraction_636, hardy, by = c("CHR", "SNP", "A1", "A2"))

# Add columns for the count of homozygous and heterozygous individuals
extraction_hardy <- extraction_hardy %>%
  mutate(`C(HOM A1)` = O_HOM1,
         `C(HET)` = O_HET)

# Identify and filter out singletons and doubletons (rare variants)
single_doubletons <- extraction_hardy %>%
  filter((`C(HOM A1)` + `C(HET)`) <= 2)

# Filter out singletons and doubletons from the main dataset
extraction_remove_single_doubletons <- extraction_hardy %>%
  filter(!(`C(HOM A1)` + `C(HET)`) <= 2)

# Remove indels (insertion-deletion polymorphisms) to focus on SNPs
extraction_delete_idels <- extraction_remove_single_doubletons %>%
  filter(A1 %in% c("A", "C", "T", "G")) %>%
  filter(A2 %in% c("A", "C", "T", "G"))

# Identify multiallelic sites
position_count <- extraction_delete_idels %>% count(BP) %>% filter(n > 1)  # Positions with more than one variant

# Filter out multiallelic variants to retain biallelic variants only
extraction_multiallelic <- extraction_delete_idels %>% filter(!BP %in% position_count$BP)

# Filter variants within the PCSK9 gene region based on base pair position (BP)
PCSK9_gene <- extraction_multiallelic %>% filter(BP >= 55039445 & BP <= 55064852)

# Save filtered variants to a file (list of variants to keep)
write.table(PCSK9_gene, file = "PCSK9_gene.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(PCSK9_gene[, c("SNP")], file = "filtered.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Load annotation data and preprocess for merging
annotation_table <- read.csv("/rds/general/user/yq523/projects/eph-prokopenko-lab-silk/live/shared/zhanna_yinglai/gnomAD_v4.1.0_ENSG00000169174_2024_06_18_12_29_18.csv", header = TRUE)
annotation_table$gnomAD.ID <- gsub("-", ":", annotation_table$gnomAD.ID)  # Replace hyphens with colons

# Find columns for merging
snp_col_index <- which(names(PCSK9_gene) == "SNP")
snp_col_index_ann <- which(names(annotation_table) == "gnomAD.ID")

# Expand rows for each SNP in PCSK9_gene for correct merging
expanded_df <- do.call(rbind, lapply(1:nrow(PCSK9_gene), function(i) {
  snps <- unlist(strsplit(PCSK9_gene[i, snp_col_index], ";"))
  df_replicated <- PCSK9_gene[i, , drop = FALSE]
  df_replicated <- df_replicated[rep(1, length(snps)), ]
  df_replicated[, snp_col_index] <- snps
  return(df_replicated)
}))

# Rename columns for merging
colnames(annotation_table)[snp_col_index_ann] <- "SNP"

# Merge with annotation table by SNP
merge_anno <- left_join(expanded_df, annotation_table, by = c("SNP"))

# Count variants by their functional annotation
grouped_counts <- merge_anno %>%
  group_by(VEP.Annotation) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Create a pie chart to visualize the distribution of variant annotations
ggplot(grouped_counts, aes(x = "", y = Count, fill = VEP.Annotation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 5) +
  labs(fill = "VEP Annotation", title = "Distribution of VEP Annotations") +
  theme(legend.position = "right")

# Define a set of lighter colors for the pie chart
unified_colors <- c("#ADD8E6", "#90EE90", "#F08080", "#FFFFE0", 
                    "#E5D8BD", "#C1E1C1", "#FFD700", "#FFB6C1", "#E0FFFF", "#DDA0DD", "#FFE4B5")

# Set a threshold for displaying numbers on the pie chart
threshold <- 15

# Create the pie chart with customized colors, size adjustments, and reduced white space
plot <- ggplot(grouped_counts, aes(x = "", y = Count, fill = VEP.Annotation)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 6) +  # Increase text size
  scale_fill_manual(values = unified_colors) +  # Apply custom colors
  labs(fill = "VEP Annotation", title = "Distribution of Variants After Annotations") +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),  # Increase title font size
        legend.text = element_text(size = 10),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +  # Reduce margins
  guides(fill = guide_legend(ncol = 1))  # Adjust legend layout

# Save the pie chart to a file with specified dimensions and DPI
ggsave("piechart_compact.png", plot = plot, width = 6, height = 4.5, dpi = 400, units = "in")

# Filter for missense variants from the merged annotation data
missense_variant <- merge_anno %>% filter(VEP.Annotation == "missense_variant")

# Save the list of missense variants
write.table(missense_variant[, c("SNP")], file = "missense_variant.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
