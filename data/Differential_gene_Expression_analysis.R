setwd("C:/Users/USER/Desktop/D_BCSIR/Genecount")

# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(openxlsx)
library(DEGreport) # for additional quality control plots
library(ashr) # for improved log fold change shrinkage
library(RColorBrewer)

# Create a directory for results if it doesn't exist
if (!dir.exists("DESeq2_Results")) {
  dir.create("DESeq2_Results")
}

# Function to create formatted Excel worksheets
create_formatted_sheet <- function(wb, sheet_name, data_frame, significance_threshold = 0.05) {
  # Add data to worksheet
  writeData(wb, sheet = sheet_name, x = data_frame, rowNames = TRUE)
  
  # Get dimensions
  n_rows <- nrow(data_frame) + 1  # +1 for header
  n_cols <- ncol(data_frame) + 1  # +1 for row names
  
  # Create style for significant genes
  sig_style <- createStyle(fgFill = "#E6F3FF") # Light blue background
  
  # Identify significant rows (padj < threshold)
  if ("padj" %in% colnames(data_frame)) {
    sig_rows <- which(data_frame$padj < significance_threshold) + 1 # +1 for header row
    if (length(sig_rows) > 0) {
      addStyle(wb, sheet = sheet_name, style = sig_style, 
               rows = sig_rows, cols = 1:n_cols, gridExpand = TRUE)
    }
  }
  
  # Freeze top row and first column
  freezePane(wb, sheet = sheet_name, firstRow = TRUE, firstCol = TRUE)
  
  # Auto-adjust column widths
  setColWidths(wb, sheet = sheet_name, cols = 1:n_cols, widths = "auto")
}

###############Data Loading and Quality Control################################
# Read count data and metadata
count_data <- read.xlsx("rnaseq.xlsx", sheet = "gene", colNames = TRUE, rowNames = TRUE)
count_data <- as.matrix(count_data)
metadata <- read.xlsx("rnaseq.xlsx", sheet = "m", colNames = TRUE, rowNames = TRUE)

# Basic QC checks
cat("=== DATA QUALITY CONTROL ===\n")
cat("Count matrix dimensions:", dim(count_data), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")
cat("Samples in count data:", colnames(count_data), "\n")
cat("Samples in metadata:", rownames(metadata), "\n")

# Check for missing values
cat("Missing values in count data:", sum(is.na(count_data)), "\n")
cat("Zero counts in data:", sum(count_data == 0), "(", 
    round(sum(count_data == 0) / length(count_data) * 100, 2), "%)\n")

# Library size statistics
lib_sizes <- colSums(count_data)
cat("Library size summary:\n")
print(summary(lib_sizes))

# Gene detection statistics
genes_detected <- colSums(count_data > 0)
cat("Genes detected per sample summary:\n")
print(summary(genes_detected))

####################Create DESeq2 Object and Pre-processing###########################
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = metadata,
                              design = ~ GroupS)

# Pre-filtering: Remove genes with very low counts
# Keep genes with at least 10 counts in at least 6 samples (smallest group size)
#keep <- rowSums(counts(dds) >= 10) >= 6
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "/", nrow(count_data), 
    "(", round(nrow(dds)/nrow(count_data)*100, 1), "%)\n")

# Set reference level
dds$GroupS <- relevel(dds$GroupS, ref = "Control")

dds <- DESeq(dds)



############### CPM Normalization and Export ################################

# Function to calculate CPM
calculate_cpm <- function(count_matrix) {
  # Add pseudocount to avoid log(0) issues if needed for downstream analysis
  # For raw CPM, we don't add pseudocount
  lib_sizes <- colSums(count_matrix)
  cpm_matrix <- t(t(count_matrix) / lib_sizes) * 1e6
  return(cpm_matrix)
}

# Calculate CPM
cpm_data <- calculate_cpm(count_data)

# Log2 transform CPM (commonly used for visualization)
log2_cpm_data <- log2(cpm_data + 1)  # Adding 1 to avoid log2(0)

# Create a summary data frame
cpm_summary <- data.frame(
  Sample = colnames(count_data),
  Total_Reads = lib_sizes,
  Genes_Detected = genes_detected,
  Mean_CPM = colMeans(cpm_data),
  Median_CPM = apply(cpm_data, 2, median),
  SD_CPM = apply(cpm_data, 2, sd)
)

# Create an Excel workbook with multiple sheets
cpm_wb <- createWorkbook()

# Sheet 1: Raw CPM values
addWorksheet(cpm_wb, "Raw_CPM")
writeData(cpm_wb, sheet = "Raw_CPM", x = cpm_data, rowNames = TRUE)

# Sheet 2: Log2 transformed CPM
addWorksheet(cpm_wb, "Log2_CPM")
writeData(cpm_wb, sheet = "Log2_CPM", x = log2_cpm_data, rowNames = TRUE)

# Sheet 3: CPM summary statistics
addWorksheet(cpm_wb, "CPM_Summary")
writeData(cpm_wb, sheet = "CPM_Summary", x = cpm_summary)

# Sheet 4: Sample metadata (for reference)
addWorksheet(cpm_wb, "Metadata")
writeData(cpm_wb, sheet = "Metadata", x = metadata, rowNames = TRUE)

# Sheet 5: Original counts (for comparison)
addWorksheet(cpm_wb, "Original_Counts")
writeData(cpm_wb, sheet = "Original_Counts", x = count_data, rowNames = TRUE)

# Apply formatting to all sheets
for(sheet in c("Raw_CPM", "Log2_CPM", "Original_Counts")) {
  # Freeze top row and first column
  freezePane(cpm_wb, sheet = sheet, firstRow = TRUE, firstCol = TRUE)
  
  # Auto-adjust column widths
  n_cols <- ncol(count_data) + 1
  setColWidths(cpm_wb, sheet = sheet, cols = 1:n_cols, widths = "auto")
}

# Format summary sheet
freezePane(cpm_wb, sheet = "CPM_Summary", firstRow = TRUE)
setColWidths(cpm_wb, sheet = "CPM_Summary", cols = 1:7, widths = "auto")

# Add conditional formatting for CPM summary
mean_style <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
high_style <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")

conditionalFormatting(cpm_wb, "CPM_Summary", 
                      cols = 4, rows = 2:(nrow(cpm_summary) + 1),
                      rule = "< 10", style = mean_style)
conditionalFormatting(cpm_wb, "CPM_Summary", 
                      cols = 4, rows = 2:(nrow(cpm_summary) + 1),
                      rule = "> 100", style = high_style)

# Save the Excel file
cpm_filename <- "DESeq2_Results/CPM_Normalized_Data.xlsx"
saveWorkbook(cpm_wb, cpm_filename, overwrite = TRUE)

cat("\n=== CPM NORMALIZATION COMPLETE ===\n")
cat("CPM Excel file saved as:", cpm_filename, "\n")
cat("Sheets included:\n")
cat("  1. Raw_CPM: Counts Per Million values\n")
cat("  2. Log2_CPM: Log2 transformed CPM values\n")
cat("  3. CPM_Summary: Summary statistics per sample\n")
cat("  4. Metadata: Sample metadata\n")
cat("  5. Original_Counts: Raw count data\n")

# Create visualization of CPM distributions
pdf("DESeq2_Results/CPM_Quality_Plots.pdf", width = 12, height = 8)

# Boxplot of log2 CPM distributions
par(mfrow = c(1, 2))
boxplot(log2_cpm_data, 
        main = "Log2 CPM Distribution per Sample",
        xlab = "Samples", ylab = "Log2(CPM + 1)",
        las = 2, col = "lightblue")
abline(h = median(log2_cpm_data), col = "red", lty = 2)

# Density plot of CPM distributions
plot(density(log2_cpm_data[,1]), 
     main = "Density of Log2 CPM Values",
     xlab = "Log2(CPM + 1)", ylab = "Density",
     col = "black", lwd = 2, xlim = c(0, 20))

for(i in 2:ncol(log2_cpm_data)) {
  lines(density(log2_cpm_data[,i]), 
        col = rainbow(ncol(log2_cpm_data))[i])
}

# Add legend for first few samples
legend("topright", 
       legend = colnames(log2_cpm_data)[1:min(5, ncol(log2_cpm_data))],
       col = rainbow(min(5, ncol(log2_cpm_data))),
       lty = 1, lwd = 2, cex = 0.7)

# Correlation heatmap of samples based on CPM
cor_matrix <- cor(log2_cpm_data, method = "spearman")

pheatmap(cor_matrix,
         main = "Sample Correlation (Spearman, log2 CPM)",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 8)

# Library size vs CPM mean
plot(lib_sizes, cpm_summary$Mean_CPM,
     main = "Library Size vs Mean CPM",
     xlab = "Total Reads", ylab = "Mean CPM",
     pch = 19, col = "blue", cex = 1.5)
text(lib_sizes, cpm_summary$Mean_CPM, 
     labels = cpm_summary$Sample, pos = 3, cex = 0.7)

dev.off()

cat("CPM quality plots saved as: DESeq2_Results/CPM_Quality_Plots.pdf\n")

# Create a separate CSV file for easy import into other tools
write.csv(cpm_data, "DESeq2_Results/CPM_Normalized_Data.csv")
write.csv(log2_cpm_data, "DESeq2_Results/Log2_CPM_Data.csv")

cat("CPM data also saved as CSV files for easy import.\n")
cat("Files created:\n")
cat("  - DESeq2_Results/CPM_Normalized_Data.csv\n")
cat("  - DESeq2_Results/Log2_CPM_Data.csv\n")




#######PcoA###################################
#################### PCoA Analysis ###################################

# Load additional required libraries
library(ggplot2)
library(ggrepel)
library(ape)
library(viridis)
library(vegan)

# Variance stabilizing transformation for PCoA
vsd <- vst(dds, blind = FALSE)

# Calculate distance matrix (Bray-Curtis distance)
#sample_dists <- dist(t(assay(vsd)), method = "euclidean")
# Using Bray-Curtis distance (commonly used in ecology)
sample_dists_bray <- vegdist(t(assay(vsd)), method = "bray")
pcoa_result <- pcoa(sample_dists_bray)

# Or Jaccard distance for presence/absence
#sample_dists_jaccard <- vegdist(t(assay(vsd)), method = "jaccard")
#pcoa_result_jaccard <- pcoa(sample_dists_jaccard)


# Alternative distance measures you can try:
# sample_dists <- dist(t(assay(vsd)), method = "manhattan")
# sample_dists <- vegdist(t(assay(vsd)), method = "bray") # from vegan package

# Perform PCoA (Principal Coordinates Analysis)
#pcoa_result <- pcoa(sample_dists)

# Extract coordinates
pcoa_scores <- as.data.frame(pcoa_result$vectors)
colnames(pcoa_scores) <- paste0("PCo", 1:ncol(pcoa_scores))

# Add sample metadata
pcoa_scores$Sample <- rownames(pcoa_scores)
pcoa_scores$Group <- metadata[rownames(pcoa_scores), "GroupS"]
pcoa_scores$GroupFull <- metadata[rownames(pcoa_scores), "Group"]

# Calculate variance explained by each axis
variance_explained <- round(pcoa_result$values$Relative_eig[1:5] * 100, 2)

cat("Variance explained by PCoA axes:\n")
for(i in 1:5) {
  cat(paste0("PCo", i, ": ", variance_explained[i], "%\n"))
}

# Create color palette
group_colors <- c("Control" = "#1f77b4", 
                  "CDF" = "#ff7f0e", 
                  "DHF" = "#2ca02c", 
                  "DSS" = "#d62728")

# Create PCoA plot
R1 <- ggplot(pcoa_scores, aes(x = PCo1, y = PCo2, color = Group, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(name = "Group", 
                     values = group_colors,
                     labels = c("Control", "Classical Dengue", "Dengue Hemorrhagic", "Dengue Shock")) +
  labs(title = "PCoA Plot of RNA-seq Samples",
       x = paste0("PCo1 (", variance_explained[1], "%)"),
       y = paste0("PCo2 (", variance_explained[2], "%)"),
       caption = "Based on Bray distance of VST-transformed counts") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  ) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  scale_fill_manual(values = group_colors)

Fig1 <- R1 +  theme_bw()

Fig1


# Save the plot
ggsave("DESeq2_Results/PCoA_plot.png", pcoa_plot, width = 10, height = 8, dpi = 300)
ggsave("DESeq2_Results/PCoA_plot.pdf", pcoa_plot, width = 10, height = 8)

# Create additional PCoA plot without sample labels (cleaner version)
pcoa_plot_clean <- ggplot(pcoa_scores, aes(x = PCo1, y = PCo2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(name = "Group", 
                     values = group_colors,
                     labels = c("Control", "Classical Dengue", "Dengue Hemorrhagic", "Dengue Shock")) +
  labs(title = "PCoA Plot of RNA-seq Samples",
       x = paste0("PCo1 (", variance_explained[1], "%)"),
       y = paste0("PCo2 (", variance_explained[2], "%)")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    aspect.ratio = 1
  ) +
  stat_ellipse(level = 0.68, linetype = 2, size = 0.5)

print(pcoa_plot_clean)
ggsave("DESeq2_Results/PCoA_plot_clean.png", pcoa_plot_clean, width = 10, height = 8, dpi = 300)

# PCoA with different axes (PCo1 vs PCo3)
if(ncol(pcoa_scores) >= 3) {
  pcoa_plot_13 <- ggplot(pcoa_scores, aes(x = PCo1, y = PCo3, color = Group, label = Sample)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(name = "Group", 
                       values = group_colors,
                       labels = c("Control", "Classical Dengue", "Dengue Hemorrhagic", "Dengue Shock")) +
    labs(title = "PCoA Plot - PCo1 vs PCo3",
         x = paste0("PCo1 (", variance_explained[1], "%)"),
         y = paste0("PCo3 (", variance_explained[3], "%)")) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      aspect.ratio = 1
    )
  
  print(pcoa_plot_13)
  ggsave("DESeq2_Results/PCoA_plot_PC1_PC3.png", pcoa_plot_13, width = 10, height = 8, dpi = 300)
}

# Save PCoA coordinates and results
write.csv(pcoa_scores, "DESeq2_Results/PCoA_coordinates.csv")
cat("PCoA coordinates saved to: DESeq2_Results/PCoA_coordinates.csv\n")

# Create a summary of PCoA results
pcoa_summary <- data.frame(
  Axis = paste0("PCo", 1:5),
  VarianceExplained = variance_explained[1:5],
  CumulativeVariance = cumsum(variance_explained[1:5])
)

print(pcoa_summary)
write.csv(pcoa_summary, "DESeq2_Results/PCoA_variance_summary.csv")

cat("\n=== PCoA ANALYSIS COMPLETED ===\n")
cat("Plots saved in: DESeq2_Results/\n")







#################### Statistical Analysis for PCoA ####################

# 1. PERMANOVA (Permutational Multivariate Analysis of Variance)
cat("\n=== PERMANOVA TEST ===\n")
permanova_result <- adonis2(sample_dists_bray ~ Group, data = pcoa_scores, permutations = 999)
print(permanova_result)

# 2. Pairwise PERMANOVA between groups
cat("\n=== PAIRWISE PERMANOVA ===\n")
# Helper function for pairwise PERMANOVA
# Alternative simpler version if the above still has issues:
# Corrected pairwise PERMANOVA function
pairwise.adonis2 <- function(dist_mat, group_var, data, permutations = 999) {
  groups <- data[[group_var]]
  group_levels <- unique(groups)
  results <- list()
  
  for (i in 1:(length(group_levels)-1)) {
    for (j in (i+1):length(group_levels)) {
      group1 <- group_levels[i]
      group2 <- group_levels[j]
      
      # Get indices for these two groups
      idx <- which(groups %in% c(group1, group2))
      
      # Subset distance matrix
      dist_subset <- as.dist(as.matrix(dist_mat)[idx, idx])
      
      # Create data subset
      data_subset <- data[idx, , drop = FALSE]
      
      # Create a new grouping variable with only 2 levels
      data_subset$PairwiseGroup <- factor(data_subset[[group_var]])
      
      # Run PERMANOVA - FIXED FORMULA
      result <- adonis2(dist_subset ~ PairwiseGroup, data = data_subset, permutations = permutations)
      
      # Store result with meaningful name
      comp_name <- paste(group1, "vs", group2)
      results[[comp_name]] <- result
    }
  }
  return(results)
}


pairwise_permanova <- pairwise.adonis2(sample_dists_bray, "Group", pcoa_scores, permutations = 999)
print(pairwise_permanova)


# 3. Homogeneity of dispersion test (PERMDISP)
cat("\n=== HOMOGENEITY OF DISPERSIONS TEST ===\n")
dispersion_test <- betadisper(sample_dists_bray, pcoa_scores$Group)
permdisp_result <- permutest(dispersion_test, permutations = 999)
print(permdisp_result)

# 4. Calculate group centroids and distances
group_centroids <- betadisper(sample_dists_bray, pcoa_scores$Group)
centroid_distances <- as.matrix(dist(group_centroids$centroids))
cat("\n=== DISTANCES BETWEEN GROUP CENTROIDS ===\n")
print(centroid_distances)

# 5. Correlation between PCoA axes and original variables
cat("\n=== CORRELATION WITH PCoA AXES ===\n")
# If you have additional metadata variables, you can add them here
# For example, if you had Age or Sex in metadata:
# env_data <- metadata[, c("Age", "Sex"), drop = FALSE]
# envfit_result <- envfit(pcoa_result, env_data, permutations = 999)
# print(envfit_result)

# 6. Save statistical results to files
# Create a list to store all results
stat_results <- list(
  PERMANOVA = permanova_result,
  Pairwise_PERMANOVA = pairwise_permanova,
  Dispersion_test = permdisp_result,
  Centroid_distances = centroid_distances,
  PCoA_variance = pcoa_summary,
  Coordinates = pcoa_scores
)

# Save detailed results to RData file
saveRDS(stat_results, "DESeq2_Results/PCoA_statistical_results.rds")

# Save summary results to text file
sink("DESeq2_Results/PCoA_statistical_summary.txt")
cat("PCoA STATISTICAL ANALYSIS SUMMARY\n")
cat("=================================\n\n")

cat("1. PERMANOVA Results:\n")
cat("-------------------\n")
print(permanova_result)

cat("\n2. Variance Explained by PCoA Axes:\n")
cat("----------------------------------\n")
print(pcoa_summary)

cat("\n3. Homogeneity of Dispersions Test:\n")
cat("----------------------------------\n")
print(permdisp_result)

cat("\n4. Distances Between Group Centroids:\n")
cat("------------------------------------\n")
print(centroid_distances)
sink()

# 7. Create enhanced PCoA plot with statistical annotations
# Extract p-value from PERMANOVA
permanova_p <- permanova_result$`Pr(>F)`[1]

# Create enhanced plot with statistical info
Fig1_enhanced <- Fig1 + 
  labs(subtitle = paste0("PERMANOVA: F = ", round(permanova_result$F[1], 3), 
                         ", p = ", format.pval(permanova_p, digits = 3),
                         "\nDispersion test p = ", format.pval(permdisp_result$tab$`Pr(>F)`[1], digits = 3))) +
  theme(plot.subtitle = element_text(size = 10, hjust = 0.5))

print(Fig1_enhanced)
ggsave("DESeq2_Results/PCoA_plot_with_stats.png", Fig1_enhanced, width = 10, height = 8, dpi = 300)
Fig1E <- Fig1_enhanced
Fig1E

# 8. Create centroid plot
# Corrected centroid plot code
centroid_data <- as.data.frame(group_centroids$centroids)
centroid_data$Group <- rownames(centroid_data)

# Check the actual column names in centroid_data
cat("Column names in centroid_data:\n")
print(colnames(centroid_data))

# Rename the columns to match pcoa_scores if needed
# The centroids typically have columns like PC1, PC2, etc.
if(all(c("PC1", "PC2") %in% colnames(centroid_data))) {
  colnames(centroid_data)[colnames(centroid_data) %in% c("PC1", "PC2")] <- c("PCo1", "PCo2")
}

# Alternative approach: manually extract the correct coordinates
centroid_data <- as.data.frame(group_centroids$centroids[, 1:2])
colnames(centroid_data) <- c("PCo1", "PCo2")
centroid_data$Group <- rownames(centroid_data)

# Now create the plot
centroid_plot <- ggplot() +
  geom_point(data = pcoa_scores, aes(x = PCo1, y = PCo2, color = Group), alpha = 0.6, size = 2) +
  geom_point(data = centroid_data, aes(x = PCo1, y = PCo2, color = Group), size = 6, shape = 17) +
  geom_text_repel(data = centroid_data, aes(x = PCo1, y = PCo2, label = Group), 
                  size = 5, fontface = "bold") +
  scale_color_manual(values = group_colors) +
  labs(title = "PCoA Plot with Group Centroids",
       x = paste0("PCo1 (", variance_explained[1], "%)"),
       y = paste0("PCo2 (", variance_explained[2], "%)")) +
  theme_bw() +
  theme(legend.position = "none", aspect.ratio = 1)

print(centroid_plot)
ggsave("DESeq2_Results/PCoA_centroids_plot.png", centroid_plot, width = 10, height = 8, dpi = 300)
Fig2 <- centroid_plot
Fig2

# 9. Create results summary table for writing
results_summary <- data.frame(
  Test = c("Overall PERMANOVA", "Dispersion Test"),
  F_value = c(round(permanova_result$F[1], 3), round(permdisp_result$tab$F[1], 3)),
  P_value = c(format.pval(permanova_result$`Pr(>F)`[1], digits = 3), 
              format.pval(permdisp_result$tab$`Pr(>F)`[1], digits = 3)),
  Significance = c(ifelse(permanova_p < 0.001, "***", 
                          ifelse(permanova_p < 0.01, "**",
                                 ifelse(permanova_p < 0.05, "*", "ns"))),
                   ifelse(permdisp_result$tab$`Pr(>F)`[1] < 0.001, "***", 
                          ifelse(permdisp_result$tab$`Pr(>F)`[1] < 0.01, "**",
                                 ifelse(permdisp_result$tab$`Pr(>F)`[1] < 0.05, "*", "ns"))))
)

write.csv(results_summary, "DESeq2_Results/PCoA_statistical_summary_table.csv", row.names = FALSE)

cat("\n=== STATISTICAL ANALYSIS COMPLETED ===\n")
cat("Results saved in: DESeq2_Results/\n")
cat("- PCoA_statistical_results.rds (complete R object)\n")
cat("- PCoA_statistical_summary.txt (text summary)\n")
cat("- PCoA_statistical_summary_table.csv (CSV table)\n")
cat("- PCoA_plot_with_stats.png (enhanced plot)\n")
cat("- PCoA_centroids_plot.png (centroid plot)\n")

library(ggpubr)
Fig1 <- ggarrange(Fig1E, Fig2, 
                  labels = c("A", "B"),
                  ncol = 2, nrow = 1,
                  legend = "none")

Fig1


##########################################################################
#########################################################################
#################### Differential Expression Analysis ####################
#################### ############################### ####################

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# Set up contrast groups
contrast_groups <- list(
  CDF_vs_Control = c("GroupS", "CDF", "Control"),
  DHF_vs_Control = c("GroupS", "DHF", "Control"), 
  DSS_vs_Control = c("GroupS", "DSS", "Control"),
  DHF_vs_CDF = c("GroupS", "DHF", "CDF"),
  DSS_vs_CDF = c("GroupS", "DSS", "CDF"),
  DSS_vs_DHF = c("GroupS", "DSS", "DHF")
)

# Function to perform DE analysis for a contrast
perform_DE_analysis <- function(contrast, dds_object, alpha = 0.05, lfc_threshold = 0) {
  contrast_name <- paste0(contrast[2], "_vs_", contrast[3])
  cat("Analyzing contrast:", contrast_name, "\n")
  
  # Extract results
  res <- results(dds_object, 
                 contrast = contrast,
                 alpha = alpha,
                 lfcThreshold = lfc_threshold,
                 cooksCutoff = TRUE,
                 independentFiltering = TRUE)
  
  # Apply log fold change shrinkage for visualization
  resLFC <- lfcShrink(dds_object, 
                      contrast = contrast, 
                      type = "ashr",
                      res = res)
  
  # Summary of results
  cat("Summary for", contrast_name, ":\n")
  print(summary(res))
  
  # Create results dataframe
  res_df <- as.data.frame(resLFC)
  res_df$gene_id <- rownames(res_df)
  res_df$significant <- ifelse(res_df$padj < alpha & !is.na(res_df$padj), "YES", "NO")
  res_df$regulation <- ifelse(res_df$significant == "YES", 
                              ifelse(res_df$log2FoldChange > 0, "UP", "DOWN"), "NS")
  
  # Order by significance and fold change
  res_df <- res_df[order(res_df$padj, abs(res_df$log2FoldChange), decreasing = c(FALSE, TRUE)), ]
  
  return(list(results = res_df, 
              contrast_name = contrast_name,
              raw_results = res,
              shrunken_results = resLFC))
}

# Perform DE analysis for all contrasts
de_results <- list()
for (contrast_name in names(contrast_groups)) {
  de_results[[contrast_name]] <- perform_DE_analysis(contrast_groups[[contrast_name]], dds)
}

#################### Visualization of DE Results ####################

cat("\n=== CREATING DE VISUALIZATIONS ===\n")

# Create directory for DE results
if (!dir.exists("DESeq2_Results/DE_Analysis")) {
  dir.create("DESeq2_Results/DE_Analysis", recursive = TRUE)
}

# 1. MA Plots for each contrast
ma_plots <- list()
for (contrast_name in names(de_results)) {
  resLFC <- de_results[[contrast_name]]$shrunken_results
  
  ma_plot <- ggplot(as.data.frame(resLFC), aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey", "red"), name = "Significant\n(|LFC|>1, padj<0.05)") +
    scale_x_log10() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
    labs(title = paste("MA Plot:", contrast_name),
         x = "Mean Expression (log10 scale)",
         y = "Log2 Fold Change") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ma_plots[[contrast_name]] <- ma_plot
  ggsave(paste0("DESeq2_Results/DE_Analysis/MA_plot_", contrast_name, ".png"), 
         ma_plot, width = 8, height = 6, dpi = 300)
}

# 2. Volcano Plots for each contrast
volcano_plots <- list()
for (contrast_name in names(de_results)) {
  res_df <- de_results[[contrast_name]]$results
  
  volcano_plot <- EnhancedVolcano(res_df,
                                  lab = rownames(res_df),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = paste('Volcano Plot:', contrast_name),
                                  subtitle = NULL,
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 2.0,
                                  labSize = 3.0,
                                  colAlpha = 0.6,
                                  legendPosition = 'right',
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5,
                                  max.overlaps = 15)
  
  volcano_plots[[contrast_name]] <- volcano_plot
  ggsave(paste0("DESeq2_Results/DE_Analysis/Volcano_plot_", contrast_name, ".png"), 
         volcano_plot, width = 10, height = 8, dpi = 300)
}

# 3. Heatmap of top differentially expressed genes
# Get top DE genes across all comparisons
get_top_genes <- function(de_results, n_top = 50) {
  all_sig_genes <- c()
  for (contrast_name in names(de_results)) {
    res_df <- de_results[[contrast_name]]$results
    sig_genes <- res_df %>% 
      filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
      arrange(padj) %>%
      head(n_top) %>%
      rownames()
    all_sig_genes <- c(all_sig_genes, sig_genes)
  }
  return(unique(all_sig_genes))
}

top_genes <- get_top_genes(de_results, n_top = 30)

if (length(top_genes) > 0) {
  # Extract normalized counts for heatmap
  vsd <- vst(dds, blind = FALSE)
  top_gene_counts <- assay(vsd)[top_genes, ]
  
  # Annotation data
  annotation_df <- as.data.frame(colData(dds)[, "GroupS", drop = FALSE])
  colnames(annotation_df) <- "Group"
  
  # Create heatmap
  heatmap_plot <- pheatmap(top_gene_counts,
                           annotation_col = annotation_df,
                           scale = "row",
                           show_rownames = TRUE,
                           show_colnames = TRUE,
                           cluster_rows = TRUE,
                           cluster_cols = TRUE,
                           treeheight_row = 0,
                           treeheight_col = 20,
                           fontsize_row = 8,
                           fontsize_col = 8,
                           main = "Top Differentially Expressed Genes\n(Z-score normalized)")
  
  # Save heatmap
  png("DESeq2_Results/DE_Analysis/DE_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
  print(heatmap_plot)
  dev.off()
}

#################### Results Summary and Export ####################

cat("\n=== EXPORTING DE RESULTS ===\n")

# Create Excel workbook for results
wb <- createWorkbook()

# Summary sheet with statistics
summary_data <- data.frame(
  Contrast = character(),
  Total_Genes = integer(),
  Significant_Genes = integer(),
  Upregulated = integer(),
  Downregulated = integer(),
  stringsAsFactors = FALSE
)

for (contrast_name in names(de_results)) {
  res_df <- de_results[[contrast_name]]$results
  
  sig_up <- sum(res_df$regulation == "UP", na.rm = TRUE)
  sig_down <- sum(res_df$regulation == "DOWN", na.rm = TRUE)
  sig_total <- sig_up + sig_down
  
  summary_data <- rbind(summary_data, data.frame(
    Contrast = contrast_name,
    Total_Genes = nrow(res_df),
    Significant_Genes = sig_total,
    Upregulated = sig_up,
    Downregulated = sig_down
  ))
  
  # Create formatted sheet for each contrast
  sheet_name <- substr(contrast_name, 1, 31) # Excel sheet name limit
  addWorksheet(wb, sheet_name)
  
  # Prepare data for export - REMOVE 'stat' column since it doesn't exist after shrinkage
  export_df <- res_df %>%
    select(gene_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, significant, regulation) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  
  writeData(wb, sheet = sheet_name, x = export_df)
  
  # Add conditional formatting
  sig_style <- createStyle(fgFill = "#E6F3FF")
  up_style <- createStyle(fgFill = "#FFE6E6") # Light red for upregulated
  down_style <- createStyle(fgFill = "#E6FFE6") # Light green for downregulated
  
  # Apply styles
  sig_rows <- which(export_df$significant == "YES") + 1
  up_rows <- which(export_df$regulation == "UP") + 1
  down_rows <- which(export_df$regulation == "DOWN") + 1
  
  if (length(sig_rows) > 0) {
    addStyle(wb, sheet = sheet_name, style = sig_style, rows = sig_rows, cols = 1:8, gridExpand = TRUE)
  }
  if (length(up_rows) > 0) {
    addStyle(wb, sheet = sheet_name, style = up_style, rows = up_rows, cols = 1:8, gridExpand = TRUE)
  }
  if (length(down_rows) > 0) {
    addStyle(wb, sheet = sheet_name, style = down_style, rows = down_rows, cols = 1:8, gridExpand = TRUE)
  }
  
  # Freeze panes and auto-fit columns
  freezePane(wb, sheet = sheet_name, firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = 1:8, widths = "auto")
}

# Add summary sheet
addWorksheet(wb, "Summary")
writeData(wb, sheet = "Summary", x = summary_data)

# Save workbook
saveWorkbook(wb, "DESeq2_Results/DE_Analysis/Differential_Expression_Results.xlsx", overwrite = TRUE)

# Save individual CSV files
for (contrast_name in names(de_results)) {
  res_df <- de_results[[contrast_name]]$results
  write.csv(res_df, 
            paste0("DESeq2_Results/DE_Analysis/DE_results_", contrast_name, ".csv"), 
            row.names = FALSE)
}

# Save R object for later use
saveRDS(de_results, "DESeq2_Results/DE_Analysis/DE_results_complete.rds")

#################### Functional Analysis Preparation ####################

# Prepare files for functional enrichment analysis
cat("\n=== PREPARING FILES FOR FUNCTIONAL ANALYSIS ===\n")

# Create directories for functional analysis
if (!dir.exists("DESeq2_Results/Functional_Analysis")) {
  dir.create("DESeq2_Results/Functional_Analysis", recursive = TRUE)
}

# Prepare gene lists for GSEA and ORA
for (contrast_name in names(de_results)) {
  res_df <- de_results[[contrast_name]]$results
  
  # For GSEA: ranked list by log2FoldChange instead of stat
  gsea_df <- res_df %>%
    filter(!is.na(log2FoldChange)) %>%
    select(gene_id, log2FoldChange) %>%
    arrange(desc(log2FoldChange))
  
  write.table(gsea_df, 
              paste0("DESeq2_Results/Functional_Analysis/", contrast_name, "_GSEA_ranked.rnk"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # For ORA: significant up and down regulated genes
  sig_up <- res_df %>% 
    filter(regulation == "UP") %>% 
    pull(gene_id)
  sig_down <- res_df %>% 
    filter(regulation == "DOWN") %>% 
    pull(gene_id)
  
  writeLines(sig_up, 
             paste0("DESeq2_Results/Functional_Analysis/", contrast_name, "_UP_genes.txt"))
  writeLines(sig_down, 
             paste0("DESeq2_Results/Functional_Analysis/", contrast_name, "_DOWN_genes.txt"))
}
#################### Final Summary ####################

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED ===\n")
cat("Results saved in: DESeq2_Results/DE_Analysis/\n")
cat("\nSummary of significant DE genes:\n")
print(summary_data)

cat("\nFiles created:\n")
cat("- Differential_Expression_Results.xlsx (Complete results with formatting)\n")
cat("- Individual CSV files for each contrast\n")
cat("- MA plots and Volcano plots for each contrast\n")
cat("- Heatmap of top DE genes\n")
cat("- Files for functional analysis (GSEA and ORA)\n")
cat("- DE_results_complete.rds (R object for later use)\n")

# Print session info for reproducibility
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())



#################### Create Single PDF with All Volcano Plots_v3 ####################

cat("\n=== CREATING COMBINED VOLCANO PLOTS PDF ===\n")

library(grid)
library(gridExtra)

# Create a single PDF with all 6 volcano plots on ONE PAGE
pdf("DESeq2_Results/DE_Analysis/4All_Volcano_Plots.pdf", width = 16, height = 20)

# Create a list of all volcano plots
volcano_plot_list <- list()

# Define the labels for each plot
plot_labels <- c("A", "B", "C", "D", "E", "F")

for (i in seq_along(names(de_results))) {
  contrast_name <- names(de_results)[i]
  label <- plot_labels[i]
  
  res_df <- de_results[[contrast_name]]$results
  
  # Get contrast info for subtitle
  total_genes <- nrow(res_df)
  sig_genes <- sum(res_df$significant == "YES", na.rm = TRUE)
  up_genes <- sum(res_df$regulation == "UP", na.rm = TRUE)
  down_genes <- sum(res_df$regulation == "DOWN", na.rm = TRUE)
  
  # Create enhanced volcano plot with panel label
  volcano_plot <- EnhancedVolcano(res_df,
                                  lab = rownames(res_df),
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = paste0(label, ": ", contrast_name),
                                  subtitle = paste('Total:', total_genes, 
                                                   '| Sig:', sig_genes,
                                                   '| Up:', up_genes,
                                                   '| Down:', down_genes),
                                  pCutoff = 0.05,
                                  FCcutoff = 1,
                                  pointSize = 1.5,
                                  labSize = 2.5,
                                  colAlpha = 0.6,
                                  legendPosition = 'right',
                                  legendLabSize = 12,
                                  legendIconSize = 4,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.3,
                                  max.overlaps = 8,
                                  caption = NULL,
                                  labFace = 'bold',
                                  boxedLabels = FALSE,
                                  xlab = bquote(~Log[2]~ "fold change"),
                                  ylab = bquote(~-Log[10]~ "adjusted p-value"))
  
  # Remove unnecessary padding from individual plots
  volcano_plot <- volcano_plot + 
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5))
  
  volcano_plot_list[[contrast_name]] <- volcano_plot
}

# Arrange all plots in a single page
grid.arrange(
  grobs = volcano_plot_list,
  ncol = 2,
  nrow = 3,
  layout_matrix = matrix(1:6, nrow = 3, byrow = TRUE),
  top = textGrob("Differential Gene Expression Analysis: Dengue Severity Comparisons",
                 gp = gpar(fontsize = 20, fontface = "bold", col = "black")),
  bottom = textGrob(paste("Analysis date:", format(Sys.Date(), "%B %d, %Y"), 
                          "| Significance thresholds: |log2FC| > 1, padj < 0.05",
                          "\nSamples: Control (n=6), CDF (n=6), DHF (n=6), DSS (n=6)"),
                    gp = gpar(fontsize = 10, fontface = "italic")),
  padding = unit(0.5, "line")
)

dev.off()

cat("✓ Combined PDF with all volcano plots saved as: DESeq2_Results/DE_Analysis/3All_Volcano_Plots.pdf\n")
cat("✓ Page dimensions: 16 x 20 inches\n")
cat("✓ Layout: 3 rows x 2 columns (6 plots total)\n")
cat("✓ All plots arranged on a single page\n")
cat("✓ Panel labels added: A, B, C, D, E, F\n")