# Set working directory
setwd("C:/Users/USER/Desktop/D_BCSIR/CIBERSORTx")

# Load required libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(gridExtra)
library(ggplot2)
library(ggthemes)


# Read data
count_data <- read.xlsx("CIBERSORTx.xlsx", sheet = "cell", colNames = TRUE, rowNames = TRUE)
metadata <- read.xlsx("CIBERSORTx.xlsx", sheet = "m", colNames = TRUE, rowNames = TRUE)

# Check the structure of metadata
cat("Metadata column names:\n")
print(colnames(metadata))
cat("\nMetadata structure:\n")
print(head(metadata))
cat("\n")

# Convert count data to matrix if needed
count_data_matrix <- as.matrix(count_data)

# 1. Prepare data for visualization
# Add sample names if not already as rownames
if (!"Sample" %in% colnames(count_data)) {
  count_data$Sample <- rownames(count_data)
}

# Melt data for ggplot
melted_data <- melt(count_data, id.vars = "Sample", 
                    variable.name = "CellType", 
                    value.name = "Proportion")

# Merge with metadata - check all possible group column names
cat("Available columns in metadata:\n")
print(colnames(metadata))
cat("\n")

# Debug: Check what group column actually exists
if ("Group" %in% colnames(metadata)) {
  cat("Found 'Group' column in metadata\n")
} else if ("GroupS" %in% colnames(metadata)) {
  cat("Found 'GroupS' column in metadata\n")
} else if ("Group." %in% colnames(metadata)) {
  cat("Found 'Group.' column in metadata\n")
} else {
  cat("No standard group column found. Available columns:\n")
  print(colnames(metadata))
}

# Merge the data
melted_data <- merge(melted_data, metadata, by.x = "Sample", by.y = "row.names")

# Debug: Check what group columns are available after merge
cat("\nColumns in melted_data after merge:\n")
print(colnames(melted_data))
cat("\n")

# Check for the group column - try different possibilities
if ("Group" %in% colnames(melted_data)) {
  cat("Using 'Group' column\n")
  melted_data$GroupFactor <- melted_data$Group
} else if ("GroupS" %in% colnames(melted_data)) {
  cat("Using 'GroupS' column\n")
  melted_data$GroupFactor <- melted_data$GroupS
} else if ("Group." %in% colnames(melted_data)) {
  cat("Using 'Group.' column\n")
  melted_data$GroupFactor <- melted_data$Group.
} else {
  # Try to find any column containing "Group"
  group_cols <- grep("Group", colnames(melted_data), value = TRUE, ignore.case = TRUE)
  if (length(group_cols) > 0) {
    cat(paste("Using first group column found:", group_cols[1], "\n"))
    melted_data$GroupFactor <- melted_data[[group_cols[1]]]
  } else {
    # If no group column found, create a dummy group
    cat("No group column found. Creating default groups.\n")
    melted_data$GroupFactor <- "Unknown"
  }
}

# Check unique values in GroupFactor
cat("\nUnique values in GroupFactor:\n")
print(unique(melted_data$GroupFactor))
cat("\n")

# Clean up group names (remove any trailing dots or spaces)
melted_data$GroupFactor <- trimws(gsub("\\.$", "", as.character(melted_data$GroupFactor)))

# Define the correct order based on your data
group_levels <- c("Control", "CDF", "DHF", "DSS")
melted_data$GroupFactor <- factor(melted_data$GroupFactor, levels = group_levels)

# Create a custom color palette with distinct colors for all cell types
cell_types <- unique(melted_data$CellType)
n_colors <- length(cell_types)

# Create a color palette
if (n_colors <= 12) {
  cell_colors <- brewer.pal(n_colors, "Set3")
} else {
  set3_colors <- brewer.pal(12, "Set3")
  set2_colors <- brewer.pal(8, "Set2")
  set1_colors <- brewer.pal(9, "Set1")
  cell_colors <- c(set3_colors, set2_colors, set1_colors)[1:n_colors]
}
names(cell_colors) <- cell_types

# 2. Sample-wise stacked bar plot with legend on right side
sample_barplot <- ggplot(melted_data, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  facet_grid(~ Group, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    panel.spacing = unit(0.2, "lines"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",  # Legend on right side
    legend.box = "vertical",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.border = element_rect(color = "grey80", fill = NA, size = 0.5)
  ) +
  guides(fill = guide_legend(ncol = 1)) +  # Single column legend
  labs(title = "Sample-wise Immune Cell Composition",
       x = "Sample", 
       y = "Cell Proportion",
       fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = cell_colors)

print(sample_barplot)

# 3. Group-wise average bar plot with legend on right side
group_summary <- melted_data %>%
  group_by(Group, CellType) %>%
  summarise(Mean_Proportion = mean(Proportion, na.rm = TRUE),
            SE_Proportion = sd(Proportion, na.rm = TRUE)/sqrt(n())) %>%
  ungroup()

group_barplot <- ggplot(group_summary, aes(x = Group, y = Mean_Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "right",  # Legend on right side
    legend.box = "vertical",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 1)) +  # Single column legend
  labs(title = "Group-wise Average Immune Cell Composition",
       x = "Disease Group", 
       y = "Mean Cell Proportion",
       fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = cell_colors)

print(group_barplot)

# 4. Save both bar plots to a single PDF
pdf("immune_cell_barplots.pdf", width = 18, height = 14)
# Arrange plots with appropriate spacing
grid.arrange(
  sample_barplot + 
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm")),  # Extra space on right for legend
  
  group_barplot + 
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm")),  # Extra space on right for legend
  
  nrow = 2,
  heights = c(1.2, 0.8)
)
dev.off()

cat("Bar plots saved to 'immune_cell_barplots.pdf'\n")



# 5. Continue with other visualizations_V1
# Group-wise box plots for major cell types
top_cells <- melted_data %>%
  group_by(CellType) %>%
  summarise(Mean_Prop = mean(Proportion)) %>%
  arrange(desc(Mean_Prop)) %>%
  slice_head(n = 22) %>%
  pull(CellType)

melted_data_top <- melted_data %>% filter(CellType %in% top_cells)

# Define group colors
group_colors <- c("Control" = "#66c2a5", 
                  "CDF" = "#fc8d62", 
                  "DHF" = "#8da0cb", 
                  "DSS" = "#e78ac3")

cell_boxplot <- ggplot(melted_data_top, aes(x = GroupS, y = Proportion, fill = GroupS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Distribution of Top 10 Immune Cell Types by Group",
       x = "Disease Group", 
       y = "Cell Proportion",
       fill = "Group") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(method = "kruskal.test", label = "p.format", 
                     label.y.npc = "top", size = 3)

print(cell_boxplot)



top_cells <- melted_data %>%
  group_by(CellType) %>%
  summarise(Mean_Prop = mean(Proportion)) %>%
  arrange(desc(Mean_Prop)) %>%
  slice_head(n = 22) %>%
  pull(CellType)

melted_data_top <- melted_data %>% filter(CellType %in% top_cells)

# Define group colors
group_colors <- c("Control" = "#66c2a5", 
                  "CDF" = "#fc8d62", 
                  "DHF" = "#8da0cb", 
                  "DSS" = "#e78ac3")

# Create boxplot with statistical significance asterisks
cell_boxplot <- ggplot(melted_data_top, aes(x = GroupS, y = Proportion, fill = GroupS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Distribution of Top 22 Immune Cell Types by Group",
       x = "Disease Group", 
       y = "Cell Proportion",
       fill = "Group") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = group_colors) +
  # Add pairwise comparisons with asterisks for significance
  stat_compare_means(
    method = "kruskal.test", 
    label = "p.signif",  # Use p.signif for asterisks
    hide.ns = TRUE,      # Hide non-significant comparisons
    label.y.npc = "top", 
    size = 4,
    vjust = 0.5,
    step.increase = 0.1  # Space between comparison brackets
  )

print(cell_boxplot)



################ V2  ##########################################
# Load required libraries
library(ggthemes)
library(ggpubr)
library(dplyr)
library(tidyr)

# Create pairwise comparisons for all groups
comparison_pairs <- list(
  c("Control", "CDF"),
  c("Control", "DHF"), 
  c("Control", "DSS"),
  c("CDF", "DHF"),
  c("CDF", "DSS"),
  c("DHF", "DSS")
)

# PART 1: Perform pairwise tests and export results
# Create a function to perform Wilcoxon test with error handling
perform_wilcox_test <- function(data, group1, group2) {
  group1_data <- data %>% filter(GroupS == group1) %>% pull(Proportion)
  group2_data <- data %>% filter(GroupS == group2) %>% pull(Proportion)
  
  # Check if we have enough data
  if (length(group1_data) < 2 || length(group2_data) < 2) {
    return(list(p_value = NA, statistic = NA))
  }
  
  # Perform Wilcoxon test with error handling
  result <- tryCatch({
    wilcox.test(group1_data, group2_data, exact = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(result)) {
    return(list(p_value = NA, statistic = NA))
  }
  
  return(list(p_value = result$p.value, statistic = result$statistic))
}

# Perform all pairwise tests for each cell type and store results
pairwise_results <- list()

for (cell in unique(melted_data_top$CellType)) {
  cell_data <- melted_data_top %>% filter(CellType == cell)
  
  cell_results <- data.frame()
  
  for (pair in comparison_pairs) {
    test_result <- perform_wilcox_test(cell_data, pair[1], pair[2])
    
    result_row <- data.frame(
      CellType = cell,
      Group1 = pair[1],
      Group2 = pair[2],
      Comparison = paste(pair[1], "vs", pair[2]),
      W_statistic = ifelse(is.na(test_result$statistic), NA, test_result$statistic),
      P_value = test_result$p_value,
      Significant = ifelse(is.na(test_result$p_value), FALSE, test_result$p_value < 0.05),
      Significance_level = case_when(
        is.na(test_result$p_value) ~ "NA",
        test_result$p_value < 0.001 ~ "***",
        test_result$p_value < 0.01 ~ "**",
        test_result$p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      )
    )
    
    cell_results <- rbind(cell_results, result_row)
  }
  
  pairwise_results[[cell]] <- cell_results
}

# Combine all results into one dataframe
all_pairwise_results <- do.call(rbind, pairwise_results)

# Filter for significant results only
significant_results <- all_pairwise_results %>%
  filter(Significant == TRUE) %>%
  arrange(CellType, P_value)

# Export pairwise test results
write.table(all_pairwise_results, "pairwise_wilcoxon_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(significant_results, "significant_pairwise_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Also export as CSV for easier use
write.csv(all_pairwise_results, "pairwise_wilcoxon_results.csv", row.names = FALSE)
write.csv(significant_results, "significant_pairwise_results.csv", row.names = FALSE)

cat("\n=== Pairwise Test Results Exported ===\n")
cat("1. pairwise_wilcoxon_results.txt - All pairwise comparisons\n")
cat("2. significant_pairwise_results.txt - Significant comparisons only\n")
cat("3. CSV versions also created for easier import\n")

# PART 2: Create box plots showing only significant comparisons

# OPTION A: All significant comparisons
# First, get unique significant comparisons across all cell types
all_sig_comparisons <- unique(significant_results$Comparison)

# Create a mapping of comparison string to vector
sig_comparison_list <- list()
for (comp in all_sig_comparisons) {
  groups <- strsplit(comp, " vs ")[[1]]
  sig_comparison_list[[comp]] <- c(groups[1], groups[2])
}

# Convert to list format for stat_compare_means
sig_comparisons <- unname(sig_comparison_list)

# Create plot with ONLY significant comparisons
if (length(sig_comparisons) > 0) {
  cell_boxplot_all_sig <- ggplot(melted_data_top, aes(x = GroupS, y = Proportion, fill = GroupS)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 16) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
    theme_few() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(title = "Distribution of Immune Cell Types by Disease Group",
         subtitle = "Showing only statistically significant pairwise comparisons (Wilcoxon test)",
         x = "Disease Group", 
         y = "Cell Proportion",
         fill = "Group") +
    scale_y_continuous(labels = scales::percent_format(scale = 100, accuracy = 0.1)) +
    scale_fill_manual(values = group_colors) +
    # Add only significant comparisons
    stat_compare_means(
      method = "wilcox.test",
      comparisons = sig_comparisons,
      label = "p.signif",
      hide.ns = TRUE,
      size = 4.5,
      tip.length = 0.02,
      bracket.size = 0.6,
      step.increase = 0.08,
      vjust = 0.5,
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      )
    )
  
  print(cell_boxplot_all_sig)
} else {
  cat("\nNo significant pairwise comparisons found across any cell types.\n")
  # Create plot without comparisons
  cell_boxplot_all_sig <- ggplot(melted_data_top, aes(x = GroupS, y = Proportion, fill = GroupS)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 16) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
    theme_few() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "red"),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(title = "Distribution of Immune Cell Types by Disease Group",
         subtitle = "No statistically significant pairwise comparisons found (p < 0.05)",
         x = "Disease Group", 
         y = "Cell Proportion",
         fill = "Group") +
    scale_y_continuous(labels = scales::percent_format(scale = 100, accuracy = 0.1)) +
    scale_fill_manual(values = group_colors)
  
  print(cell_boxplot_all_sig)
}

# OPTION B: Create individual plots for each cell type showing their specific significant comparisons
# This ensures each facet only shows comparisons that are significant for that specific cell type

# Create a list to store plots for each cell type
cell_specific_plots <- list()

for (cell in unique(melted_data_top$CellType)) {
  # Get significant comparisons for this cell type
  cell_sig_results <- significant_results %>% 
    filter(CellType == cell)
  
  # Convert to list format for stat_compare_means
  cell_sig_comparisons <- list()
  if (nrow(cell_sig_results) > 0) {
    for (i in 1:nrow(cell_sig_results)) {
      cell_sig_comparisons[[i]] <- c(cell_sig_results$Group1[i], cell_sig_results$Group2[i])
    }
  }
  
  # Get cell data
  cell_data <- melted_data_top %>% filter(CellType == cell)
  
  # Create plot
  p <- ggplot(cell_data, aes(x = GroupS, y = Proportion, fill = GroupS)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7, shape = 16) +
    theme_few() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(title = cell,
         y = "Proportion") +
    scale_y_continuous(labels = scales::percent_format(scale = 100, accuracy = 0.1)) +
    scale_fill_manual(values = group_colors)
  
  # Add significant comparisons if any
  if (length(cell_sig_comparisons) > 0) {
    p <- p + stat_compare_means(
      comparisons = cell_sig_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = TRUE,
      size = 4.5,
      tip.length = 0.02,
      bracket.size = 0.6,
      step.increase = 0.1,
      vjust = 0.5
    )
  }
  
  cell_specific_plots[[cell]] <- p
}

# Combine all cell-specific plots
library(patchwork)

if (length(cell_specific_plots) > 0) {
  combined_cell_plots <- wrap_plots(cell_specific_plots, ncol = 5) +
    plot_annotation(
      title = "Cell Type-Specific Significant Comparisons",
      subtitle = "Each panel shows only comparisons significant for that specific cell type",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )
  
  print(combined_cell_plots)
}

# OPTION C: Create a summary plot showing only the MOST significant comparison per cell type
# This reduces clutter

top_sig_comparisons <- significant_results %>%
  group_by(CellType) %>%
  slice_min(P_value, n = 1) %>%
  ungroup()

# Convert to list format
top_sig_list <- list()
if (nrow(top_sig_comparisons) > 0) {
  for (i in 1:nrow(top_sig_comparisons)) {
    top_sig_list[[i]] <- c(top_sig_comparisons$Group1[i], top_sig_comparisons$Group2[i])
  }
  
  cell_boxplot_top_sig <- ggplot(melted_data_top, aes(x = GroupS, y = Proportion, fill = GroupS)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, shape = 16) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 5) +
    theme_few() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(title = "Distribution of Immune Cell Types by Disease Group",
         subtitle = "Showing only the most significant pairwise comparison per cell type",
         x = "Disease Group", 
         y = "Cell Proportion",
         fill = "Group") +
    scale_y_continuous(labels = scales::percent_format(scale = 100, accuracy = 0.1)) +
    scale_fill_manual(values = group_colors) +
    # Add only the top significant comparison per cell type
    stat_compare_means(
      method = "wilcox.test",
      comparisons = top_sig_list,
      label = "p.signif",
      hide.ns = TRUE,
      size = 4.5,
      tip.length = 0.02,
      bracket.size = 0.6,
      step.increase = 0.08,
      vjust = 0.5
    )
  
  print(cell_boxplot_top_sig)
}

# Save all plots
if (exists("cell_boxplot_all_sig")) {
  ggsave("celltype_boxplots_all_significant.pdf", cell_boxplot_all_sig, 
         width = 22, height = 18, dpi = 300, bg = "white")
}

if (exists("combined_cell_plots")) {
  ggsave("celltype_boxplots_cell_specific.pdf", combined_cell_plots, 
         width = 24, height = 20, dpi = 300, bg = "white")
}

if (exists("cell_boxplot_top_sig")) {
  ggsave("celltype_boxplots_top_significant.pdf", cell_boxplot_top_sig, 
         width = 22, height = 18, dpi = 300, bg = "white")
}

# Save all plots in a PDF
pdf("celltype_significant_comparisons.pdf", width = 16, height = 12)
if (exists("cell_boxplot_all_sig")) print(cell_boxplot_all_sig)
if (exists("combined_cell_plots")) print(combined_cell_plots)
if (exists("cell_boxplot_top_sig")) print(cell_boxplot_top_sig)
dev.off()

# Create a summary statistics table
summary_stats <- melted_data_top %>%
  group_by(CellType, GroupS) %>%
  summarise(
    Mean = mean(Proportion, na.rm = TRUE),
    Median = median(Proportion, na.rm = TRUE),
    SD = sd(Proportion, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    Min = min(Proportion, na.rm = TRUE),
    Max = max(Proportion, na.rm = TRUE),
    N = n(),
    .groups = 'drop'
  )

# Create a comprehensive significance summary table
sig_summary <- sig_results %>%
  mutate(
    Significance = case_when(
      Kruskal_P < 0.001 ~ "*** (p < 0.001)",
      Kruskal_P < 0.01 ~ "** (p < 0.01)",
      Kruskal_P < 0.05 ~ "* (p < 0.05)",
      TRUE ~ "Not Significant"
    ),
    `FDR Adjusted P` = p.adjust(Kruskal_P, method = "BH")
  ) %>%
  select(
    `Cell Type` = CellType,
    `Kruskal-Wallis P` = Kruskal_P,
    `FDR Adjusted P` = `FDR Adjusted P`,
    Significance,
    `Significance Level` = Sig_Level
  ) %>%
  arrange(`Kruskal-Wallis P`)

# Save significance summary
write.csv(sig_summary, "celltype_significance_summary1.csv", row.names = FALSE)


# Export summary statistics
write.table(summary_stats, "celltype_summary_statistics2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(summary_stats, "celltype_summary_statistics.csv", row.names = FALSE)

# Print summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\nFiles created:\n")
cat("1. pairwise_wilcoxon_results.txt - All pairwise test results\n")
cat("2. significant_pairwise_results.txt - Significant results only\n")
cat("3. celltype_summary_statistics.txt - Summary statistics\n")
cat("4. celltype_boxplots_all_significant.pdf - All significant comparisons\n")
cat("5. celltype_boxplots_cell_specific.pdf - Cell-specific comparisons\n")
cat("6. celltype_boxplots_top_significant.pdf - Top significant comparisons\n")
cat("7. celltype_significant_comparisons.pdf - All plots in PDF\n")
cat("\nSignificant comparisons found:", nrow(significant_results), "\n")
cat("Cell types with significant differences:", length(unique(significant_results$CellType)), "\n")
cat("\nSignificance codes:\n")
cat("*** p < 0.001\n")
cat("**  p < 0.01\n")
cat("*   p < 0.05\n")
cat("NS  Not significant (not shown on plots)\n")







# 6. Heatmap visualization
# Prepare data for heatmap
heatmap_data <- as.matrix(count_data[, -which(colnames(count_data) == "Sample")])
rownames(heatmap_data) <- rownames(count_data)

# Create annotation for samples - use GroupFactor from metadata
# First, let's create a mapping from Sample to GroupFactor
group_mapping <- unique(melted_data[, c("Sample", "Group")])
rownames(group_mapping) <- group_mapping$Sample

# Create annotation for samples
annotation_df <- data.frame(
  Group = group_mapping[rownames(heatmap_data), "Group"],
  Age = metadata[rownames(heatmap_data), "Age"],
  Sex = metadata[rownames(heatmap_data), "Sex"]
)
rownames(annotation_df) <- rownames(heatmap_data)

# Create heatmap
pdf("immune_cell_heatmap1.pdf", width = 14, height = 10)
pheatmap(t(heatmap_data),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = annotation_df,
         annotation_colors = list(GroupS = group_colors),
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         fontsize = 10,
         main = "Immune Cell Composition Heatmap (Z-score normalized)")
dev.off()

cat("Heatmap saved to 'immune_cell_heatmap.pdf'\n")






# 7. Correlation heatmap between cell types
cor_matrix <- cor(heatmap_data, use = "complete.obs")

pdf("cell_correlation_heatmap.pdf", width = 12, height = 10)
pheatmap(cor_matrix,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         fontsize_col = 8,
         fontsize = 10,
         main = "Correlation Between Immune Cell Types")
dev.off()

cat("Correlation heatmap saved to 'cell_correlation_heatmap.pdf'\n")

# 9. Save individual plots as PNG
ggsave("sample_stacked_barplot.png", sample_barplot, width = 18, height = 10, dpi = 300)
ggsave("group_stacked_barplot.png", group_barplot, width = 14, height = 8, dpi = 300)
ggsave("celltype_boxplots.png", cell_boxplot, width = 16, height = 10, dpi = 300)

# Save boxplot as PDF
pdf("celltype_boxplots.pdf", width = 14, height = 10)
print(cell_boxplot)
dev.off()


# 10. Create a summary statistics table
summary_stats <- melted_data %>%
  group_by(GroupFactor, CellType) %>%
  summarise(
    Mean = mean(Proportion, na.rm = TRUE),
    SD = sd(Proportion, na.rm = TRUE),
    Min = min(Proportion, na.rm = TRUE),
    Max = max(Proportion, na.rm = TRUE),
    N = n(),
    .groups = 'drop'
  ) %>%
  arrange(GroupFactor, desc(Mean))

# Save summary statistics
write.csv(summary_stats, "immune_cell_summary_statistics.csv", row.names = FALSE)
cat("Summary statistics saved to 'immune_cell_summary_statistics.csv'\n")

# Print session info for reproducibility
cat("\n=== Analysis Complete ===\n")
cat("All visualizations have been saved to PDF and PNG files.\n")
cat("Check the current directory for:\n")
cat("1. immune_cell_barplots.pdf - Both stacked bar plots with legends on right\n")
cat("2. immune_cell_heatmap.pdf - Heatmap of cell composition\n")
cat("3. cell_correlation_heatmap.pdf - Correlation between cell types\n")
cat("4. celltype_boxplots.pdf - Box plots of top 10 cell types\n")
cat("5. pca_analysis.pdf - PCA plot\n")
cat("6. immune_cell_summary_statistics.csv - Summary statistics\n")
cat("\nDebug Info:\n")
cat("- Number of samples:", length(unique(melted_data$Sample)), "\n")
cat("- Number of cell types:", length(unique(melted_data$CellType)), "\n")
cat("- Groups found:", paste(unique(melted_data$GroupFactor), collapse = ", "), "\n")









# 3. Group-wise average bar plot with legend on right side
# Make sure we're using GroupS (not Group)
group_summary <- melted_data %>%
  group_by(GroupS, CellType) %>%  # Changed from Group to GroupS
  summarise(Mean_Proportion = mean(Proportion, na.rm = TRUE),
            SE_Proportion = sd(Proportion, na.rm = TRUE)/sqrt(n()),
            N = n()) %>%
  ungroup()

# Export the group summary percentages
group_percentages <- group_summary %>%
  mutate(
    Percentage = paste0(round(Mean_Proportion * 100, 2), "%"),
    Mean_Percent = round(Mean_Proportion * 100, 2),
    SE_Percent = round(SE_Proportion * 100, 2)
  ) %>%
  select(
    Group = GroupS,
    CellType,
    Mean = Mean_Percent,
    SE = SE_Percent,
    Percentage,
    N
  ) %>%
  arrange(Group, desc(Mean))

# Print the summary
cat("\n=== Group-wise Average Immune Cell Composition (Percentages) ===\n")
print(group_percentages, n = Inf)

# Export to files
write.csv(group_percentages, "group_immune_cell_percentages.csv", row.names = FALSE)
write.table(group_percentages, "group_immune_cell_percentages.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create a formatted table for publication
library(knitr)

group_table_formatted <- group_percentages %>%
  mutate(
    Mean_SE = sprintf("%.1f ± %.1f", Mean, SE),
    Group = factor(Group, levels = c("Control", "CDF", "DHF", "DSS"))
  ) %>%
  select(Group, CellType, `Mean ± SE (%)` = Mean_SE, N) %>%
  pivot_wider(
    names_from = Group, 
    values_from = c(`Mean ± SE (%)`, N),
    names_sep = "_"
  )

# Check the column names that were actually created
cat("Column names after pivot_wider:\n")
print(colnames(group_table_formatted))

# Now select the columns based on what was actually created
# The column names will be in the format: "Mean ± SE (%)_GroupName" and "N_GroupName"
group_table_formatted <- group_table_formatted %>%
  select(
    CellType,
    `Control (%)` = `Mean ± SE (%)_Control`,
    `Control (N)` = `N_Control`,
    `CDF (%)` = `Mean ± SE (%)_CDF`,
    `CDF (N)` = `N_CDF`,
    `DHF (%)` = `Mean ± SE (%)_DHF`,
    `DHF (N)` = `N_DHF`,
    `DSS (%)` = `Mean ± SE (%)_DSS`,
    `DSS (N)` = `N_DSS`
  )

# Alternative simpler approach if the above doesn't work:
group_table_formatted_simple <- group_percentages %>%
  mutate(
    Mean_SE = sprintf("%.1f ± %.1f", Mean, SE),
    Group = factor(Group, levels = c("Control", "CDF", "DHF", "DSS"))
  ) %>%
  select(Group, CellType, `Mean ± SE (%)` = Mean_SE, N) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(`Mean ± SE (%)`, N)
  )

# Check column names
cat("\nColumn names in simple version:\n")
print(colnames(group_table_formatted_simple))

# Rename columns
if ("Mean ± SE (%)_Control" %in% colnames(group_table_formatted_simple)) {
  # Use underscore format
  group_table_formatted_simple <- group_table_formatted_simple %>%
    rename(
      `Control (%)` = `Mean ± SE (%)_Control`,
      `Control (N)` = `N_Control`,
      `CDF (%)` = `Mean ± SE (%)_CDF`,
      `CDF (N)` = `N_CDF`,
      `DHF (%)` = `Mean ± SE (%)_DHF`,
      `DHF (N)` = `N_DHF`,
      `DSS (%)` = `Mean ± SE (%)_DSS`,
      `DSS (N)` = `N_DSS`
    )
} else if ("Control_Mean ± SE (%)" %in% colnames(group_table_formatted_simple)) {
  # Use different naming pattern
  group_table_formatted_simple <- group_table_formatted_simple %>%
    rename(
      `Control (%)` = `Control_Mean ± SE (%)`,
      `Control (N)` = `Control_N`,
      `CDF (%)` = `CDF_Mean ± SE (%)`,
      `CDF (N)` = `CDF_N`,
      `DHF (%)` = `DHF_Mean ± SE (%)`,
      `DHF (N)` = `DHF_N`,
      `DSS (%)` = `DSS_Mean ± SE (%)`,
      `DSS (N)` = `DSS_N`
    )
}

# Export formatted table
write.csv(group_table_formatted_simple, "group_immune_cell_table_formatted.csv", row.names = FALSE)

cat("\nFormatted table saved to: group_immune_cell_table_formatted.csv\n")
cat("\nFirst few rows of the formatted table:\n")
print(head(group_table_formatted_simple))



# Create the bar plot (using GroupS)
group_barplot <- ggplot(group_summary, aes(x = GroupS, y = Mean_Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    legend.position = "right",  # Legend on right side
    legend.box = "vertical",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  guides(fill = guide_legend(ncol = 1)) +  # Single column legend
  labs(title = "Group-wise Average Immune Cell Composition",
       subtitle = "Based on GroupS classification",
       x = "Disease Group", 
       y = "Mean Cell Proportion",
       fill = "Cell Type") +
  scale_y_continuous(labels = scales::percent_format(scale = 100),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("Control", "CDF", "DHF", "DSS")) +
  scale_fill_manual(values = cell_colors)

print(group_barplot)

# Save the plot
ggsave("group_immune_cell_composition.png", group_barplot, width = 12, height = 8, dpi = 300, bg = "white")

# Additional: Create a summary table with percentages that sum to 100% for each group
group_totals <- group_summary %>%
  group_by(GroupS) %>%
  summarise(Total_Proportion = sum(Mean_Proportion))

# Check if totals are close to 1 (100%)
cat("\n=== Checking Group Totals ===\n")
print(group_totals)

# Create a percentage distribution table
percentage_distribution <- group_percentages %>%
  select(Group, CellType, Mean) %>%
  pivot_wider(names_from = Group, values_from = Mean, values_fill = 0) %>%
  arrange(desc(Control))  # Sort by Control percentage

# Export percentage distribution
write.csv(percentage_distribution, "immune_cell_percentage_distribution.csv", row.names = FALSE)

cat("\nPercentage distribution saved to: immune_cell_percentage_distribution.csv\n")

# Create a heatmap of percentage distribution
library(pheatmap)

# Prepare data for heatmap
heatmap_matrix <- as.matrix(percentage_distribution[, -1])
rownames(heatmap_matrix) <- percentage_distribution$CellType

# Create heatmap
heatmap_plot <- pheatmap(heatmap_matrix,
                         scale = "none",
                         clustering_distance_rows = "euclidean",
                         clustering_distance_cols = "euclidean",
                         clustering_method = "complete",
                         color = colorRampPalette(c("white", "yellow", "red"))(100),
                         display_numbers = TRUE,
                         number_format = "%.1f",
                         fontsize_number = 8,
                         main = "Immune Cell Percentage Distribution by Group",
                         fontsize_row = 8,
                         fontsize_col = 10,
                         angle_col = 45)

# Save heatmap
ggsave("immune_cell_percentage_heatmap.png", heatmap_plot, width = 10, height = 12, dpi = 300, bg = "white")

# Create a summary of top 5 cell types for each group
top_cells_per_group <- group_percentages %>%
  group_by(Group) %>%
  arrange(desc(Mean)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(Group, CellType, Mean, Percentage)

cat("\n=== Top 5 Immune Cell Types per Group ===\n")
print(top_cells_per_group, n = Inf)

# Export top cells
write.csv(top_cells_per_group, "top_immune_cells_per_group.csv", row.names = FALSE)

# Create a summary statistics file
summary_stats <- melted_data %>%
  group_by(GroupS) %>%
  summarise(
    Total_Samples = n_distinct(Sample),
    Mean_Cells_per_Sample = mean(tapply(Proportion, Sample, sum, na.rm = TRUE), na.rm = TRUE),
    Total_Cell_Types = n_distinct(CellType)
  )

cat("\n=== Group Summary Statistics ===\n")
print(summary_stats)

# Export all statistics
write.csv(summary_stats, "group_summary_statistics.csv", row.names = FALSE)

# Create a comprehensive report
sink("immune_cell_analysis_report.txt")
cat("===========================================\n")
cat("IMMUNE CELL COMPOSITION ANALYSIS REPORT\n")
cat("===========================================\n\n")

cat("Date: ", as.character(Sys.Date()), "\n")
cat("Total Samples: ", n_distinct(melted_data$Sample), "\n")
cat("Total Cell Types: ", n_distinct(melted_data$CellType), "\n\n")

cat("GROUP DISTRIBUTION:\n")
cat("-------------------\n")
for (group in c("Control", "CDF", "DHF", "DSS")) {
  n_samples <- melted_data %>% filter(GroupS == group) %>% distinct(Sample) %>% nrow()
  cat(sprintf("%-10s: %2d samples\n", group, n_samples))
}
cat("\n")

cat("AVERAGE CELL COMPOSITION BY GROUP (Top 3 per group):\n")
cat("----------------------------------------------------\n")
for (group in c("Control", "CDF", "DHF", "DSS")) {
  cat("\n", group, ":\n")
  top3 <- group_percentages %>% 
    filter(Group == group) %>% 
    arrange(desc(Mean)) %>% 
    slice_head(n = 3)
  
  for (i in 1:nrow(top3)) {
    cat(sprintf("  %-25s: %6.2f%%\n", 
                top3$CellType[i], 
                top3$Mean[i]))
  }
}

cat("\n\nFILES GENERATED:\n")
cat("----------------\n")
cat("1. group_immune_cell_percentages.csv - Complete percentage data\n")
cat("2. group_immune_cell_table_formatted.csv - Formatted table for publication\n")
cat("3. immune_cell_percentage_distribution.csv - Wide format distribution\n")
cat("4. top_immune_cells_per_group.csv - Top 5 cells per group\n")
cat("5. group_summary_statistics.csv - Group-level statistics\n")
cat("6. group_immune_cell_composition.png - Stacked bar plot\n")
cat("7. immune_cell_percentage_heatmap.png - Heatmap visualization\n")
cat("8. immune_cell_analysis_report.txt - This report\n")

sink()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Report saved to: immune_cell_analysis_report.txt\n")
cat("All percentage data exported to CSV files.\n")

