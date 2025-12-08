# Load libraries
library(IOBR)
library(ggplot2)
library(corrplot)
library(data.table)
library(reshape2)
library(forcats)

# Set working directory (optional)
setwd(".")

# === Load Your Dataset ===

df <- fread("C:/Users/ASUS/Downloads/df_merged.csv")
df <- as.data.frame(df)

# Add sample ID as a column
df$ID <- paste0("Sample_", seq_len(nrow(df)))

# Save target and ID separately, then extract expression matrix
group_info <- df[, c("ID", "target")]
expr_mat <- df[, !(names(df) %in% c("ID", "target"))]

# Move ID as rownames
rownames(expr_mat) <- df$ID

# === Immune Cell Estimation with CIBERSORT ===
cibersort <- deconvo_tme(eset = t(expr_mat), 
                         method = "cibersort", 
                         arrays = TRUE, 
                         perm = 1000)

# === Add back group info ===
cibersort <- merge(group_info, cibersort, by = "ID")
cibersort$target <- gsub(0, "Normal", cibersort$target)
cibersort$target <- gsub(1, "CRC", cibersort$target)
cibersort$target <- as.factor(cibersort$target)
colnames(cibersort) <- gsub("_CIBERSORT", "", colnames(cibersort))

# === Violin Plot of Immune Cell Fractions ===
cibersort_long <- melt(cibersort, id.vars = "target", 
                       variable.name = "Cell", 
                       value.name = "Fraction")
non_cell_cols <- c("ID", "P-value", "Correlation", "RMSE")
cibersort_long <- cibersort_long[!cibersort_long$Cell %in% non_cell_cols, ]


cibersort_long$Fraction <- as.numeric(cibersort_long$Fraction)

library(ggplot2)

# Create stars based on significance level
annotations$stars <- cut(annotations$p_value,
                         breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
                         labels = c("***", "**", "*", "ns"))

# Base plot
p <- ggplot(cibersort_long, aes(x = Cell, y = Fraction, fill = target)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = c("Normal" = "skyblue", "CRC" = "tomato")) +
  theme_bw() +
  labs(x = "", y = "Estimated Fraction", fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12))

# Add p-values and stars
p <- p + 
  geom_text(data = annotations,
            aes(x = Cell, y = y,
                label = ifelse(p_value < 0.001, "p < 0.001", paste0("p = ", round(p_value, 3)))),
            size = 3.2, inherit.aes = FALSE, vjust = -0.2) +
  geom_text(data = annotations,
            aes(x = Cell, y = y - 0.02, label = stars),
            size = 4, inherit.aes = FALSE)

# Display plot
print(p)

# Save plot
ggsave("Results/Analysis/violin_plot_with_pvalues.png", plot = p, width = 12, height = 7, dpi = 300)

# Load required package
library(corrplot)

# Subset only immune cell fractions (exclude metadata)
immune_cells <- cibersort[, !(colnames(cibersort) %in% c("ID", "target", "P-value", "Correlation", "RMSE"))]

# Convert all to numeric
immune_cells <- as.data.frame(sapply(immune_cells, as.numeric))
rownames(immune_cells) <- cibersort$ID

# Pearson correlation
cor_matrix <- cor(immune_cells, method = "pearson", use = "pairwise.complete.obs")

# Color palette
color_palette <- colorRampPalette(c("blue", "white", "tomato"))(200)

# Create output directory if missing
if (!dir.exists("Results/Analysis")) dir.create("Results/immuneAnalysis", recursive = TRUE)

# Save heatmap
png("Results/Analysis/immune_cell_correlation_heatmap.png", width = 3000, height = 3000, res = 300)

corrplot(cor_matrix,
         method = "square",        # square-shaped cells
         col = color_palette,
         type = "full",
         order = "hclust",         # hierarchical clustering
         addgrid.col = "grey85",   # subtle box spacing
         tl.col = "black",         # label text color
         tl.srt = 45,              # rotate text labels
         tl.cex = 0.9,             # smaller text size
         number.cex = 0.65,        # smaller correlation numbers
         addCoef.col = "black",    # color of correlation values
         cl.ratio = 0.15,          # legend width
         cl.cex = 1.1,             # legend text size
         mar = c(1, 1, 3, 1),      # outer margins
         diag = TRUE               # keep diagonal
)

dev.off()





# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# === Prepare Data ===
# Expression data (6 genes of interest)
genes <- c("ADH1B", "CDH3", "GUCA2B", "PDE9A", "ABCA8", "GUCA2A")

# Extract gene expression and immune cells
gene_expr <- df[, genes]
immune_cells <- cibersort[, !(names(cibersort) %in% c("ID", "target", "P-value", "Correlation", "RMSE"))]

# Make sure rows align by sample
gene_expr <- gene_expr[match(cibersort$ID, df$ID), ]
rownames(gene_expr) <- cibersort$ID
rownames(immune_cells) <- cibersort$ID

# === Compute Correlation Matrix Between Genes and Immune Cells ===
results <- data.frame()

for (g in genes) {
  for (ic in colnames(immune_cells)) {
    test <- cor.test(gene_expr[[g]], immune_cells[[ic]], method = "pearson")
    results <- rbind(results, data.frame(
      Gene = g,
      Immune_Cell = ic,
      Correlation = test$estimate,
      P_value = test$p.value
    ))
  }
}

# === Prepare for Plot ===
results <- results %>%
  mutate(
    Signif = ifelse(P_value < 0.05, "<0.05", round(P_value, 3)),
    DotSize = abs(Correlation),
    PvalClass = cut(P_value, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                    labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")),
    TextColor = ifelse(P_value < 0.05, "red", "black")
  )

# === Plot ===
ggplot(results, aes(x = Correlation, y = fct_reorder(Immune_Cell, Correlation))) +
  geom_point(aes(size = DotSize, color = PvalClass)) +
  geom_text(aes(label = Signif), hjust = -0.3, color = results$TextColor, size = 3) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("0-0.2" = "red", "0.2-0.4" = "#FDAE61", "0.4-0.6" = "#FEE08B",
                                "0.6-0.8" = "#D9EF8B", "0.8-1" = "grey70")) +
  scale_size(range = c(2, 6)) +
  labs(title = "Correlation between Diagnostic Genes and Immune Cells",
       x = "Correlation Coefficient", y = "Immune Cell Type", 
       color = "P-value", size = "|Correlation|") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 9))
ggsave("Results/Analysis/gene_immune_correlation_dotplot.png", width = 12, height = 8, dpi = 300)






# Load libraries
library(ggplot2)
library(dplyr)
library(forcats)

# === Define genes of interest ===
genes <- c("ADH1B", "CDH3", "GUCA2B", "PDE9A", "ABCA8", "GUCA2A")

# === Extract expression + immune data ===
# 'df' = gene expression dataframe (must include 'ID' and gene columns)
# 'cibersort' = immune cell fraction dataframe (must include 'ID' and immune columns)

# Subset gene expression
gene_expr <- df[, c("ID", genes)]

# Subset immune cells, removing metadata columns
immune_cells <- cibersort[, !(names(cibersort) %in% c("target", "P-value", "Correlation", "RMSE"))]

# Match samples by ID
matched_indices <- match(cibersort$ID, gene_expr$ID)
if (any(is.na(matched_indices))) stop("Some sample IDs in cibersort do not match gene_expr.")
gene_expr <- gene_expr[matched_indices, ]
rownames(gene_expr) <- gene_expr$ID
rownames(immune_cells) <- cibersort$ID

# Remove ID column to keep numeric data only
gene_expr <- gene_expr[, genes]
immune_cells <- immune_cells[, !(names(immune_cells) %in% c("ID"))]

# Create output directory
dir.create("Results/Analysis/gene_individual_plots", recursive = TRUE, showWarnings = FALSE)

# Custom color scale for P-value bins
pval_colors <- c(
  "0-0.2" = "#d73027",
  "0.2-0.4" = "#fc8d59",
  "0.4-0.6" = "#fee08b",
  "0.6-0.8" = "#d9ef8b",
  "0.8-1" = "#91cf60"
)

# === Loop over each gene ===
for (gene in genes) {
  res <- data.frame()
  
  for (ic in colnames(immune_cells)) {
    test <- cor.test(gene_expr[[gene]], immune_cells[[ic]], method = "pearson")
    res <- rbind(res, data.frame(
      Gene = gene,
      Immune_Cell = ic,
      Correlation = test$estimate,
      P_value = test$p.value
    ))
  }
  
  # Annotate results for plotting
  res <- res %>%
    mutate(
      P_label = ifelse(P_value < 0.001, "<0.001", sprintf("%.3f", P_value)),
      TextColor = ifelse(P_value < 0.05, "red", "black"),
      Size = abs(Correlation),
      Correlation_Label = sprintf("%.2f", Correlation),
      Hjust_Label = ifelse(Correlation >= 0, -0.1, 1.1),
      PvalueClass = cut(P_value, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1"))
    )
  
  # Plot: Lollipop style
  p <- ggplot(res, aes(x = Correlation, y = fct_reorder(Immune_Cell, Correlation))) +
    geom_segment(aes(x = 0, xend = Correlation, yend = Immune_Cell),
                 color = "black", size = 1) +
    geom_point(aes(size = Size, color = PvalueClass)) +
    geom_text(
      aes(label = Correlation_Label, x = Correlation + ifelse(Correlation >= 0, 0.05, -0.05)),
      hjust = ifelse(res$Correlation >= 0, 0, 1),
      color = "black", size = 3.5
    ) +
    scale_color_manual(values = pval_colors) +
    scale_size_continuous(range = c(3, 8)) +
    coord_cartesian(xlim = c(min(res$Correlation) - 0.3, max(res$Correlation) + 0.3)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = "#91cf60"),
      panel.background = element_rect(fill = "white", color = "#91cf60"),
      strip.background = element_rect(fill = "white", color = "#91cf60"),
      legend.background = element_rect(fill = "white", color = "#91cf60"),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    ) +
    labs(
      title = paste("", gene),
      x = "Correlation Coefficient",
      y = NULL,
      color = "P value",
      size = "abs(Correlation)"
    )
  
  # Save plot
  ggsave(paste0("Results/AAAnalysis/gene_individual_plots/", gene, "_correlation.png"),
         plot = p, width = 7, height = 7, dpi = 800, bg = "white")
}
dev.off()


  










