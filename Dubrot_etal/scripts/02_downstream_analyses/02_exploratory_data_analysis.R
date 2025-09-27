### UNIFIED TRANSCRIPTOMICS QC SCRIPT ####
# Date: Aug 1, 2025
# Author: MA Yaqin
# Description: EDA, normalization QC, and PCA visualization

#======== Libraries 
library(tximport)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(org.Mm.eg.db)
library(RColorBrewer)
library(matrixStats)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(paletteer)
library(grid)

#======== Parameters 
# directories
PROCDIR <- paste0(getwd(),"/processed/")
METADIR <- paste0(getwd(), "/metadata/")
RESDIR <- paste0(getwd(),"/result/")

#======== Data import 
txi <- readRDS(paste0(RESDIR, "transcript_abundance.rds"))
tx2gene <- readRDS(paste0(RESDIR, "tx2gene_mm10.rds"))
metadata <- readRDS(paste0(RESDIR, "metadata.rds"))

#======== Prepare DESeqDataSet 
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ treatment)

#======== Plot raw count distribution (ggplot)
raw_log <- log2(assay(dds))
raw_long <- as.data.frame(raw_log) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Log10Count")

bdist <- ggplot(raw_long, aes(x = Sample, y = Log10Count)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Raw Read Counts Distribution", x = "Sample", y = "Log2 of raw read counts")
bdist

#======== Manual RLE plot BEFORE normalization (ggplot)
rle_before_mat <- sweep(assay(dds), 1, rowMedians(assay(dds), na.rm = TRUE, useNames = FALSE))
rle_before_long <- as.data.frame(rle_before_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Deviation")

rlebef <- ggplot(rle_before_long, aes(x = Sample, y = Deviation)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, coef = 2) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "RLE Before Normalization (VST)", y = "Deviation from gene median")
rlebef

#======== PCA from unnormalized data
log_counts <- log2(counts(dds, normalized = FALSE) + 0.01)
log_counts_filtered <- log_counts[apply(log_counts, 1, var) != 0, ]

pca <- prcomp(t(log_counts_filtered), scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data$treatment <- colData(dds)$treatment
pca_data$cell_line <- colData(dds)$cell_line
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

pca_before <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment, pch = cell_line)) +
  geom_point(size = 4) +
  labs(title = "PCA of raw (unnormalized) counts",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()
pca_before

#======== Run normalization
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

#======== Plot normalized expression distribution (ggplot)
vsd_long <- as.data.frame(assay(vsd)) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "VST")

adist <- ggplot(vsd_long, aes(x = Sample, y = VST)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = median(assay(vsd)), color = "blue", linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Normalized Distributions (VST)", y = "Log2 normalized expression (VST)")
adist

#======== Manual RLE plot AFTER normalization (ggplot)
rle_after_mat <- sweep(assay(vsd), 1, rowMedians(assay(vsd), na.rm = TRUE, useNames = FALSE))
rle_after_long <- as.data.frame(rle_after_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Deviation")

rleaf <- ggplot(rle_after_long, aes(x = Sample, y = Deviation)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 0.5, coef = 2) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "RLE After Normalization (VST)", y = "Deviation from gene median")
rleaf

#======== PCA from normalized data
pca_data <- plotPCA(vsd, intgroup = "treatment", returnData = TRUE)
pca_data$cell_line <- colData(dds)$cell_line
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_after <- ggplot(pca_data, aes(x = PC1, y = PC2, col = treatment, pch = cell_line)) +
  geom_point(size = 4) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "PCA of VST-normalized data",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))
pca_after

#======== Heatmap
mat <- assay(vsd)
topVarGenes <- head(order(rowVars(mat, useNames = FALSE), decreasing = TRUE), 500)
mat_top <- mat[topVarGenes, ]
mat_scaled <- t(scale(t(mat_top)))

ann_col <- as.data.frame(colData(vsd)[, c("cell_line", "treatment","replicate")])
ann_col$cell_line <- as.factor(ann_col$cell_line)
ann_col$replicate <- as.factor(ann_col$replicate)
ann_col$treatment <- as.factor(ann_col$treatment)

treatment_colors <- setNames(brewer.pal(3, "Dark2")[1:length(unique(ann_col$treatment))],
                             levels(ann_col$treatment))
cell_colors <- setNames(as.character(paletteer_d("ggsci::category10_d3")[1:length(unique(ann_col$cell_line))]),
                        levels(ann_col$cell_line))
replicate_colors <- setNames(brewer.pal(length(unique(ann_col$replicate)), "Set1"),
                             levels(ann_col$replicate))
annotation_colors <- list(treatment = treatment_colors,
                          cell_line = cell_colors,
                          replicate = replicate_colors)

# Create heatmap
h <- pheatmap(mat_scaled,
              annotation_col = ann_col,
              annotation_colors = annotation_colors,
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              show_rownames = FALSE,
              fontsize_col = 10,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
              main = "Top 500 Most Variable Genes (VST-scaled)")
h

p <- ggarrange(ggarrange(bdist, rlebef, pca_before, nrow = 3, labels = c("A", "C", "E")), 
               ggarrange(adist, rleaf, pca_after,nrow = 3, labels = c("B", "D", "F")))
p

for (i in c("png", "jpeg", "tiff")) {
  ggsave(plot = p,
         filename = paste0(RESDIR, "Normalization effect on data.",i),
         width = 8, height = 12, units = "in",
         device = i)  
  
  ggsave(plot = h,
         filename = paste0(RESDIR, "Heatmap.",i),
         width = 6, height = 8, units = "in",
         device = i)  
}

#======== PCA and HEATMAP per cell line
cell_name <- unique(metadata$cell_line)

for (cell in cell_name) {
  message("Processing cell line: ", cell)
  
  # Set output folder path
  cell_resdir <- file.path(RESDIR, "/cell_line/",cell)
  dir.create(cell_resdir, recursive = TRUE, showWarnings = FALSE)
  
  vsd_sub <- vsd[, colData(vsd)$cell_line == cell]
  
  #======== PCA from normalized data
  pca_data <- plotPCA(vsd_sub, intgroup = "treatment", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  pca_sub <- ggplot(pca_data, aes(x = PC1, y = PC2, col = treatment, pch = treatment)) +
    geom_point(size = 4) +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"))
  pca_sub
  
  #======== Explained variance
  top_n_genes <- 500
  rv <- rowVars(assay(vsd_sub), useNames = TRUE)
  select <- order(rv, decreasing = TRUE)[seq_len(top_n_genes)]
  vsd_top <- assay(vsd_sub)[select, ]
  pca_res <- prcomp(t(vsd_top), scale. = FALSE)
  percentVar_full <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
  
  # Create dataframe with all PCs
  var_df <- data.frame(PC = factor(paste0("PC", seq_along(percentVar_full)), 
                                   levels = paste0("PC", seq_along(percentVar_full))),
                       Variance = percentVar_full)
  
  # Plot all PCs
  p_var <- ggplot(var_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 2.5, check_overlap = TRUE) +
    scale_y_continuous(limits = c(0, max(var_df$Variance) + 5)) +
    labs(y = "Variance Explained (%)",
         x = "Principal Component") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  p_var
  
  #======== Heatmap
  mat <- assay(vsd_sub)
  topVarGenes <- head(order(rowVars(mat, useNames = FALSE), decreasing = TRUE), 500)
  mat_top <- mat[topVarGenes, ]
  mat_scaled <- t(scale(t(mat_top)))
  
  ann_col <- as.data.frame(colData(vsd)[, c("treatment","replicate")])
  ann_col$replicate <- as.factor(ann_col$replicate)
  ann_col$treatment <- as.factor(ann_col$treatment)
  
  # Create heatmap
  h_sub <- grid.grabExpr({
    pheatmap(mat_scaled,
             annotation_col = ann_col,
             annotation_colors = annotation_colors,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = FALSE,
             fontsize_col = 10,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
  })
  h_sub
  
  p_sub <- ggarrange(ggarrange(pca_sub, p_var, ncol = 1, labels = c("A","B")), h_sub, 
                     nrow=1, ncol=2, labels = c("","C")) %>%
    annotate_figure(top = text_grob(paste("Cell:", cell),
                                    face = "bold", size = 14, hjust = 0.5))
  p_sub
  
  for (n in c("png", "jpeg", "tiff")) {
    ggsave(plot = p_sub,
           filename = file.path(cell_resdir, paste0("Heatmap_and_PCA_cell_", cell, ".", n)),
           width = 12, height = 12, units = "in",
           device = n)  
    
  }

}
