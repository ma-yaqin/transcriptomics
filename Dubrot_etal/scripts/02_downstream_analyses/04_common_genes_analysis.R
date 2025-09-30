### ANALYSIS ON COMMON GENES ####
# Date: May 31st, 2025
# Author: MA Yaqin
# Description: conduct the analysis on genes with similar upregulation patterns
# in various cell lines

#======== Libraries 
library(tidyverse)
library(ggplot2)
library(paletteer)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(circlize)
library(grid)

#======== Parameters 
PROCDIR <- paste0(getwd(),"/processed")
METADIR <- paste0(getwd(), "/metadata")
RESDIR <- paste0(getwd(),"/result")

#======== Data import 
bg_list <- readRDS(file.path(RESDIR, "common_bg_genes.rds"))
upDE_list <- readRDS(file.path(RESDIR, "common_upDE_genes.rds"))
downDE_list <- readRDS(file.path(RESDIR, "common_downDE_genes.rds"))

markers_ifng <- c("Stat1","Irf1","Gbp2","Gbp5","Gbp7","Cxcl9","Cxcl10","Cxcl11",
                  "Tap1","Tap2","Tapbp","B2m","Psmb8","Psmb9","Psmb10","Cd274",
                  "Ifngr1","Ifngr2","Socs1","Socs3","Jak1","Jak2","Isg15","Icam1",
                  "Nos2","Irf2","Irf4","Irf8","Irf9","Gbp3","Cybb","Prkcd",
                  "Ptpn11","Ciita","Spi1","Igtp","Ifi47","Ifi202")


#======== Heatmap of upregulated genes
# Combine all background DE dataframes
bg_all <- bind_rows(lapply(names(bg_list), function(cell) {
  bg_list[[cell]] %>%
    mutate(cell_line = cell) %>%
    dplyr::select(symbol, log2FoldChange, padj, cell_line)
}))

# Reshape to wide format: log2FC matrix
lfc_mat <- bg_all %>%
  dplyr::select(symbol, cell_line, log2FoldChange) %>%
  pivot_wider(names_from = cell_line, 
              values_from = log2FoldChange) %>%
  as.data.frame() %>%
  column_to_rownames("symbol") %>%
  replace(is.na(.),0)

# Keep only marker genes
marker_lfc_mat <- lfc_mat[rownames(lfc_mat) %in% markers_ifng, , drop = FALSE] %>%
  as.matrix()

fc_min <- min(lfc_mat, na.rm = TRUE)
fc_max <- max(lfc_mat, na.rm = TRUE)

log2fc_col_fun <- colorRamp2(c(fc_min, 0, fc_max), 
                             c("blue", "white", "red"))

# FDR matrix (subset to markers too)
fdr_mat <- bg_all %>%
  dplyr::select(symbol, cell_line, padj) %>%
  pivot_wider(names_from = cell_line, values_from = padj) %>%
  as.data.frame() %>%
  column_to_rownames("symbol")

marker_fdr_mat <- fdr_mat[rownames(fdr_mat) %in% markers_ifng, , drop = FALSE] %>%
  as.matrix()

# Heatmap with FDR overlay (multiple significance levels)
ht <- Heatmap(
  marker_lfc_mat,
  name = "log2FC",
  col = log2fc_col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),  # emphasize markers
  column_title = "IFNγ Marker Genes",
  cell_fun = function(j, i, x, y, w, h, fill) {
    fdr <- marker_fdr_mat[i, j]
    if (!is.na(fdr)) {
      if (fdr < 0.001) {
        grid.points(x, y, pch = 16, size = unit(2, "mm"), gp = gpar(col = "black"))  # circle
      } else if (fdr < 0.05) {
        grid.points(x, y, pch = 17, size = unit(2, "mm"), gp = gpar(col = "black"))  # triangle
      } else if (fdr < 0.1) {
        grid.points(x, y, pch = 15, size = unit(2, "mm"), gp = gpar(col = "black"))  # square
      }
    }
  }
)

#  Add legend for FDR significance
fdr_legend <- Legend(labels = c("FDR < 0.001", "FDR < 0.05", "FDR < 0.1"),
                     type = "points",
                     pch = c(16, 17, 15),
                     legend_gp = gpar(col = "black"),
                     title = "Significance")

# Draw heatmap with custom legend
heatmap <- draw(ht, annotation_legend_list = list(fdr_legend))


#======== Cluster comparison
# Background gene list
bg_df_list <- lapply(bg_list, function(x) x[, c("entrezid", "symbol", "log2FoldChange", "padj")])
bg_gene_list <- lapply(bg_list, function(x) x$entrezid)

# Upregulated gene list
upDE_df_list <- lapply(upDE_list, function(x) x[, c("entrezid", "symbol", "log2FoldChange", "padj")])
upDE_gene_list <- lapply(upDE_list, function(x) x$entrezid)

#==== GO BP
# ORA
compGOBP <- compareCluster(geneCluster = upDE_gene_list,
                           universe = bg_gene_list,
                           fun = enrichGO, 
                           OrgDb = org.Mm.eg.db,   
                           ont = "BP",            
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05)
up_bpora <- dotplot(compGOBP, showCategory = 10)
up_bpora

# network plot
compGOBP <- setReadable(compGOBP, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(compGOBP, showCategory = 10)

# emapplot
compGOBP_pairwise <- pairwise_termsim(compGOBP)
emapplot(compGOBP_pairwise)

#==== GO CC
# ORA
compGOCC <- compareCluster(geneCluster = upDE_gene_list, 
                           universe = bg_gene_list,
                           fun = enrichGO, 
                           OrgDb = org.Mm.eg.db,   
                           ont = "CC",            
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05)
up_ccora <- dotplot(compGOCC, showCategory = 10)
up_ccora

# network plot
compGOCC <- setReadable(compGOCC, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(compGOCC, showCategory = 10)

# emapplot
compGOCC_pairwise <- pairwise_termsim(compGOCC)
emapplot(compGOCC_pairwise)

#==== GO MF
# ORA
compGOMF <- compareCluster(geneCluster = upDE_gene_list, 
                           universe = bg_gene_list,
                           fun = enrichGO, 
                           OrgDb = org.Mm.eg.db,   
                           ont = "MF",            
                           pAdjustMethod = "BH", 
                           pvalueCutoff = 0.05)
up_mfora <- dotplot(compGOMF, showCategory = 10)
up_mfora

# network plot
compGOMF <- setReadable(compGOMF, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
cnetplot(compGOMF, showCategory = 10)

# emapplot
compGOMF_pairwise <- pairwise_termsim(compGOMF)
emapplot(compGOMF_pairwise)

#==== Pathway reactomePA
# Pathway reactomePA
compRPA <- compareCluster(geneCluster = upDE_gene_list,
                           universe = bg_gene_list,
                           fun = enrichPathway, 
                           organism = "mouse",
                           pvalueCutoff = 0.05)
up_rpaora <- dotplot(compRPA, showCategory = 10)
up_rpaora
for (n in c("png", "jpeg", "tiff")) {
  ggsave(plot = up_rpaora,
         filename = file.path(RESDIR, paste0("Upregulated genes Reactome pathways dotplot.", n)),
         width = 12, height = 12, units = "in",
         device = n)  
  
}


# emapplot
compRPA_pairwise <- pairwise_termsim(compRPA)
rpa_map <- emapplot(compRPA_pairwise)
rpa_map
for (n in c("png", "jpeg", "tiff")) {
  ggsave(plot = rpa_map,
         filename = file.path(RESDIR, paste0("Upregulated genes Reactome pathways map.", n)),
         width = 8, height = 8, units = "in",
         device = n)  
  
}


# common plot
p <- ggarrange(up_rpaora, rpa_map, ncol = 2, labels = c("A","B"))
for (n in c("png", "jpeg", "tiff")) {
  ggsave(plot = p,
         filename = file.path(RESDIR, paste0("Upregulated genes Reactome pathways.", n)),
         width = 12, height = 12, units = "in",
         device = n)  
  
}
