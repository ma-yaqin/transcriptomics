### DIFFERENTIAL ABUNDANCE ANALYSIS ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: conduct the differential abundance analysis as well as 
# Functional enrichment analysis

#======== Libraries 
library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(ggrepel)
library(ggtext)
library(ggpubr)

#======== Functions
annotate_genes <- function(df, keytype = "SYMBOL") {
  df %>%
    mutate(
      symbol = rownames(.),  # keep the original symbol
      entrezid = mapIds(org.Mm.eg.db,
                        keys = rownames(.),
                        column = "ENTREZID",
                        keytype = keytype,
                        multiVals = "first")
    ) %>%
    drop_na(entrezid) # only keep mapped genes
}

#======== Parameters 
PROCDIR <- paste0(getwd(),"/processed/")
METADIR <- paste0(getwd(), "/metadata/")
RESDIR <- paste0(getwd(),"/result/")
fdr <- 0.05
lfc <- 0.5

#======== Data import 
txi <- readRDS(paste0(RESDIR, "transcript_abundance.rds"))
tx2gene <- readRDS(paste0(RESDIR, "transcript_annotation.rds"))
metadata <- readRDS(paste0(RESDIR, "metadata.rds"))
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ cell_line)

cell_list <- unique(metadata$cell_line)

# list to be populated with up- and downregulated genes in all cancer cells
upDE_list <- list()
downDE_list <- list()
bg_list <- list()
for (cell in cell_list) {
  cell_resdir <- file.path(RESDIR, "/cell_line/",cell)
  
  #======== Data subset and preprocess
  dds_sub <- dds[,dds$cell_line == cell]
  dds_sub$cell_line <- droplevels(dds_sub$cell_line)
  dds_sub$treatment <- droplevels(dds_sub$treatment)
  design(dds_sub) <- ~ treatment
  
  #======== Differential expression analysis
  dds_sub <- DESeq(dds_sub)
  res <- results(dds_sub, alpha = fdr,
                 contrast = c("treatment", "IFNg", "NS"))
  summary(res)
  
  #======== Result exploration and visualization
  # MA Plot
  res_df <- res %>%
    as.data.frame() %>%
    annotate_genes() %>%
    rownames_to_column("gene") %>%
    mutate(highlight = ifelse(!is.na(.$padj) & .$padj < fdr & abs(.$log2FoldChange) >= lfc, 
                              "Yes", "No"))
  
  write_delim(res_df, 
              file.path(cell_resdir, paste0("DE_result", cell, ".tsv")))
  
  mp <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange, color = highlight)) +
    geom_point(alpha = 0.5) +
    scale_x_log10() +
    geom_hline(yintercept = c(-lfc, lfc), linetype = "dashed", color = "grey") +
    theme_minimal() +
    labs(title = "MA Plot", y = "Log2 Fold Change", x = "Mean Expression", color = "Significance") +
    scale_color_manual(values = c("No" = "grey", "Yes" = "firebrick")) +
    geom_text_repel(data = head(res_df[order(res_df$padj), ], 10),
                    aes(label = symbol),
                    size = 3)
  mp
  
  # Volcano plot
  res_df_clean <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(threshold = case_when(.$padj < fdr & .$log2FoldChange >= lfc ~ "IFNg",
                                 .$padj < fdr & .$log2FoldChange <= -lfc ~ "NS",
                                 TRUE ~ "Not Signif."))
  res_df_clean_label <- res_df_clean %>%
    filter(symbol %in% c("H2-T23", "B2m", "Tap1", "Tap2", "Tapbp", "Erap1"))
    
  vp <- ggplot(res_df_clean, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("IFNg" = "#1B9E77", "NS" = "#D95F02", "Not Signif." = "grey")) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "black") +
    geom_text_repel(data = res_df_clean_label, aes(label = symbol), 
                    size = 3, max.overlaps = 20) +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(padj)", color = "Direction") +
    theme_minimal()
  vp
  
  # Top 10 DE genes (by padj)
  res_clean <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange), !is.na(symbol)) %>%
    distinct(symbol, .keep_all = TRUE)
  
  top_combined <- bind_rows(
    res_clean %>%
      filter(log2FoldChange >= lfc, padj <= fdr) %>%
      arrange(padj) %>%
      slice_head(n = 10),
    res_clean %>%
      filter(log2FoldChange <= -lfc, padj <= fdr) %>%
      arrange(padj) %>%
      slice_head(n = 10)
  ) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "IFNg", "NS"),
           symbol_italic = paste0("*", symbol, "*")) %>%
    mutate(symbol_md = case_when(
        direction == "IFNg" ~ paste0("<span style='color:#1b9e77'><i>", symbol, "</i></span>"),
        direction == "NS" ~ paste0("<span style='color:#d95f02'><i>", symbol, "</i></span>"))
    )
  
  dp <- ggplot(top_combined, aes(x = reorder(symbol_md, log2FoldChange), 
                                 y = log2FoldChange, 
                                 fill = direction)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(IFNg = "#1b9e77", NS = "#d95f02")) +
    labs(title = "Top Differentially Expressed Genes",
         x = "Gene", y = "Log2 Fold Change",
         fill = "Direction") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.y = element_markdown(size = 10))
  dp
  
  #======== Enrichment analysis
  # background genes
  bg_res <- res %>%
    as.data.frame %>%
    annotate_genes()
  bg_list[[cell]] <- bg_res
  bg_entrez <- bg_res$entrezid

  # upregulated genes
  upDE <- res %>%
    subset(log2FoldChange >= lfc & padj <= fdr) %>%
    as.data.frame() %>%
    annotate_genes()
  upDE_list[[cell]] <- upDE
  upDE_entrez <- upDE$entrezid
  
  
  # downregulated genes
  downDE <- res %>%
    subset(log2FoldChange <= -lfc & padj <= fdr) %>%
    as.data.frame() %>%
    annotate_genes()
  downDE_list[[cell]] <- downDE
  downDE_entrez <- downDE$entrezid
  
  # GO biological process
  GO_up <- enrichGO(upDE_entrez,
                    universe = bg_entrez,
                    OrgDb = org.Mm.eg.db,
                    pAdjustMethod = "BH",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
  write_delim(as.data.frame(GO_up), 
              file.path(cell_resdir, paste0("ORA_upregulated_genes", cell, ".tsv")))
  upora <- clusterProfiler::dotplot(GO_up, showCategory = 10, font.size=8, title = "GO BP: IFNg")
  upora
  
  GO_down <- enrichGO(downDE_entrez,
                      universe = bg_entrez,
                      OrgDb = "org.Mm.eg.db",
                      pAdjustMethod = "BH",
                      ont = "BP",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
  write_delim(as.data.frame(GO_down), 
              file.path(cell_resdir, paste0("ORA_downregulated_genes", cell, ".tsv")))
  downora <- clusterProfiler::dotplot(GO_down, showCategory = 10, font.size=8, title = "GO BP: NS")
  downora
  
  # Summarized plotting
  plot_all <- ggarrange(ggarrange(vp, dp, labels = c("A","B"),
                                  common.legend = TRUE, legend = "bottom"),
                        ggarrange(upora, downora, labels = c("C","D")),
                        nrow = 2)
  plot_all
  
  for (n in c("png", "jpeg", "tiff")) {
    ggsave(plot = plot_all,
           filename = file.path(cell_resdir, paste0("DE_result", cell, ".", n)),
           width = 12, height = 12, units = "in",
           device = n)  
    
  }
}

#======== Export the common genes
saveRDS(bg_list, file.path(RESDIR, "common_bg_genes.rds"))
saveRDS(upDE_list, file.path(RESDIR, "common_upDE_genes.rds"))
saveRDS(downDE_list, file.path(RESDIR, "common_downDE_genes.rds"))
