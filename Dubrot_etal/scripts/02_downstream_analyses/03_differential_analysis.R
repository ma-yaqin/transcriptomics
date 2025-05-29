### TRANSCRIPT DATA IMPORT AND PREPROCESSING ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: preprocess

#======== Libraries 
library(DESeq2)
library(org.Mm.eg.db)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)

#======== Parameters 
PROCDIR <- paste0(getwd(),"/processed/")
METADIR <- paste0(getwd(), "/metadata/")
RESDIR <- paste0(getwd(),"/result/")

#======== Data import 
txi <- readRDS(paste0(RESDIR, "transcript_abundance.rds"))
tx2gene <- readRDS(paste0(RESDIR, "transcript_annotation.rds"))
metadata <- readRDS(paste0(RESDIR, "metadata.rds"))

#======== Data preprocesss
dds <- DESeqDataSetFromTximport(txi,
                                colData = metadata,
                                design = ~ cell_line)
boxplot(log10(assay(dds)), 
        las = 3, 
        xlab = "Sample", 
        ylab = "Log10 of read counts", 
        main = "Read Counts Distribution")
dds <- estimateSizeFactors(dds)
vsd <- vst(dds,blind=TRUE)
boxplot(assay(vsd), 
        xlab="", 
        ylab="Log2 counts per million",
        las=3,
        main="Normalised Distributions")
abline(h=median(assay(vsd)), col="blue")

PCA <- plotPCA(vsd,intgroup=c("cell_line"), returnData = TRUE)
percentVar <- round(100 * attr(PCA, "percentVar"))
ggplot(PCA, aes(x = PC1, y=PC2, col = cell_line, pch = cell_line)) + 
  geom_point(size = 5) + 
  theme_light()+
  scale_color_brewer(palette = "Dark2") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))


#======== Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05, lfcThreshold = 1, 
               contrast = c("cell_line", "IFNg", "NS"))
summary(res)
plotMA(res, ylim = c(-5, 5), main = "MA-Plot of NS vs INFg")
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("cell_line"))

#======== Enrichment analysis
# background genes
bg_res <- res %>%
  as.data.frame() %>%
  mutate(symbol = mapIds(org.Mm.eg.db,
                         keys=rownames(.),
                         column="SYMBOL",
                         keytype="REFSEQ",
                         multiVals="first"),
         entrezid = mapIds(org.Mm.eg.db,
                           keys=rownames(.),
                           column="ENTREZID",
                           keytype="REFSEQ",
                           multiVals="first")) %>%
  drop_na()
bg_entrez <- bg_res$entrezid

# Significant genes
upDE <- res %>%
  subset(log2FoldChange > 1 & padj < 0.05) %>%
  as.data.frame() %>%
  mutate(symbol = mapIds(org.Mm.eg.db,
                         keys=rownames(.),
                         column="SYMBOL",
                         keytype="REFSEQ",
                         multiVals="first"),
         entrezid = mapIds(org.Mm.eg.db,
                           keys=rownames(.),
                           column="ENTREZID",
                           keytype="REFSEQ",
                           multiVals="first")) %>%
  drop_na()
upDE_entrez <- upDE$entrezid
downDE <- res %>%
  subset(log2FoldChange < -1 & padj < 0.05) %>%
  as.data.frame() %>%
  mutate(symbol = mapIds(org.Mm.eg.db,
                         keys=rownames(.),
                         column="SYMBOL",
                         keytype="REFSEQ",
                         multiVals="first"),
         entrezid = mapIds(org.Mm.eg.db,
                           keys=rownames(.),
                           column="ENTREZID",
                           keytype="REFSEQ",
                           multiVals="first")) %>%
  drop_na()
downDE_entrez <- downDE$entrezid
de_entrez <- list(upDE = upDE_entrez, 
                  downDE = downDE_entrez)
# GO biological process
GO_up <- compareCluster(upDE_entrez,
                           universe = bg_entrez,
                           fun = "enrichGO",
                           OrgDb = org.Mm.eg.db,
                           pAdjustMethod = "BH",
                           ont = "BP",
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05)
clusterProfiler::dotplot(GO_up, font.size=10)

GO_down <- enrichGO(downDE_entrez,
                    universe = bg_entrez,
                    OrgDb = "org.Mm.eg.db",
                    pAdjustMethod = "BH",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
clusterProfiler::dotplot(GO_down, font.size=10)
