### EXPLORATORY DATA ANALYSIS ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: explore characteristics of the data

#======== Libraries 
library(tximport)
library(tidyverse)
library(ggplot2)
library(org.Mm.eg.db)
library(edgeR)
library(RUVSeq)
library(EDASeq)

#======== Parameters 
PROCDIR <- paste0(getwd(),"/processed/")
METADIR <- paste0(getwd(), "/metadata/")
RESDIR <- paste0(getwd(),"/result/")

#======== Data import 
txi <- readRDS(paste0(RESDIR, "transcript_abundance.rds"))
tx2gene <- readRDS(paste0(RESDIR, "transcript_annotation.rds"))
metadata <- readRDS(paste0(RESDIR, "metadata.rds"))

#======== Data preprocess
counts <- txi$counts
counts <- round(counts)
x <- metadata$cell_line
set <- newSeqExpressionSet(as.matrix(counts),
                           phenoData = data.frame(x,
                                                  row.names=colnames(counts)))

colors <- RColorBrewer::brewer.pal(3, "Set2")
plotRLE(set, 
        outline=FALSE, ylim=c(-1, 1), col=colors[x], las = 3, 
        main = "Relative log expression")
plotPCA(set, col=colors[x], cex=1.2, 
        main = "PCA")
