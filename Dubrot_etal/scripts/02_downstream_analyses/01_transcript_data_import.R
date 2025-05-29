### TRANSCRIPT DATA IMPORT AND PREPROCESSING ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: preprocess

#======== Libraries 
library(tximport)
library(edgeR)
library(tidyverse)
library(readxl)
library(ggplot2)
library(biomaRt)
library(org.Mm.eg.db)

#======== Parameters 
PROCDIR <- paste0(getwd(),"/processed/")
METADIR <- paste0(getwd(), "/metadata/")
RESDIR <- paste0(getwd(),"/result/")

#======== Data import 
# metadata
metadata <- read_xlsx(paste0(METADIR, "dubrot_metadata.xlsx")) %>%
  mutate(cell_line = factor(cell_line, levels = c("NS", "IFNb", "IFNg")),
         replicate = factor(replicate, levels = c(1:3)))
samnames <- metadata$SeqID

# abundance file
abd_files <- file.path(paste0(PROCDIR,samnames,"/",
                              samnames,"_abundance.tsv"))
names(abd_files) <- samnames
txi <- tximport(abd_files, type = "kallisto", txOut = TRUE)
rownames(txi$abundance) <- gsub("\\.\\d+$", "", rownames(txi$abundance))
rownames(txi$counts) <- gsub("\\.\\d+$", "", rownames(txi$counts))
rownames(txi$length) <- gsub("\\.\\d+$", "", rownames(txi$length))
saveRDS(txi, paste0(RESDIR, "transcript_abundance.rds"))
# annotation
tx_ids <- rownames(txi$abundance)
tx2gene <- AnnotationDbi::select(org.Mm.eg.db,
                                 keys = tx_ids,
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))
saveRDS(tx2gene, paste0(RESDIR, "transcript_annotation.rds"))