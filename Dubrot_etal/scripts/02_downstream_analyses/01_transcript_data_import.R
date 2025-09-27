### TRANSCRIPT DATA IMPORT AND PREPROCESSING ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: preprocess

#======== Libraries 
library(tximport)
library(GenomicFeatures)
library(tidyverse)
library(readxl)
library(ggplot2)
library(biomaRt)
library(org.Mm.eg.db)

#======== Parameters 
# directories
PROCDIR <- file.path(getwd(),"processed")
METADIR <- file.path(getwd(), "metadata")
RESDIR <- file.path(getwd(),"result")

#======== Database building
# gene database import and preprocess
gtf_file <- paste0(getwd(),"/database/mm10/mm10.refGene.gtf.gz")
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
saveRDS(tx2gene, file.path(RESDIR, "tx2gene_mm10.rds"))

#======== Data import 
# metadata
metadata <- read_xlsx(file.path(METADIR, "dubrot_metadata.xlsx")) %>%
  mutate(treatment = factor(treatment, levels = c("NS", "IFNg")),
         replicate = factor(replicate, levels = c(1:3)))
samnames <- metadata$SeqID
saveRDS(metadata, file.path(RESDIR, "metadata.rds"))

# abundance file
abd_files <- file.path(PROCDIR,samnames,
                    paste0(samnames,"_abundance.tsv"))
names(abd_files) <- samnames
missing_files <- abd_files[!file.exists(abd_files)]
if (length(missing_files) > 0) {
  stop(paste("Missing abundance files for these samples:\n",
             paste(names(missing_files), "→", missing_files, collapse = "\n")))
}

txi <- tximport(abd_files, type = "kallisto", tx2gene = tx2gene, txOut = FALSE)
saveRDS(txi, file.path(RESDIR, "transcript_abundance.rds"))
