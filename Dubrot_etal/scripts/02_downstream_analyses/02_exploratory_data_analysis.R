### TRANSCRIPT DATA IMPORT AND PREPROCESSING ####
# Date: May 29th, 2025
# Author: MA Yaqin
# Description: preprocess

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

#======== Data preprocess
counts <- txi$counts
counts <- round(counts)
genes <- rownames(counts)[grep("^ENS", rownames(counts))]

