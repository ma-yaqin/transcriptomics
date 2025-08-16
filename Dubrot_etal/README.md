# Bulk RNA-seq Analysis of Cancer Cell Lines Under IFNγ Stimulation

## Overview
This project investigates the transcriptional dynamics of three cancer cell lines following IFNγ stimulation using publicly available RNA-seq data by (SRA accession: PRJNA850930, Dubrot et al. 2022: https://doi.org/10.1038/s41590-022-01315-x).  

The workflow includes:
- Preprocessing and quality trimming (Trim Galore)  
- Transcript quantification via pseudoalignment (Kallisto)  
- Exploratory Data Analysis (EDA) on raw and VST-normalized counts  
- Differential expression analysis (DESeq2)  
- Gene Ontology Biological Process (GO BP) enrichment (clusterProfiler)  

## Data
- **Source:** Publicly available RNA-seq data from SRA  
- **Cell lines:** Three cancer cell lines  
- **Conditions:** Control vs IFNγ stimulation  

## Methods
1. **Preprocessing**: Quality trimming using `Trim Galore`  
2. **Quantification**: Transcript abundance estimation with `Kallisto`  
3. **EDA**: Performed on both raw counts and VST-normalized counts  
4. **Differential Expression**: Using `DESeq2`  
   - Significant DEGs: |log2FC| > 0.5 and FDR < 0.05  
5. **Enrichment Analysis**: Overrepresentation analysis of GO BP terms using `clusterProfiler`  

## Tools and Versions
- Trim Galore (0.6.10)  
- Kallisto (0.51.1)  
- DESeq2 (1.34.0)  
- clusterProfiler (4.2.2)  
- R (4.1.1)   

## Repository Structure
