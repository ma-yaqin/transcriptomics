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
- **Cell lines:** Three cancer cell lines (KPC.2, CT26.2, YUMMER1.7)
- **Conditions:** No stimulation (NS) vs IFNγ-stimulated cell lines

## Methods
1. **Preprocessing**: Quality trimming using `Trim Galore`  
2. **Quantification**: Transcript abundance estimation with `Kallisto`  
3. **EDA**: Performed on both raw counts and VST-normalized counts  
4. **Differential Expression**: Using `DESeq2`  
   - Significant DEGs: |log2FC| > 0.5 and FDR < 0.05  
5. **Enrichment Analysis**: Overrepresentation analysis of GO BP terms using `clusterProfiler`  

---

## 🔧 Requirements

Before running the scripts, ensure you have the following tools installed and available in your `PATH`:

- [SRA Toolkit](https://github.com/ncbi/sra-tools) (with `prefetch` and `fastq-dump`)  
- `wget`  
- `bash` (tested on v4+)  
- Trim Galore (0.6.10)  
- Kallisto (0.51.1)  
- DESeq2 (1.34.0)  
- clusterProfiler (4.2.2)  
- R (4.1.1)     

---

## 📂 Directory Structure
├── database/ # Reference transcriptome and Kallisto index </br>
├── raw/ </br>
│ ├── metadata.tsv # Metadata file containing SRA IDs </br>
│ ├── fastq/ # FASTQ output directory </br>
│ └── *.sra # (Temporary) downloaded SRA files </br>
└── scripts/ # Analysis scripts


---
## ⚡ Preparation
### 1. Database Preprocess

**Script:** [`scripts/00_preprocess/01_database_generation.sh`](scripts/00_preprocess/01_database_generation.sh)  

This script downloads the **Mus musculus transcriptome (mm10)** reference from UCSC and builds a **Kallisto index**.

**Usage:**

```bash
bash scripts/00_preprocess/01_database_generation.sh
```
Output:

Reference FASTA → database/mm10/refMrna.fa.gz
Kallisto index → database/mm10/mm10_mrna.idx

---
### 2. Public Data Fetch
**Script:** [`scripts/00_preprocess/02_public_data_fetch.sh`](scripts/00_preprocess/02_public_data_fetch.sh)

This script downloads RNA-seq datasets from SRA (listed in raw/metadata.tsv) and converts them into compressed FASTQ files.

**Metadata format** (metadata.tsv):
A simple tab-separated file where the first column contains SRA IDs (additional columns are ignored).

Example:

SRR123456 </br>
SRR123457 </br>
SRR123458


Usage:

```bash
bash scripts/00_preprocess/02_public_data_fetch.sh
```

Output:

Paired FASTQ files → raw/SRR123456_1.fastq.gz, raw/SRR123456_2.fastq.gz
Temporary .sra files are removed automatically after conversion.

---

### 3. Reads Mapping & Quantification

**Script:** [`scripts/01_reads_mapping/01_reads_mapping.sh`](scripts/01_reads_mapping/01_reads_mapping.sh)  

This script performs two main steps for each sample:  

1. **Preprocessing with Trim Galore** (adapter trimming and quality filtering).  
2. **Quantification with Kallisto** using the previously built index.  

**Workflow:**

- Input: Paired-end FASTQ files from `raw/`  
- Output:  
  - Quantification files (`processed/<sample>/<sample>_abundance.tsv`)  
  - Trimmed FASTQ files (`processed/<sample>/trimmed/`)  
  - Runtime logs and error summaries  

**Usage:**

```bash
bash scripts/01_reads_mapping/reads_mapping.sh
```

**Outputs per sample:** </br>
- **processed/sample/sample_abundance.tsv** – Transcript abundance estimates from Kallisto </br>
- **processed/sample/trimmed/** – Quality-trimmed FASTQ files </br>
- **processed/sample/sample_runtime.txt** – Runtime (in seconds) for quantification </br>
- **processed/sample/** – Renamed Kallisto output files (prefixed with sample name) </br>

**Additional outputs:** </br>
- **processed/kallisto_runtime_summary.tsv** – Summary table of runtime per sample </br>
- **kallisto_errors.log** – Log of errors encountered during processing </br>

**Notes:** </br>
- The script automatically skips samples if results already exist. </br>
- If paired FASTQ files are missing, the sample is skipped with a warning. </br>
- .sra downloads are expected to be converted beforehand using the Public Data Fetch step.

---
## Main analysis
### 1. Data import
**Script:** [`scripts/02_downstream_analyses/01_transcript_data_import.R`](scripts/02_downstream_analyses/01_transcript_data_import.R)  

**Database building**
To initiate the analysis, transcript abundance were annotated with their corresponding genes using the following script

```R
#======== Database building
# gene database import and preprocess
gtf_file <- paste0(getwd(),"/database/mm10/mm10.refGene.gtf.gz")
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
saveRDS(tx2gene, file.path(RESDIR, "tx2gene_mm10.rds"))
```

What it does:
- Defines the path to the GTF file (mm10.refGene.gtf.gz) containing transcript annotation for mouse (mm10).
- Creates a TxDb object (txdb) from the GTF, which allows querying transcripts, exons, and gene relationships.
- Extracts all transcript IDs (TXNAME) as keys.
- Uses select() to map transcript IDs (TXNAME) to gene IDs (GENEID) → this creates the tx2gene mapping.
- Saves the mapping to result/tx2gene_mm10.rds for later use in tximport.

**Metadata import**
Next, samples are linked with their corresponding information stored in metadata. 

```R
# metadata
metadata <- read_xlsx(file.path(METADIR, "dubrot_metadata.xlsx")) %>%
  mutate(treatment = factor(treatment, levels = c("NS", "IFNg")),
         replicate = factor(replicate, levels = c(1:3)))
samnames <- metadata$SeqID
saveRDS(metadata, file.path(RESDIR, "metadata.rds"))
```

**Abundance import**
Finally, transcript abundance estimates were summarized at the gene level:

```R
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
```
What it does:
- Imports Kallisto abundance files using the tximport package.
- Converts transcript-level quantifications into gene-level counts and abundance using the tx2gene mapping created earlier.
- Uses txOut = FALSE to ensure results are summarized by gene rather than transcript.
- Saves the imported data to result/transcript_abundance.rds for downstream analyses (e.g., DESeq2).

---
### 2. Exploratory data analysis
The following script performs EDA (Exploratory Data Analysis), normalization QC, and PCA/heatmap visualization of RNA-seq data. It uses DESeq2’s VST for variance stabilization and produces plots to assess the effect of normalization. This analysis ensures that the RNA-seq dataset is high-quality, normalized, and free of major technical biases, providing confidence for downstream differential expression or other transcriptomic analyses.

**Script:** [`scripts/02_downstream_analyses/02_exploratory_data_analysis.R`](scripts/02_downstream_analyses/02_exploratory_data_analysis.R) 

```R
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
```
This script performs exploratory data analysis (EDA) and quality control of RNA-seq count data:
1. DESeq2 dataset creation: Converts transcript abundance data into a DESeqDataSet for downstream analysis.
2. Raw count visualization:
  - Boxplots of log2-transformed raw counts per sample.
  - RLE (Relative Log Expression) plots to check deviations from gene medians before normalization.
3. PCA on raw counts: Identifies sample clustering patterns and potential outliers.
4. Normalization using VST:
  - Variance Stabilizing Transformation (VST) applied to stabilize variance across genes.
  - Generates post-normalization QC plots including boxplots, RLE plots, and PCA to assess normalization effects.
5. Heatmap of top variable genes: Shows expression patterns and sample clustering using the 500 most variable genes.
6. Combined visualization: Produces multi-panel figures summarizing raw vs. normalized counts and PCA results.
7. Output formats: Figures saved in PNG, JPEG, and TIFF for reporting.

**Additional per-cell line analysis**
This part of the script performs cell line-specific exploratory analyses to better understand the transcriptional variation within each cell line:
1. Subset normalized data: For each unique cell_line, a subset of VST-normalized expression data is created.
2. PCA within each cell line:
  - Principal Component Analysis (PCA) is performed using the top variable genes to examine variation and clustering by treatment.
  - Both the first two PCs and the variance contribution of all PCs are visualized.
3. Heatmap of top variable genes:
4. For each cell line, a heatmap of the 500 most variable genes is generated.
5. Samples are annotated by treatment and replicate to highlight clustering patterns.

**Output:**
Combined figures (PCA, explained variance, and heatmap) are saved for each cell line in PNG, JPEG, and TIFF formats.

```R
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
```



