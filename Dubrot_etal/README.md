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

**Script:** [`scripts/00_preprocess/database_preprocess.sh`](scripts/00_preprocess/database_preprocess.sh)  

This script downloads the **Mus musculus transcriptome (mm10)** reference from UCSC and builds a **Kallisto index**.

**Usage:**

```bash
bash scripts/00_preprocess/database_preprocess.sh
```
Output:

Reference FASTA → database/mm10/refMrna.fa.gz
Kallisto index → database/mm10/mm10_mrna.idx

---
### 2. Public Data Fetch

**Script:** scripts/00_preprocess/public_data_fetch.sh

This script downloads RNA-seq datasets from SRA (listed in raw/metadata.tsv) and converts them into compressed FASTQ files.

**Metadata format** (metadata.tsv):
A simple tab-separated file where the first column contains SRA IDs (additional columns are ignored).

Example:

SRR123456 </br>
SRR123457 </br>
SRR123458


Usage:

```bash
bash scripts/00_preprocess/public_data_fetch.sh
```

Output:

Paired FASTQ files → raw/SRR123456_1.fastq.gz, raw/SRR123456_2.fastq.gz
Temporary .sra files are removed automatically after conversion.

---

### 3. Reads Mapping & Quantification

**Script:** [`scripts/01_reads_mapping/reads_mapping.sh`](scripts/01_reads_mapping/reads_mapping.sh)  

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

Outputs per sample:

processed/<sample>/<sample>_abundance.tsv – Transcript abundance estimates from Kallisto

processed/<sample>/trimmed/ – Quality-trimmed FASTQ files

processed/<sample>/<sample>_runtime.txt – Runtime (in seconds) for quantification

processed/<sample>/ – Renamed Kallisto output files (prefixed with sample name)

Additional outputs:

processed/kallisto_runtime_summary.tsv – Summary table of runtime per sample

kallisto_errors.log – Log of errors encountered during processing

Notes:

The script automatically skips samples if results already exist.

If paired FASTQ files are missing, the sample is skipped with a warning.

.sra downloads are expected to be converted beforehand using the Public Data Fetch step.


