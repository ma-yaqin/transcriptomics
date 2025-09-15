# RNA-seq Analysis of IFNγ-stimulated Cells

This repository contains scripts to perform RNA-seq analysis for studying transcriptional changes in **IFNγ-stimulated cells**. The workflow includes:

1. **Database preparation** – fetching the reference transcriptome and building an index for quantification.  
2. **Public data retrieval** – downloading SRA datasets and converting them into FASTQ files.  
3. **Downstream analysis** – alignment-free quantification, differential expression, and biological interpretation.  

---

## 🔧 Requirements

Before running the scripts, ensure you have the following tools installed and available in your `PATH`:

- [Kallisto](https://pachterlab.github.io/kallisto/) ≥ v0.48  
- [SRA Toolkit](https://github.com/ncbi/sra-tools) (with `prefetch` and `fastq-dump`)  
- `wget`  
- `bash` (tested on v4+)  

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

**Script:** [`scripts/database_preprocess.sh`](scripts/database_preprocess.sh)  

This script downloads the **Mus musculus transcriptome (mm10)** reference from UCSC and builds a **Kallisto index**.

**Usage:**

```bash
bash scripts/database_preprocess.sh
```
Output:

Reference FASTA → database/mm10/refMrna.fa.gz
Kallisto index → database/mm10/mm10_mrna.idx

---
### 2. Public Data Fetch

**Script:** scripts/public_data_fetch.sh

This script downloads RNA-seq datasets from SRA (listed in raw/metadata.tsv) and converts them into compressed FASTQ files.

**Metadata format** (metadata.tsv):
A simple tab-separated file where the first column contains SRA IDs (additional columns are ignored).

Example:

SRR123456 </br>
SRR123457 </br>
SRR123458


Usage:

```bash
bash scripts/public_data_fetch.sh
```

Output:

Paired FASTQ files → raw/SRR123456_1.fastq.gz, raw/SRR123456_2.fastq.gz
Temporary .sra files are removed automatically after conversion.

---

### 3. Reads Mapping & Quantification

**Script:** [`scripts/reads_mapping.sh`](scripts/reads_mapping.sh)  

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
bash scripts/reads_mapping.sh
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

