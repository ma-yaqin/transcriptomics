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
├── database/ # Reference transcriptome and Kallisto index
├── raw/
│ ├── metadata.tsv # Metadata file containing SRA IDs
│ ├── fastq/ # FASTQ output directory
│ └── *.sra # (Temporary) downloaded SRA files
└── scripts/ # Analysis scripts


---
## ⚡ Preparation
### 1. Database Preprocess

**Script:** [`scripts/database_preprocess.sh`](scripts/database_preprocess.sh)  

This script downloads the **Mus musculus transcriptome (mm10)** reference from UCSC and builds a **Kallisto index**.

**Usage:**

```bash
bash scripts/database_preprocess.sh

Output:

Reference FASTA → database/mm10/refMrna.fa.gz
Kallisto index → database/mm10/mm10_mrna.idx

---
### 2. Public Data Fetch

**Script:** scripts/public_data_fetch.sh

This script downloads RNA-seq datasets from SRA (listed in raw/metadata.tsv) and converts them into compressed FASTQ files.

Metadata format (metadata.tsv):
A simple tab-separated file where the first column contains SRA IDs (additional columns are ignored).

Example:

SRR123456
SRR123457
SRR123458


Usage:

```bash
bash scripts/public_data_fetch.sh


Output:

Paired FASTQ files → raw/SRR123456_1.fastq.gz, raw/SRR123456_2.fastq.gz

Temporary .sra files are removed automatically after conversion.

