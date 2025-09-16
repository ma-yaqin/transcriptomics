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



