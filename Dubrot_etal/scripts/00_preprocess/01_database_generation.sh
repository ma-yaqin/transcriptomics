#!/bin/bash
# DATABASE PREPROCESS
# Downloads Mus musculus transcriptome (mm10) from UCSC and builds Kallisto index.

set -euo pipefail

# ===== Parameters =====
database="mm10"
db_link="https://hgdownload.soe.ucsc.edu/goldenPath/${database}/bigZips/refMrna.fa.gz"
db_dir="database/${database}"
index_file="${db_dir}/${database}_mrna.idx"
fasta_file="${db_dir}/refMrna.fa.gz"
threads=4

# ===== Create directory =====
mkdir -p "${db_dir}"

# ===== Download =====
echo "[INFO] Downloading ${database} transcriptome..."
wget -c "${db_link}" -O "${fasta_file}"

# ===== Build Kallisto index =====
echo "[INFO] Building Kallisto index..."
kallisto index -i "${index_file}" -k 31 --threads=${threads} "${fasta_file}"

echo "[INFO] Index built at: ${index_file}"
