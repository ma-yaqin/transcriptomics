#!/bin/bash
# DATABASE PREPROCESS
# It fetches Mus musculus genome to be compared with the RNA-Seq data.

# PARAMS
database="GRCm39"
db_link="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_"$database"/GCF_000001635.27_"$database"_rna.fna.gz"
threads=4

# DATA FETCH
wget $db_link 

# INDEX BUILDING
#gzip -d $database".fa.gz"
kallisto index --index=$database --threads=$threads "GCF_000001635.27_"$database"_rna.fna.gz"