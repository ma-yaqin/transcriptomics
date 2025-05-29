#!/bin/bash
# PUBLIC DATA FETCH
# This script downloads SRA files and splits them into FASTQ files for paired-end sequencing

# PARAMS
METADATA="raw/metadata.tsv"

# Create a directory for fastq files if it doesn't exist
mkdir -p raw/fastq

# SRA download and conversion
while IFS=$'\t' read -r sra_id _ || [[ -n "$sra_id" ]]; do
    # Skip empty lines or lines without SRA IDs
    if [[ -z "$sra_id" ]]; then
        continue
    fi

    # Define the FASTQ file names
    fastq1="raw/${sra_id}_1.fastq.gz"
    fastq2="raw/${sra_id}_2.fastq.gz"

    # Check if the FASTQ files already exist
    if [[ -f "$fastq1" && -f "$fastq2" ]]; then
        echo "FASTQ files already exist for $sra_id. Skipping download and conversion."
        continue
    fi

    echo "Fetching data for SRA ID: $sra_id"
    
    # Run prefetch to download SRA file
    if prefetch -o raw/"$sra_id" "$sra_id"; then
        echo "Successfully fetched: $sra_id"
        
        # Split into paired FASTQ files using fastq-dump
        echo "Converting SRA to FASTQ for $sra_id"
        if fastq-dump --split-files --gzip -O raw/ raw/"$sra_id"; then
            echo "FASTQ files created for: $sra_id"
        # Remove the original SRA file after successful conversion
            rm -f raw/"$sra_id"
        else
            echo "Error converting SRA to FASTQ for: $sra_id" >&2
        fi
    else
        echo "Error fetching: $sra_id" >&2
    fi
done < "$METADATA"
