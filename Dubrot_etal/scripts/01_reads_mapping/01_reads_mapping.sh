#!/bin/bash
# READS MAPPING
# trim_galore preprocessing + kallisto quantification (sequential)

# ===== Parameters =====
database="mm10"
db_dir="database/${database}"
index_file="${db_dir}/${database}_mrna.idx"
threads=4
log_file="kallisto_errors.log"
summary_file="processed/kallisto_runtime_summary.tsv"

# ===== Directory and logging creation =====
# Ensure output root directory exists
mkdir -p processed

# Prepare summary and error log files
echo -e "sample\truntime_seconds" > "$summary_file"
> "$log_file"  # Clear log file

# ===== Mapping =====
for r1 in raw/*_1.fastq.gz; do
  r2="${r1/_1.fastq.gz/_2.fastq.gz}"
  sample=$(basename "$r1" _1.fastq.gz)
  trim_dir="processed/${sample}/trimmed"
  output_dir="processed/${sample}"
  quant_file="${output_dir}/${sample}_abundance.tsv"

  # Check if paired file exists
  if [[ ! -f "$r2" ]]; then
    echo "[$sample] Missing paired file: $r2. Skipping..." | tee -a "$log_file"
    continue
  fi

  # Skip if already completed
  if [[ -f "$quant_file" ]]; then
    echo "[$sample] Output already exists. Skipping..."
    continue
  fi

  mkdir -p "$trim_dir"

  # Run Trim Galore if trimmed reads don't exist
  trimmed_r1="${trim_dir}/${sample}_1_val_1.fq.gz"
  trimmed_r2="${trim_dir}/${sample}_2_val_2.fq.gz"
  if [[ ! -f "$trimmed_r1" || ! -f "$trimmed_r2" ]]; then
    echo "[$sample] Running Trim Galore..."

    if ! trim_galore --paired -o "$trim_dir" "$r1" "$r2" 2>>"$log_file"; then
      echo "[$sample] Trim Galore failed. Skipping sample." | tee -a "$log_file"
      continue
    fi
  else
    echo "[$sample] Trimmed files already exist. Skipping trimming."
  fi

  # Run kallisto quant
  echo "[$sample] Running kallisto..."
  start=$(date +%s)

  if kallisto quant -i ${index_file} -o "$output_dir" -t $threads "$trimmed_r1" "$trimmed_r2" 2>>"$log_file"; then
    end=$(date +%s)
    runtime=$((end - start))
    echo -e "$sample\t$runtime" >> "$summary_file"
    echo "$runtime seconds" > "${output_dir}/${sample}_runtime.txt"

    # Rename result files to include sample name
    for file in "$output_dir"/*; do
        base=$(basename "$file")
        [[ "$base" == ${sample}_* ]] && continue
        mv "$file" "$output_dir/${sample}_${base}"
    done

    echo "[$sample] Finished in $runtime seconds."
  else
    echo "[$sample] Kallisto quant failed. See $log_file for details." | tee -a "$log_file"
  fi

done

echo "All samples processed."
echo "Summary saved to: $summary_file"
