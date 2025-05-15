#!/bin/bash
# READS MAPPING - kallisto for all paired-end FASTQ files in raw/ (parallel)

# PARAMS
database="GRCm39"
threads=4

# Ensure output root directory exists
mkdir -p processed

# MAPPING AND QUANTIFICATION FOR ALL READ PAIRS (in parallel)
for r1 in raw/*_1.fastq.gz; do
  r2="${r1/_1.fastq.gz/_2.fastq.gz}"
  sample=$(basename "$r1" _1.fastq.gz)
  output_dir="processed/${sample}"
  mkdir -p "$output_dir"

  {
    echo "Running kallisto for $sample..."
    start=$(date +%s)

    kallisto quant -i database/$database -o "$output_dir" -t $threads "$r1" "$r2"

    end=$(date +%s)
    runtime=$((end - start))
    echo "$runtime seconds" > "${output_dir}/${sample}_runtime.txt"

    # Rename result files to include sample name
    for file in "$output_dir"/*; do
      base=$(basename "$file")
      mv "$file" "$output_dir/${sample}_${base}"
    done

    echo "Finished $sample in $runtime seconds."
  } &
done

# Wait for all background jobs to finish
wait
echo "All kallisto jobs completed."