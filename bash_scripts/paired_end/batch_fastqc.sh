source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

for file in *.fq.gz; do
    [ -f "$file" ] || continue
    echo "Processing ${file}"
    date
    time fastqc --threads 8 -o fastqc_output ${file}
done