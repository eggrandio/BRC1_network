for file in *.fastq.gz; do
    [ -f "$file" ] || continue
    echo "Processing $file"
    fastp -i $file -o trimmed_$file -A --max_len1 50 --html fastp_$file.html
    rm fastp.json
done