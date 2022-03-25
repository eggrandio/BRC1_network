source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

for file in *_1.fq.gz; do
    [ -f "$file" ] || continue
    echo "Processing ${file%%_*}_1.fq.gz and ${file%%_*}_2.fq.gz"
    #fastp -i "${file%%_*}_1.fq.gz" -I "${file%%_*}_2.fq.gz" -o ${file%%_*}_1_trimmed.fq.gz -O ${file%%_*}_2.fq.gz --disable_adapter_trimming --html fastp_$file.html
    #rm fastp.json
done