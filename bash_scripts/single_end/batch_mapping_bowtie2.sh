# Run this to log output to file:
# bash batch_bowtie2_sam_bam_bigwig.sh 2>&1 | tee batch_mapping.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

for file in *.fastq.gz; do
    [ -f "$file" ] || continue
    FQ="${file%%.*}.fastq.gz"
    SAM_FILENAME="${file%%.*}_hisat2.sam"
    BAM_FILENAME="${file%%.*}_hisat2.bam"
    SORTED_BAM_FILENAME="${file%%.*}_hisat2_sorted.bam"
    BIGWIG_FILENAME="${file%%.*}_hisat2.bw"

    echo "Processing ${FQ} to ${SAM_FILENAME}"
    date
    time bowtie2 -p 4 --very-sensitive-local -t -x ~/mapping/Ath_TAIR10_genome/athaliana -U "${FQ}" -S "${SAM_FILENAME}"
    echo "Converting ${SAM_FILENAME} to ${BAM_FILENAME}"
    date
    time samtools view -b -@ 4 "${SAM_FILENAME}" > "${BAM_FILENAME}"
    [ -f "${BAM_FILENAME}" ] && rm ${SAM_FILENAME}
    echo "Sorting ${BAM_FILENAME} to ${SORTED_BAM_FILENAME}"
    date
    time samtools sort -@ 4 "${BAM_FILENAME}" -o "${SORTED_BAM_FILENAME}"
    [ -f "${SORTED_BAM_FILENAME}" ] && rm ${BAM_FILENAME}
    echo "Indexing ${SORTED_BAM_FILENAME}"
    date
    time samtools index -@ 4 "${SORTED_BAM_FILENAME}"
    echo "Converting ${SORTED_BAM_FILENAME} to ${BIGWIG_FILENAME}"
    date
    time bamCoverage --effectiveGenomeSize 119481543 -p 4 --binSize 5 -b "${SORTED_BAM_FILENAME}" -o "${BIGWIG_FILENAME}"
    done
