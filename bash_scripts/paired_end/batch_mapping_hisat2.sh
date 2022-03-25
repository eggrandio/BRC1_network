# Run this to log output to file:
# bash batch_hisat2_sam_bam_bigwig.sh 2>&1 | tee batch_hisat2_mapping.log
# --max-intronlen comes from "agat_sp_manage_introns.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3", where max annotated intron length was 11602 bp
# to check for strandness, align without the option and then run infer_experiment.py -r ~/mapping/Ath_TAIR10_genome/Arabidopsis_thaliana.TAIR10.52.bed -i XXX_sorted.bam
# dUTP protocols generate reverse stranded reads, so we have to use option RF (see https://www.biostars.org/p/295344/ for more details)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

for file in *_1.fq.gz; do
    [ -f "$file" ] || continue
    FQ_1="${file%%_*}_1.fq.gz"
    FQ_2="${file%%_*}_2.fq.gz"
    SAM_FILENAME="${file%%_*}_hisat2.sam"
    BAM_FILENAME="${file%%_*}_hisat2.bam"
    SORTED_BAM_FILENAME="${file%%_*}_hisat2_sorted.bam"
    BIGWIG_FILENAME="${file%%_*}_hisat2.bw"

    echo "Processing ${FQ_1} and ${FQ_2} to ${SAM_FILENAME}"
    date
    time hisat2 -p 4 --very-sensitive --rna-strandness RF --max-intronlen 11700 -t -x ~/mapping/Ath_TAIR10_genome/athaliana -1 "${FQ_1}" -2 "${FQ_2}" -S "${SAM_FILENAME}"
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
