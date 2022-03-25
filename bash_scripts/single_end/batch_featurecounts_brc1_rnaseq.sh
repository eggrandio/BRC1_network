BAM_FILES=(*_hisat2_sorted.bam)
TODAY=$(date +"%Y_%m_%d")
GTF_FILE=~/mapping/Ath_TAIR10_genome/Arabidopsis_thaliana.TAIR10.52.gtf
# In our case, sequencing is "reversely stranded" so we use option "-s 2" 

date
echo "Processing ${BAM_FILES[@]} to ${TODAY}_feature_counts_gene_hisat2.txt"
time featureCounts -a ${GTF_FILE} -o ${TODAY}_feature_counts_gene_hisat2.txt -t gene -g gene_id -O -T 4 ${BAM_FILES[@]}

# -O Assign reads to all their overlapping meta-features
# -s 2 reversely stranded
# -C Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.

