BAM_FILES=(*_hisat2_sorted.bam)
EXPERIMENT_1=${BAM_FILES[@]:0:18}
EXPERIMENT_2=${BAM_FILES[@]:18:30}
TODAY=$(date +"%Y_%m_%d")
GTF_FILE=~/mapping/Ath_TAIR10_genome/Arabidopsis_thaliana.TAIR10.52.gtf
# In our case, sequencing is "reversely stranded" so we use option "-s 2" 

date
echo "Processing ${EXPERIMENT_1} to ${TODAY}_feature_counts_experiment_1_gene_hisat2.txt"
time featureCounts -p -a ${GTF_FILE} -o ${TODAY}_feature_counts_experiment_1_gene_hisat2.txt -t gene -g gene_id -O -s 2 -C -T 4 ${EXPERIMENT_1}

date
echo "Processing ${EXPERIMENT_2} to ${TODAY}_feature_counts_experiment_2_gene_hisat2.txt"
time featureCounts -p -a ${GTF_FILE} -o ${TODAY}_feature_counts_experiment_2_gene_hisat2.txt -t gene -g gene_id -O -s 2 -C -T 4 ${EXPERIMENT_2}

# # -O Assign reads to all their overlapping meta-features
# # -s 2 reversely stranded
# # -C Do not count read pairs that have their two ends mapping to different chromosomes or mapping to same chromosome but on different strands.

# # test if unstranded makes difference

# echo "Processing ${EXPERIMENT_1} to ${TODAY}_feature_counts_experiment_1_gene_hisat2.txt"
# time featureCounts -p -a ${GTF_FILE} -o ${TODAY}_feature_counts_experiment_1_gene_hisat2_unstranded.txt -t gene -g gene_id -O -C -T 4 ${EXPERIMENT_1}

# date
# echo "Processing ${EXPERIMENT_2} to ${TODAY}_feature_counts_experiment_2_gene_hisat2.txt"
# time featureCounts -p -a ${GTF_FILE} -o ${TODAY}_feature_counts_experiment_2_gene_hisat2_unstranded.txt -t gene -g gene_id -O -C -T 4 ${EXPERIMENT_2}

# # test if -t exon -g gene_id make a difference in counts

# date
# echo "Processing ${EXPERIMENT_1}"
# time featureCounts -p -a $GTF_FILE -o ${TODAY}_feature_counts_experiment_1_default.txt -O -s 2 -C -T 4 $EXPERIMENT_1

# date
# echo "Processing ${EXPERIMENT_2}"
# time featureCounts -p -a $GTF_FILE -o ${TODAY}_feature_counts_experiment_2_default.txt -O -s 2 -C -T 4 $EXPERIMENT_2