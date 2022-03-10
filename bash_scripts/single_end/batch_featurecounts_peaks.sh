# to be run on /mnt/d/PERSONAL_SACO/Cubas lab/BRC1 ChIPseq/sam_output
# to generate .gff check R script
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

GTF_FILE=/mnt/d/PERSONAL_SACO/Cubas\ lab/R/output/macs2_consensus.gff
BAM_FILES=(*_sorted.bam)
TODAY=$(date +"%Y_%m_%d")

date
echo "Processing ${BAM_FILES} to ${TODAY}_feature_counts_macs2_consensus.txt"
time featureCounts -a "${GTF_FILE}" -o ${TODAY}_feature_counts_macs2_consensus.txt -t gene -g gene_id -O -T 4 ${BAM_FILES[@]}
