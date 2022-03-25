# fragment size was estimated with R csaw correlateReads(.bam) %>% maximizeCcf()
# TAIR10 effective genome size is 119481543
source ~/miniconda3/etc/profile.d/conda.sh
conda activate mapping

# # Default qvalue threshold (0.05)
# Q_VALUE=0.05
# time macs2 callpeak -t ChIP2_S5_R1_001_sorted.bam -c Input2_S2_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP2_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP3_S6_R1_001_sorted.bam -c Input3_S3_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP3_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP7_S4_R1_001_sorted.bam -c Input7_S1_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP7_q${Q_VALUE}" --nomodel --extsize 208 --call-summits

# # qvalue threshold (0.01)
# Q_VALUE=0.01
# time macs2 callpeak -t ChIP2_S5_R1_001_sorted.bam -c Input2_S2_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP2_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP3_S6_R1_001_sorted.bam -c Input3_S3_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP3_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP7_S4_R1_001_sorted.bam -c Input7_S1_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP7_q${Q_VALUE}" --nomodel --extsize 208 --call-summits

# # qvalue threshold (0.001)
# Q_VALUE=0.001
# time macs2 callpeak -t ChIP2_S5_R1_001_sorted.bam -c Input2_S2_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP2_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP3_S6_R1_001_sorted.bam -c Input3_S3_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP3_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
# time macs2 callpeak -t ChIP7_S4_R1_001_sorted.bam -c Input7_S1_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP7_q${Q_VALUE}" --nomodel --extsize 208 --call-summits

# high qvalue threshold to use with edgeR (1)
Q_VALUE=1
time macs2 callpeak -t ChIP2_S5_R1_001_sorted.bam -c Input2_S2_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP2_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
time macs2 callpeak -t ChIP3_S6_R1_001_sorted.bam -c Input3_S3_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP3_q${Q_VALUE}" --nomodel --extsize 208 --call-summits
time macs2 callpeak -t ChIP7_S4_R1_001_sorted.bam -c Input7_S1_R1_001_sorted.bam -g 119481543 --outdir macs2 -q ${Q_VALUE} -n "BRC1_ChIP7_q${Q_VALUE}" --nomodel --extsize 208 --call-summits