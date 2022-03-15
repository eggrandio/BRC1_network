# Install required tools through conda:
conda activate mapping
conda install bowtie2
conda install fastp
conda install fastqc
conda install hisat2
conda install macs2
conda install agat
(...)

# Retrieve Arabidopsis genome from ensemble:
cd ~/mapping/Ath_TAIR10_genome/
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.*

# Build bowtie2 and HISAT2 genome indexes
# Build bowtie2 index:
fasta_files=$(find *.fa | paste -sd,)
bowtie2-build --threads 4 -f $fasta_files athaliana

# Retrieve ensemble plants 52 annotation file (.gff3) for hisat2 index
cd ~/mapping/Ath_TAIR10_genome/
wget http://ftp.ensemblgenomes.org/pub/plants/release-52/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.52.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.52.gff3.gz

# hisat2 scripts require .gtf format
agat_convert_sp_gff2gtf.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 -o Arabidopsis_thaliana.TAIR10.52.gtf

hisat2_extract_splice_sites.py Arabidopsis_thaliana.TAIR10.52.gtf > Arabidopsis_thaliana.TAIR10.52_splice_sites
hisat2_extract_exons.py Arabidopsis_thaliana.TAIR10.52.gtf > Arabidopsis_thaliana.TAIR10.52_exons

# hisat2 index requires unzipped fasta files
gunzip -k *.fa.gz
fasta_files=$(find *.fa | paste -sd,)
hisat2-build --threads 4 --ss Arabidopsis_thaliana.TAIR10.52_splice_sites --exon Arabidopsis_thaliana.TAIR10.52_exons -f $fasta_files athaliana

# To check strandness of fastq (if not known), align one file with hisat2 without setting --rna-strandness, then use rseqc to check strandness (requires .bed file for annotation). To interpret results, check https://www.biostars.org/p/295344/
# In our case we know it is reverse stranded

agat_convert_sp_gff2bed.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 -o Arabidopsis_thaliana.TAIR10.52.bed

infer_experiment.py -r ~/mapping/Ath_TAIR10_genome/Arabidopsis_thaliana.TAIR10.52.bed -i A01_hisat2_sorted.bam

# Check intron size from .gff to set limit in hisat2

agat_sp_manage_introns.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 > Arabidopsis_thaliana.TAIR10.52_exon_intron_info.txt # to write info
agat_sp_manage_introns.pl --gff Arabidopsis_thaliana.TAIR10.52.gff3 --plot --break 200 # to visualize info

# Check quality of sequences:
fasta_files=*.fastq.gz
fastqc -o fastqc_output -t 8 ${fasta_files}

# Align using bowtie2/HISAT2

# To generate .gff file from MACS2 consensus peaks
GFF_FILE=/mnt/d/PERSONAL_SACO/Cubas\ lab/R/output/test.gff
BED_FILE=/mnt/d/PERSONAL_SACO/Cubas\ lab/R/output/test.bed

agat_convert_sp_gff2gtf.pl --gff "${GFF_FILE}" -o macs2_consensus.gtf

agat_convert_bed2gff.pl --bed "${BED_FILE}" -o macs2_consensus.gff