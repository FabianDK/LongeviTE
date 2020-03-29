#Carnes2015 RNA-seq analysis 
mkdir -p Carnes/RNAseq/raw_reads

#1) Download data and required programs
#We obtained RNAseq fastq files directly from the authors (Wen Huang, Trudy Mackay), please see active download links in the extra file "Carnes2015_RNAseq_links.txt". Save fastq files into: Carnes/RNAseq/raw_reads
#Download STAR mapping tool: https://github.com/alexdobin/STAR/releases 
#We used: STAR-2.7.0e
#Download subread package (which includes featureCounts): http://subread.sourceforge.net/
#We used: version 1.6.3
#Download IGV: https://software.broadinstitute.org/software/igv/download
#We used: IGVSource_2.4.10

#3) FastQC check
sh required_scripts/fastq_check.sh Carnes/RNAseq/raw_reads

#4) Trim reads and remove adapters
mkdir trim
for file in *fastq.gz
do
cutadapt -j 20 -a AGATCGGAAGAGC --minimum-length 75 -q 20 -o $file.trim $file
done
mv *.trim Carnes/RNAseq/raw_reads/trim/

#3) FastQC check to see if trimming worked
sh required_scripts/fastq_check_decompr.sh Carnes/RNAseq/raw_reads/trim/

#4) rename files
cd trim
for file in *.fastq.gz.trim; do
    mv "$file" "$(basename "$file" .fastq.gz.trim).trim.fastq"
done

###################################################################################################################
## MAPPING AGAINST GENOME, BACTERIA AND DEVIATE LIBRARY
###################################################################################################################

Carnes/RNAseq/raw_reads/trim/

#5) Generate genome index
cd .. 
mkdir reference/genome_dir_add
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir reference/genome_dir_add --genomeFastaFiles reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta --sjdbGTFfile reference/dmel-all-r6.27_TE_BACT_ADD.gtf --sjdbOverhang 74

#6) Map using STAR
cd Carnes/RNAseq/raw_reads/trim
for file in *.trim.fastq
do
STAR --runThreadN 10 --genomeDir reference/genome_dir_add --readFilesIn $file --alignIntronMin 23 --alignIntronMax 268107 -–outFilterMultimapNmax 10 --outReadsUnmapped FastX --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${file%.trim.fastq}
done

#7) Index bam files
for file in *.bam
do
samtools index $file
done

#8) Load bam files and gtf into IGV: e.g., check if most reads map to exons

#9) Count reads using featureCounts with all bams as input
featureCounts -a reference/dmel-all-r6.27_TE_BACT_ADD.gtf -o readCounts_carnes.txt -t exon -g gene_id --extraAttributes gene_symbol -T 12 B1_W1_F1_TAGCTT_L002_c5ubdanxxAligned.sortedByCoord.out.bam B1_W1_F2_TCGAAG_L002_c5ubdanxxAligned.sortedByCoord.out.bam B1_W1_M1_TACAGC_L003_c5ubdanxxAligned.sortedByCoord.out.bam B1_W1_M2_CGGAAT_L005_c5ubdanxxAligned.sortedByCoord.out.bam B1_W5_F1_TAGCTT_L004_c5ubdanxxAligned.sortedByCoord.out.bam B1_W5_F2_TCGAAG_L006_c5ubdanxxAligned.sortedByCoord.out.bam B1_W5_M1_CGTACG_L004_c5ubdanxxAligned.sortedByCoord.out.bam B1_W5_M2_CTCAGA_L008_c5ubdanxxAligned.sortedByCoord.out.bam B2_W1_F1_GGCTAC_L002_c5ubdanxxAligned.sortedByCoord.out.bam B2_W1_F2_CGGAAT_L003_c5ubdanxxAligned.sortedByCoord.out.bam B2_W1_M1_TATAAT_L003_c5ubdanxxAligned.sortedByCoord.out.bam B2_W1_M2_CTAGCT_L005_c5ubdanxxAligned.sortedByCoord.out.bam B2_W5_F1_GGCTAC_L004_c5ubdanxxAligned.sortedByCoord.out.bam B2_W5_F2_CGGAAT_L007_c5ubdanxxAligned.sortedByCoord.out.bam B2_W5_M1_GAGTGG_L004_c5ubdanxxAligned.sortedByCoord.out.bam B2_W5_M2_GCGCTA_L008_c5ubdanxxAligned.sortedByCoord.out.bam B3_W1_F1_TCGGCA_L002_c5ubdanxxAligned.sortedByCoord.out.bam B3_W1_F2_TAATCG_L005_c5ubdanxxAligned.sortedByCoord.out.bam B3_W1_M1_GTGGCC_L003_c5ubdanxxAligned.sortedByCoord.out.bam B3_W1_M2_CTATAC_L005_c5ubdanxxAligned.sortedByCoord.out.bam B3_W5_F1_TCGGCA_L004_c5ubdanxxAligned.sortedByCoord.out.bam B3_W5_F2_CTAGCT_L007_c5ubdanxxAligned.sortedByCoord.out.bam B3_W5_M1_GGTAGC_L004_c5ubdanxxAligned.sortedByCoord.out.bam B3_W5_M2_TAATCG_L008_c5ubdanxxAligned.sortedByCoord.out.bam B4_W1_F1_GTAGAG_L002_c5ubdanxxAligned.sortedByCoord.out.bam B4_W1_F2_CTAGCT_L003_c5ubdanxxAligned.sortedByCoord.out.bam B4_W1_M1_CGTACG_L003_c5ubdanxxAligned.sortedByCoord.out.bam B4_W1_M2_CTCAGA_L005_c5ubdanxxAligned.sortedByCoord.out.bam B4_W5_F1_GTAGAG_L006_c5ubdanxxAligned.sortedByCoord.out.bam B4_W5_F2_CTATAC_L007_c5ubdanxxAligned.sortedByCoord.out.bam B4_W5_M1_ACTGAT_L004_c5ubdanxxAligned.sortedByCoord.out.bam B4_W5_M2_TACAGC_L008_c5ubdanxxAligned.sortedByCoord.out.bam B5_W1_F1_GTCCGC_L002_c5ubdanxxAligned.sortedByCoord.out.bam B5_W1_F2_CTATAC_L003_c5ubdanxxAligned.sortedByCoord.out.bam B5_W1_M1_GAGTGG_L003_c5ubdanxxAligned.sortedByCoord.out.bam B5_W1_M2_GCGCTA_L005_c5ubdanxxAligned.sortedByCoord.out.bam B5_W5_F1_GTCCGC_L006_c5ubdanxxAligned.sortedByCoord.out.bam B5_W5_F2_CTCAGA_L007_c5ubdanxxAligned.sortedByCoord.out.bam B5_W5_M1_CCAACA_L004_c5ubdanxxAligned.sortedByCoord.out.bam B5_W5_M2_TATAAT_L008_c5ubdanxxAligned.sortedByCoord.out.bam O1_W1_F1_CGATGT_L002_c5ubdanxxAligned.sortedByCoord.out.bam O1_W1_F2_GTTTCG_L002_c5ubdanxxAligned.sortedByCoord.out.bam O1_W1_M1_TACAGC_L005_c5ubdanxxAligned.sortedByCoord.out.bam O1_W1_M2_GGTAGC_L003_c5ubdanxxAligned.sortedByCoord.out.bam O1_W5_F1_CGATGT_L004_c5ubdanxxAligned.sortedByCoord.out.bam O1_W5_F2_GTTTCG_L006_c5ubdanxxAligned.sortedByCoord.out.bam O1_W5_M1_GCGCTA_L007_c5ubdanxxAligned.sortedByCoord.out.bam O1_W5_M2_TCATTC_L004_c5ubdanxxAligned.sortedByCoord.out.bam O2_W1_F1_TGACCA_L002_c5ubdanxxAligned.sortedByCoord.out.bam O2_W1_F2_CAACTA_L002_c5ubdanxxAligned.sortedByCoord.out.bam O2_W1_M1_TATAAT_L005_c5ubdanxxAligned.sortedByCoord.out.bam O2_W1_M2_ACTGAT_L003_c5ubdanxxAligned.sortedByCoord.out.bam O2_W5_F1_TGACCA_L004_c5ubdanxxAligned.sortedByCoord.out.bam O2_W5_F2_CAACTA_L006_c5ubdanxxAligned.sortedByCoord.out.bam O2_W5_M1_TAATCG_L007_c5ubdanxxAligned.sortedByCoord.out.bam O2_W5_M2_TCCCGA_L004_c5ubdanxxAligned.sortedByCoord.out.bam O3_W1_F1_ACAGTG_L002_c5ubdanxxAligned.sortedByCoord.out.bam O3_W1_F2_CAAAAG_L002_c5ubdanxxAligned.sortedByCoord.out.bam O3_W1_M1_CTCAGA_L003_c5ubdanxxAligned.sortedByCoord.out.bam O3_W1_M2_CCAACA_L003_c5ubdanxxAligned.sortedByCoord.out.bam O3_W5_F1_ACAGTG_L004_c5ubdanxxAligned.sortedByCoord.out.bam O3_W5_F2_CAAAAG_L006_c5ubdanxxAligned.sortedByCoord.out.bam O3_W5_M1_TACAGC_L007_c5ubdanxxAligned.sortedByCoord.out.bam O3_W5_M2_CGGAAT_L008_c5ubdanxxAligned.sortedByCoord.out.bam O4_W1_F1_CAGATC_L002_c5ubdanxxAligned.sortedByCoord.out.bam O4_W1_F2_CATGGC_L002_c5ubdanxxAligned.sortedByCoord.out.bam O4_W1_M1_GCGCTA_L003_c5ubdanxxAligned.sortedByCoord.out.bam O4_W1_M2_TCATTC_L003_c5ubdanxxAligned.sortedByCoord.out.bam O4_W5_F1_CAGATC_L004_c5ubdanxxAligned.sortedByCoord.out.bam O4_W5_F2_CATGGC_L006_c5ubdanxxAligned.sortedByCoord.out.bam O4_W5_M1_TATAAT_L007_c5ubdanxxAligned.sortedByCoord.out.bam O4_W5_M2_CTAGCT_L008_c5ubdanxxAligned.sortedByCoord.out.bam O5_W1_F1_ATCACG_L002_c5ubdanxxAligned.sortedByCoord.out.bam O5_W1_F2_CATTTT_L002_c5ubdanxxAligned.sortedByCoord.out.bam O5_W1_M1_TAATCG_L003_c5ubdanxxAligned.sortedByCoord.out.bam O5_W1_M2_TCCCGA_L003_c5ubdanxxAligned.sortedByCoord.out.bam O5_W5_F1_ATCACG_L004_c5ubdanxxAligned.sortedByCoord.out.bam O5_W5_F2_CATTTT_L006_c5ubdanxxAligned.sortedByCoord.out.bam O5_W5_M1_GTGGCC_L004_c5ubdanxxAligned.sortedByCoord.out.bam O5_W5_M2_CTATAC_L008_c5ubdanxxAligned.sortedByCoord.out.bam
