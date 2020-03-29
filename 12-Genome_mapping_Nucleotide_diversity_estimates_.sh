#Mapping to reference genome to find out nucleotide diversity & check of Wolbachia presence

#1) Download and install: 
#PicardTools: https://broadinstitute.github.io/picard/ We used: version 2.18.27
#PoPoolation: https://sourceforge.net/projects/popoolation/

#2) Add additional sequences to reference genome
cd reference
#Additional fasta sequences and gtfs were obtained from ENSEMBL (https://www.ensembl.org/), accession numbers in file name and in methods: Wolbachia_NC_002978.6 Lactobacillus_CP000416.1 Acetobacter_AP011121.1
#Save into reference folder
#We then edited the headers of fasta files so that they are simplified 
#E.g. we chagned ">NC_002978.6 Wolbachia endosymbiont of Drosophila melanogaster, complete genome" to ">Wolbachia"

#Create reference fasta
cat rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta Wolbachia_NC_002978.6_edit.fasta Lactobacillus_CP000416.1_edit.fasta Acetobacter_AP011121.1_edit.fasta | sed 's/>Wolbachia/\n>Wolbachia/g' > rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta 
grep "^>" rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta | less

#3) Create cutom gtf in excel with all entries as "exon": This is to count all reads in mapping to TEs in the downstream analysis because some annotations are not completed (e.g. DMAURA only has UTRs annotated)
#Modify "te_features.gff" in the deviaTE library (https://github.com/W-L/deviaTE/tree/master/deviaTE/lib) and save as "custom_deviaTE_gtf_for_RNAseq.gtf" in reference folder
cat custom_deviaTE_gtf_for_RNAseq.gtf | sed -e "s/[[:space:]]\+/\t/g" | sed 's/\t"/ "/g;s/;\t/; /g' > custom_deviaTE_gtf_for_RNAseq_edit.gtf 

#We also edited gtf files of bacteria:
cat Wolbachia_NC_002978.6.gtf | grep -v "^#" | sed 's/Chromosome/Wolbachia/g' > Wolbachia_NC_002978.6_edit.gtf
cat Acetobacter_AP011121.1.gtf | grep -v "^#" | sed 's/Chromosome/Acetobacter/g' > Acetobacter_AP011121.1_edit.gtf 
cat Lactobacillus_CP000416.1.gtf | grep -v "^#" | sed 's/Chromosome/Lactobacillus/g' > Lactobacillus_CP000416.1_edit.gtf

#merge GTF files 
cat Wolbachia_NC_002978.6_edit.gtf Acetobacter_AP011121.1_edit.gtf Lactobacillus_CP000416.1_edit.gtf custom_deviaTE_gtf_for_RNAseq_edit.gtf > dmel-all-r6.27_TE_BACT_ADD.gtf

#If downstream errors occur, add new line at the end of the file

#4) Index with bwa
bwa index rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta

#Before next steps don't forget to check quality encoding of fastq sequences!
#If the quality scores contain character 0 it is phred+33. If the quality scores do not contain 0, it is phred+64. 
#FASTQC: Sanger encoding for Remolina, Carnes, Hoedjes
#Fabian: Illumina 1.5

#5) MAP USING BWA: FABIAN
cd Fabian/reads/unzipped/trimmed/
mkdir bwamem_mapped
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_Ra1.trim.fastq Cont_Ra2.trim.fastq > bwamem_mapped/Cont_Ra.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_Rb1.trim.fastq Cont_Rb2.trim.fastq > bwamem_mapped/Cont_Rb.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_La1.trim.fastq Sel_La2.trim.fastq > bwamem_mapped/Sel_La.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_Lb1.trim.fastq Sel_Lb2.trim.fastq > bwamem_mapped/Sel_Lb.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_2La1.trim.fastq Sel_2La2.trim.fastq > bwamem_mapped/Sel_2La.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_2Lb1.trim.fastq Sel_2Lb2.trim.fastq > bwamem_mapped/Sel_2Lb.sam

#6) MAP USING BWA: CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir bwamem_mapped

bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B1_r1.trim.fastq Cont_B1_r2.trim.fastq > bwamem_mapped/Cont_B1.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B2_r1.trim.fastq Cont_B2_r2.trim.fastq > bwamem_mapped/Cont_B2.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B3_r1.trim.fastq Cont_B3_r2.trim.fastq > bwamem_mapped/Cont_B3.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B4_r1.trim.fastq Cont_B4_r2.trim.fastq > bwamem_mapped/Cont_B4.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B5_r1.trim.fastq Cont_B5_r2.trim.fastq > bwamem_mapped/Cont_B5.sam

bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O1_r1.trim.fastq Sel_O1_r2.trim.fastq > bwamem_mapped/Sel_O1.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O2_r1.trim.fastq Sel_O2_r2.trim.fastq > bwamem_mapped/Sel_O2.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O3_r1.trim.fastq Sel_O3_r2.trim.fastq > bwamem_mapped/Sel_O3.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O4_r1.trim.fastq Sel_O4_r2.trim.fastq > bwamem_mapped/Sel_O4.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O5_r1.trim.fastq Sel_O5_r2.trim.fastq > bwamem_mapped/Sel_O5.sam

#7) MAP USING BWA: REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir bwamem_mapped

for i in {1..6}
do
for j in {1..3}
do
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_${j}_${i}_r1.trim.fastq Cont_${j}_${i}_r2.trim.fastq > bwamem_mapped/Cont_${j}_${i}.sam
bwa mem -t 20 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_${j}_${i}_r1.trim.fastq Sel_${j}_${i}_r2.trim.fastq > bwamem_mapped/Sel_${j}_${i}.sam
done
done

#8) MAP USING BWA: HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir bwamem_mapped

for i in {1..4}
do
echo "CE${i}"
echo "CP${i}"
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta CE${i}_r1.trim.fastq CE${i}_r2.trim.fastq > bwamem_mapped/CE${i}.sam
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta CP${i}_r1.trim.fastq CP${i}_r2.trim.fastq > bwamem_mapped/CP${i}.sam

echo "LE${i}"
echo "LP${i}"
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta LE${i}_r1.trim.fastq LE${i}_r2.trim.fastq > bwamem_mapped/LE${i}.sam
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta LP${i}_r1.trim.fastq LP${i}_r2.trim.fastq > bwamem_mapped/LP${i}.sam

echo "HE${i}"
echo "HP${i}"
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta HE${i}_r1.trim.fastq HE${i}_r2.trim.fastq > bwamem_mapped/HE${i}.sam
bwa mem -t 24 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta HP${i}_r1.trim.fastq HP${i}_r2.trim.fastq > bwamem_mapped/HP${i}.sam
done

#9) SAM TO BAM, SORT, RENAME: FABIAN
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/
for file in *.sam
do 
outname="$(echo $file | sed 's/\.sam/\.bam/g')"
samtools view -Sb $file > ${outname}
done

for file in *.bam
do 
outname="$(echo $file | sed 's/\.bam/\_sorted\.bam/g')"
echo "$outname"
samtools sort $file ${outname}
done

for file in *.bam.bam; do
    mv "$file" "$(basename "$file" .bam.bam).bam"
done

#10) SAM TO BAM, SORT, RENAME: REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
for file in *.sam
do 
outname="$(echo $file | sed 's/\.sam/\.bam/g')"
samtools view -Sb $file > ${outname}
done

for file in *.bam
do 
outname="$(echo $file | sed 's/\.bam/\_sorted\.bam/g')"
echo "$outname"
samtools sort $file ${outname}
done

for file in *.bam.bam; do
    mv "$file" "$(basename "$file" .bam.bam).bam"
done


#11) SAM TO BAM, SORT, RENAME: CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
for file in *.sam
do 
outname="$(echo $file | sed 's/\.sam/\.bam/g')"
samtools view -Sb $file > ${outname}
done

for file in *.bam
do 
outname="$(echo $file | sed 's/\.bam/\_sorted\.bam/g')"
echo "$outname"
samtools sort $file ${outname}
done

for file in *.bam.bam; do
    mv "$file" "$(basename "$file" .bam.bam).bam"
done

#12) SAM TO BAM, SORT, RENAME: HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped

for file in *.sam
do 
outname="$(echo $file | sed 's/\.sam/\.bam/g')"
samtools view -Sb $file > ${outname}
done

for file in *.bam
do 
outname="$(echo $file | sed 's/\.bam/\_sorted\.bam/g')"
echo "$outname"
samtools sort $file ${outname}
done

#rename files
for file in *.bam.bam; do
    mv "$file" "$(basename "$file" .bam.bam).bam"
done

#13) REMOVE DUPLICATES: FABIAN
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/
mkdir DupRem
mkdir tmp
mkdir tmp2
mkdir tmp3

java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp REMOVE_DUPLICATES=true I=Cont_Ra_sorted.bam O=DupRem/Cont_Ra_sorted_DupRem.bam M=DupRem/Cont_Ra_sorted_DupRemMetrics.bam
java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp REMOVE_DUPLICATES=true I=Cont_Rb_sorted.bam O=DupRem/Cont_Rb_sorted_DupRem.bam M=DupRem/Cont_Rb_sorted_DupRemMetrics.bam
java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp2 REMOVE_DUPLICATES=true I=Sel_La_sorted.bam O=DupRem/Sel_La_sorted_DupRem.bam M=DupRem/Sel_La_sorted_DupRemMetrics.bam
java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp2 REMOVE_DUPLICATES=true I=Sel_2La_sorted.bam O=DupRem/Sel_2La_sorted_DupRem.bam M=DupRem/Sel_2La_sorted_DupRemMetrics.bam
java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp3 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp3 REMOVE_DUPLICATES=true I=Sel_Lb_sorted.bam O=DupRem/Sel_Lb_sorted_DupRem.bam M=DupRem/Sel_Lb_sorted_DupRemMetrics.bam
java -Xmx20g -Djava.io.tmpdir=`pwd`/tmp3 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates TMP_DIR=`pwd`/tmp3 REMOVE_DUPLICATES=true I=Sel_2Lb_sorted.bam O=DupRem/Sel_2Lb_sorted_DupRem.bam M=DupRem/Sel_2Lb_sorted_DupRemMetrics.bam

#13) REMOVE DUPLICATES: REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
mkdir DupRem

for file in *_sorted.bam
do
outname="$(echo $file | sed 's/\.bam/_DupRem\.bam/g')"
outname2="$(echo $file | sed 's/\.bam/_DupRemMetrics\.txt/g')"
java -Xmx8g -jar /PATH/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$file O=DupRem/${outname} M=DupRem/${outname2}
done

#14) REMOVE DUPLICATES: CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
mkdir DupRem

for file in *_sorted.bam
do
outname="$(echo $file | sed 's/\.bam/_DupRem\.bam/g')"
outname2="$(echo $file | sed 's/\.bam/_DupRemMetrics\.txt/g')"
java -Xmx20g -jar /PATH/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$file O=DupRem/${outname} M=DupRem/${outname2}
done

#14) REMOVE DUPLICATES: HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
mkdir tmp
mkdir tmp2
mkdir tmp3
mkdir DupRem

for file in H*_sorted.bam
do
echo "$file"
outname="$(echo $file | sed 's/\.bam/_DupRem\.bam/g')"
outname2="$(echo $file | sed 's/\.bam/_DupRemMetrics\.txt/g')"
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$file O=DupRem/${outname} M=DupRem/${outname2}
done

for file in C*_sorted.bam
do
echo "$file"
outname="$(echo $file | sed 's/\.bam/_DupRem\.bam/g')"
outname2="$(echo $file | sed 's/\.bam/_DupRemMetrics\.txt/g')"
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$file O=DupRem/${outname} M=DupRem/${outname2}
done

for file in L*_sorted.bam
do
echo "$file"
outname="$(echo $file | sed 's/\.bam/_DupRem\.bam/g')"
outname2="$(echo $file | sed 's/\.bam/_DupRemMetrics\.txt/g')"
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp3 -jar /PATH/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$file O=DupRem/${outname} M=DupRem/${outname2}
done

#15) merge files from multiple lanes - only for Remolina!
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped
samtools merge Cont_1_DupRem_merge.bam Cont_1_1_sorted_DupRem.bam Cont_1_2_sorted_DupRem.bam Cont_1_3_sorted_DupRem.bam Cont_1_4_sorted_DupRem.bam Cont_1_5_sorted_DupRem.bam Cont_1_6_sorted_DupRem.bam
samtools merge Cont_2_DupRem_merge.bam Cont_2_1_sorted_DupRem.bam Cont_2_2_sorted_DupRem.bam Cont_2_3_sorted_DupRem.bam Cont_2_4_sorted_DupRem.bam Cont_2_5_sorted_DupRem.bam Cont_2_6_sorted_DupRem.bam
samtools merge Cont_3_DupRem_merge.bam Cont_3_1_sorted_DupRem.bam Cont_3_2_sorted_DupRem.bam Cont_3_3_sorted_DupRem.bam Cont_3_4_sorted_DupRem.bam Cont_3_5_sorted_DupRem.bam Cont_3_6_sorted_DupRem.bam
samtools merge Sel_1_DupRem_merge.bam Sel_1_1_sorted_DupRem.bam Sel_1_2_sorted_DupRem.bam Sel_1_3_sorted_DupRem.bam Sel_1_4_sorted_DupRem.bam Sel_1_5_sorted_DupRem.bam Sel_1_6_sorted_DupRem.bam
samtools merge Sel_2_DupRem_merge.bam Sel_2_1_sorted_DupRem.bam Sel_2_2_sorted_DupRem.bam Sel_2_3_sorted_DupRem.bam Sel_2_4_sorted_DupRem.bam Sel_2_5_sorted_DupRem.bam Sel_2_6_sorted_DupRem.bam
samtools merge Sel_3_DupRem_merge.bam Sel_3_1_sorted_DupRem.bam Sel_3_2_sorted_DupRem.bam Sel_3_3_sorted_DupRem.bam Sel_3_4_sorted_DupRem.bam Sel_3_5_sorted_DupRem.bam Sel_3_6_sorted_DupRem.bam


#16) REMOVE READS WITH LOW MAPPING QUALITY AND WITHOUT PROPER PAIRS, CREATE PILEUP
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem
mkdir pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_Ra_sorted_DupRem.bam > pileup/Cont_Ra_filtered.pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_Rb_sorted_DupRem.bam > pileup/Cont_Rb_filtered.pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_La_sorted_DupRem.bam > pileup/Sel_La_filtered.pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_Lb_sorted_DupRem.bam > pileup/Sel_Lb_filtered.pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_2La_sorted_DupRem.bam > pileup/Sel_2La_filtered.pileup
samtools mpileup -6 -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_2Lb_sorted_DupRem.bam > pileup/Sel_2Lb_filtered.pileup

#REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem
mkdir pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_1_DupRem_merge.bam > pileup/Cont_1_DupRem_merge_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_2_DupRem_merge.bam > pileup/Cont_2_DupRem_merge_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_3_DupRem_merge.bam > pileup/Cont_3_DupRem_merge_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_1_DupRem_merge.bam > pileup/Sel_1_DupRem_merge_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_2_DupRem_merge.bam > pileup/Sel_2_DupRem_merge_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_3_DupRem_merge.bam > pileup/Sel_3_DupRem_merge_filtered.pileup

cut -f6 pileup/Sel_1_DupRem_merge_filtered.pileup | grep "0" | less

#CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem
mkdir pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B1_sorted_DupRem.bam > pileup/Cont_B1_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B2_sorted_DupRem.bam > pileup/Cont_B2_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B3_sorted_DupRem.bam > pileup/Cont_B3_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B4_sorted_DupRem.bam > pileup/Cont_B4_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Cont_B5_sorted_DupRem.bam > pileup/Cont_B5_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O1_sorted_DupRem.bam > pileup/Sel_O1_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O2_sorted_DupRem.bam > pileup/Sel_O2_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O3_sorted_DupRem.bam > pileup/Sel_O3_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O4_sorted_DupRem.bam > pileup/Sel_O4_DupRem_filtered.pileup
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta Sel_O5_sorted_DupRem.bam > pileup/Sel_O5_DupRem_filtered.pileup

#HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem
mkdir pileup
for file in *DupRem.bam
do
echo "$file"
outname="$(echo $file | sed 's/\.bam/_filtered\.pileup/g')"
echo "${outname}"
samtools mpileup -q 20 --rf 0x2 --ff 0x4 --ff 0x8 -f reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE_BACT.fasta $file > pileup/${outname}
done


#17) CALCULATE GENOME WIDE COVERAGE
#FABIAN
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir coverage
for file in *.pileup
do
grep -wE '^2L|^2R|^3L|^3R|^X|^4' $file > $file.Dmel
done

for file in *.Dmel
do
cat $file | cut -f4 | sort -nr | uniq > coverage/$file.uniqcov
cat $file | cut -f4 | sort -nr | uniq | head -1 > coverage/$file.maxcov
cat $file | awk '{sum+=\$4} END { print sum/NR}' > coverage/$file.averagecov
done

#REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir coverage
for file in *.pileup
do
grep -wE '^2L|^2R|^3L|^3R|^X|^4' $file > $file.Dmel
done

for file in *.Dmel
do
cat $file | cut -f4 | sort -nr | uniq > coverage/$file.uniqcov
cat $file | cut -f4 | sort -nr | uniq | head -1 > coverage/$file.maxcov
cat $file | awk '{sum+=\$4} END { print sum/NR}' > coverage/$file.averagecov
done

#CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir coverage
for file in *.pileup
do
grep -wE '^2L|^2R|^3L|^3R|^X|^4' $file > $file.Dmel
done

for file in *.Dmel
do
cat $file | cut -f4 | sort -nr | uniq > coverage/$file.uniqcov
cat $file | cut -f4 | sort -nr | uniq | head -1 > coverage/$file.maxcov
cat $file | awk '{sum+=\$4} END { print sum/NR}' > coverage/$file.averagecov
done

#HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir coverage
for file in *.pileup
do
grep -wE '^2L|^2R|^3L|^3R|^X|^4' $file > $file.Dmel
done

for file in *.Dmel
do
cat $file | cut -f4 | sort -nr | uniq > coverage/$file.uniqcov
cat $file | cut -f4 | sort -nr | uniq | head -1 > coverage/$file.maxcov
cat $file | awk '{sum+=\$4} END { print sum/NR}' > coverage/$file.averagecov
done

#18) Check average coverage output file
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem/pileup/coverage
cat *.averagecov | grep -v "Wolbachia" | sed 's/Average \=  //g' | awk '{sum+=$1} END { print "Average across all = " sum/NR}'
# Fabian Average across all = 161.916 

for study in "Carnes" "Hoedjes" "Remolina"
do
echo $study
cat $study/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup/coverage/*.averagecov | grep -v "Wolbachia" | sed 's/Average \=  //g' | awk '{sum+=$1} END { print "Average across all = " sum/NR}'
done
# Carnes
# Average across all = 22.5829
# Hoedjes
# Average across all = 101.431
# Remolina
# Average across all = 41.1345

#####################################################
#19) ESTIMATE PI

#FABIAN 
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir pi

#Define max-coverage cut-offs for each population
for file in *.pileup
do
cov="$(cat coverage/$file.Dmel.averagecov | sed 's/Average \=//g')"
cov2="$(echo "$cov * 2" | bc -l)" 
cov3="$(printf "%0.f" $cov2)"
echo "$file $cov * 2 = $cov2 | $cov3"
done

#Calculate pi
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_Ra_filtered.pileup --output pi/Cont_Ra_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 330 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_Rb_filtered.pileup --output pi/Cont_Rb_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 340 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_La_filtered.pileup --output pi/Sel_La_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 315 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_Lb_filtered.pileup --output pi/Sel_Lb_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 317 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2La_filtered.pileup --output pi/Sel_2La_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 286 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2Lb_filtered.pileup --output pi/Sel_2Lb_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 355 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6

cd pi
for file in *pi
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.pi//g' > $file.new
done
cat *.new > Fabian_all_pi_100kb.new

#REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir pi

#Define max-coverage cut-offs for each population
for file in *.pileup
do
cov="$(cat coverage/$file.Dmel.averagecov | sed 's/Average \=//g')"
cov2="$(echo "$cov * 2" | bc -l)" 
cov3="$(printf "%0.f" $cov2)"
echo "$file $cov * 2 = $cov2 | $cov3"
done

#Calculate pi
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_1_DupRem_merge_filtered.pileup --output pi/Cont_1_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 81 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_2_DupRem_merge_filtered.pileup --output pi/Cont_2_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 107 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_3_DupRem_merge_filtered.pileup --output pi/Cont_3_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 79 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_1_DupRem_merge_filtered.pileup --output pi/Sel_1_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 75 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2_DupRem_merge_filtered.pileup --output pi/Sel_2_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 85 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_3_DupRem_merge_filtered.pileup --output pi/Sel_3_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 67 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6

cd pi
for file in *pi
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.pi//g' > $file.new
done
cat *.new > Remolina_all_pi_100kb.new

#CARNES - SANGER ENCODING
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir pi

#Define max-coverage cut-offs for each population
for file in *.pileup
do
cov="$(cat coverage/$file.Dmel.averagecov | sed 's/Average \=//g')"
cov2="$(echo "$cov * 2" | bc -l)" 
cov3="$(printf "%0.f" $cov2)"
echo "$file $cov * 2 = $cov2 | $cov3"
done

#Calculate pi
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B1_DupRem_filtered.pileup --output pi/Cont_B1_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 34 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B2_DupRem_filtered.pileup --output pi/Cont_B2_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 48 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B3_DupRem_filtered.pileup --output pi/Cont_B3_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 32 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B4_DupRem_filtered.pileup --output pi/Cont_B4_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 20 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B5_DupRem_filtered.pileup --output pi/Cont_B5_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 26 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O1_DupRem_filtered.pileup --output pi/Sel_O1_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 49 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O2_DupRem_filtered.pileup --output pi/Sel_O2_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 83 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O3_DupRem_filtered.pileup --output pi/Sel_O3_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 74 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O4_DupRem_filtered.pileup --output pi/Sel_O4_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 49 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O5_DupRem_filtered.pileup --output pi/Sel_O5_100kb.pi --min-count 2 --min-coverage 5 --max-coverage 38 --window-size 100000 --step-size 100000 --pool-size 100 --measure pi --fastq-type sanger --min-covered-fraction 0.6

cd pi
for file in *pi
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.pi//g' > $file.new
done
cat *.new > Carnes_all_pi_100kb.new

#HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir pi

for file in *.pileup
do
cov="$(cat coverage/$file.Dmel.averagecov)"
cov2="$(echo "$cov * 2" | bc -l)" 
cov3="$(printf "%0.f" $cov2)"
echo "$cov * 2 = $cov2 | $cov3"
outname="$(echo $file | sed 's/\_sorted\_DupRem\_filtered\.pileup/\_100kb\.pi/g')"
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input $file --output pi/${outname} --min-count 2 --min-coverage 5 --max-coverage $cov3 --window-size 100000 --step-size 100000 --pool-size 250 --measure pi --fastq-type sanger --min-covered-fraction 0.6
done

cd pi
for file in *pi
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.pi//g' > $file.new
done
cat *.new > Hoedjes_all_pi_100kb.new


#####################################################
#20) ESTIMATE THETA

#FABIAN
cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir theta

#Calculate theta
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_Ra_filtered.pileup --output theta/Cont_Ra_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 330 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_Rb_filtered.pileup --output theta/Cont_Rb_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 340 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_La_filtered.pileup --output theta/Sel_La_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 315 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_Lb_filtered.pileup --output theta/Sel_Lb_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 317 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2La_filtered.pileup --output theta/Sel_2La_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 286 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2Lb_filtered.pileup --output theta/Sel_2Lb_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 355 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6

cd theta
for file in *theta
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.theta//g' > $file.new
done
cat *.new > Fabian_all_theta_100kb.new

#REMOLINA
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir theta

#Calculate theta
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_1_DupRem_merge_filtered.pileup --output theta/Cont_1_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 81 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_2_DupRem_merge_filtered.pileup --output theta/Cont_2_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 107 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_3_DupRem_merge_filtered.pileup --output theta/Cont_3_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 79 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_1_DupRem_merge_filtered.pileup --output theta/Sel_1_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 75 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_2_DupRem_merge_filtered.pileup --output theta/Sel_2_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 85 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_3_DupRem_merge_filtered.pileup --output theta/Sel_3_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 67 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6

cd theta
for file in *theta
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.theta//g' > $file.new
done
cat *.new > Remolina_all_theta_100kb.new

#CARNES
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir theta

#Calculate theta
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B1_DupRem_filtered.pileup --output thetaCont_B1_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 34 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B2_DupRem_filtered.pileup --output theta/Cont_B2_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 48 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B3_DupRem_filtered.pileup --output theta/Cont_B3_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 32 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B4_DupRem_filtered.pileup --output theta/Cont_B4_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 20 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Cont_B5_DupRem_filtered.pileup --output theta/Cont_B5_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 26 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O1_DupRem_filtered.pileup --output theta/Sel_O1_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 49 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O2_DupRem_filtered.pileup --output theta/Sel_O2_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 83 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O3_DupRem_filtered.pileup --output theta/Sel_O3_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 74 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O4_DupRem_filtered.pileup --output theta/Sel_O4_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 49 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input Sel_O5_DupRem_filtered.pileup --output theta/Sel_O5_100kb.theta --min-count 2 --min-coverage 5 --max-coverage 38 --window-size 100000 --step-size 100000 --pool-size 100 --measure theta --fastq-type sanger --min-covered-fraction 0.6

cd theta
for file in *theta
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.theta//g' > $file.new
done
cat *.new > Carnes_all_theta_100kb.new

#HOEDJES
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
mkdir theta

#Calculate theta
for file in *.pileup
do
cov="$(cat coverage/$file.Dmel.averagecov)"
cov2="$(echo "$cov * 2" | bc -l)" 
cov3="$(printf "%0.f" $cov2)"
echo "$cov * 2 = $cov2 | $cov3"
outname="$(echo $file | sed 's/\_sorted\_DupRem\_filtered\.pileup/\_100kb\.theta/g')"
perl /PATH/popoolation_1.2.2/Variance-sliding.pl --input $file --output theta/${outname} --min-count 2 --min-coverage 5 --max-coverage $cov3 --window-size 100000 --step-size 100000 --pool-size 250 --measure theta --fastq-type sanger --min-covered-fraction 0.6
done

cd theta
for file in *theta
do
awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $file | sed 's/\_100kb\.theta//g' > $file.new
done
cat *.new > Hoedjes_all_theta_100kb.new


#####################################################
#21) Check Wolbachia presence

cd Fabian/reads/unzipped/trimmed/bwamem_mapped/DupRem/pileup
for file in *.pileup
do
echo "$file"
grep 'Wolbachia' $file > $file.Wolbachia
done
for file in *.Wolbachia
do
cat $file | awk '{sum+=$4} END { print "Average (Wolbachia) = ",sum/NR}' > coverage/$file.averagecov &
done
for file in *Wolbachia*
do
echo "$file"
cat $file
done

cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
for file in *.pileup
do
echo "$file"
grep 'Wolbachia' $file > $file.Wolbachia
done
for file in *.Wolbachia
do
cat $file | awk '{sum+=$4} END { print "Average (Wolbachia) = ",sum/NR}' > coverage/$file.averagecov &
done
for file in *Wolbachia*
do
echo "$file"
cat $file
done

cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
for file in *.pileup
do
echo "$file"
grep 'Wolbachia' $file > $file.Wolbachia
done
for file in *.Wolbachia
do
cat $file | awk '{sum+=$4} END { print "Average (Wolbachia) = ",sum/NR}' > coverage/$file.averagecov &
done
for file in *Wolbachia*
do
echo "$file"
cat $file
done

cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/bwamem_mapped/DupRem/pileup
for file in *.pileup
do
echo "$file"
grep 'Wolbachia' $file > $file.Wolbachia
done
for file in *.Wolbachia
do
cat $file | awk '{sum+=$4} END { print "Average (Wolbachia) = ",sum/NR}' > coverage/$file.averagecov &
done
for file in *Wolbachia*
do
echo "$file"
cat $file
done
