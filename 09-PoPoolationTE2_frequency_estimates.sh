#POPOOLATION TE2 ANALYSIS

#We followed the recommended pipeline by the authors of Kofler et al at: https://sourceforge.net/p/popoolation-te2/wiki/Home/
#Also see: Kofler et al. 2016, MBE

#PREPARATORY WORK

#1) Downlaod PoPoolation TE2
wget -O popte2-v1.10.03.jar https://sourceforge.net/projects/popoolation-te2/files/popte2-v1.10.03.jar/download

#2) Install RepeatMasker
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-9.tar.gz
tar -zxvf RepeatMasker-open-4-0-9.tar.gz

#3) Download tandem repeat finder http://tandem.bu.edu/trf/trf.download.html
#Go into folder of TandemRepeatFinder
#Rename file to: trf
#make file executable: chmod +x+u trf

#4) Download RMblast
wget http://www.repeatmasker.org/rmblast-2.9.0+-x64-linux.tar.gz
tar zxvf rmblast-2.9.0-x64-linux.tar.gz

#5) Configure RepeatMasker
cd /PATH/RepeatMasker
perl ./configure
# Enter path to trf
# Choose 2 for RMblast
# Enter path to RMblast: /PATH/rmblast-2.9.0
# choose 5
# Done to finish

#6) Get D.melanogaster reference sequence (and gtf for later RNAseq analysis):
mkdir reference
cd reference	
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.27.fasta.gz
zcat dmel-all-chromosome-r6.27.fasta.gz > dmel-all-chromosome-r6.27.fasta

wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gtf/dmel-all-r6.27.gtf.gz
zcat dmel-all-r6.27.gtf.gz > dmel-all-r6.27.gtf

#7) Prepare TE consensus sequence library
https://github.com/W-L/deviaTE/blob/master/deviaTE/lib/te_library
#Open and delete sequences of D. melanogaster single copy genes: 
# >Dmel_rpl32
# >Dmel_RpII140
# >Dmel_Act5C
# >Dmel_piwi
# >Dmel_p53
#Save as: te_library_deviaTE_noDmel.fasta

#8) Edit te_library_deviaTE_noDmel.fasta:
sed 's/\(>.*\)/\1&/' te_library_deviaTE_noDmel.fasta | sed 's/\(>.*\)>/\1\t/g' > te_library_deviaTE_noDmel_edit.fasta 

#9) Repeat-mask reference genome
mkdir rep_masked
cd reference/rep_masked
perl /PATH/RepeatMasker -gccalc -s -cutoff 200 -no_is -nolow -norna -gff -u -pa 18 -lib te_library_deviaTE_noDmel_edit.fasta reference/dmel-all-chromosome-r6.27.fasta

#10) Move output to new folder
cd reference
mv dmel-all-chromosome-r6.27.fasta.cat.gz reference/rep_masked
mv dmel-all-chromosome-r6.27.fasta.masked reference/rep_masked
mv dmel-all-chromosome-r6.27.fasta.tbl reference/rep_masked
mv dmel-all-chromosome-r6.27.fasta.out reference/rep_masked
mv dmel-all-chromosome-r6.27.fasta.out.gff reference/rep_masked

#11) merge masked file with TE fasta file
cat reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked reference/te_library_deviaTE_noDmel_edit.fasta > reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta

#12) Create TE hierarchy
cat reference/te_library_deviaTE_noDmel_edit.fasta | grep '^>' | perl -pe 's/>//' > reference/te-hierarchy.txt
#edit in excel using information from, e.g. Casey Bergman's TE annotations: https://github.com/bergmanlab/transposons

#13) also make sure cols are tab-separated 
awk '{print $1"\t"$2"\t"$3}' reference/te-hierarchy_edit.txt > reference/te-hierarchy_edit.txt
#Filte should look like this:
# id     family  order
# 412     412     LTR
# 1360    1360    TIR
# ACCORD  ACCORD  LTR
# ...

#15) Install Cutadapt following instructions here: https://cutadapt.readthedocs.io/en/stable/
#We used: cutadapt 2.0 with Python 3.6.3

#16) Download and install bwa: https://sourceforge.net/projects/bio-bwa/files/
#We used: 0.7.17-r1198-dirty

#17) Download and install samtools: https://sourceforge.net/projects/samtools/files/
#We used: Version: 0.1.19-44428cd

#18) index reference
bwa index reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta

##################################################################################################################
#CARNES2015
##################################################################################################################

#1) Trim reads in paired end mode
cd Carnes/SRA/reads/gzip_fastq/unzipped 
mkdir trimmed

for i in {1..5}
do
cutadapt -j 15 --minimum-length 50 -q 18 -o Cont_B${i}_r1.trim.fastq -p Cont_B${i}_r2.trim.fastq Cont_B${i}_r1.fastq Cont_B${i}_r2.fastq 
cutadapt -j 15 --minimum-length 50 -q 18 -o Sel_O${i}_r1.trim.fastq -p Sel_O${i}_r2.trim.fastq Sel_O${i}_r1.fastq Sel_O${i}_r2.fastq
done

mv Carnes/SRA/reads/gzip_fastq/unzipped/*.trim.fastq Carnes/SRA/reads/gzip_fastq/unzipped/trimmed

#2) Change headers: read pair 1 end with /1 and pair 2 with /2 (otherwise error in SE2PE below)
for i in {1..5}
do
sed 's/ length=//g' Cont_B${i}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > Cont_B${i}_r1.trim.headerEdit.fastq
sed 's/ length=//g' Cont_B${i}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > Cont_B${i}_r2.trim.headerEdit.fastq
sed 's/ length=//g' Sel_O${i}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > Sel_O${i}_r1.trim.headerEdit.fastq
sed 's/ length=//g' Sel_O${i}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > Sel_O${i}_r2.trim.headerEdit.fastq
done

#3) Map using bwa
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir mapped
for file in *.trim.headerEdit.fastq
do
bwa bwasw -t 18 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta $file > mapped/$file.sam
done

#4) Restore PE info
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed
for i in {1..5}
do
java -Xmx4g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Cont_B${i}_r1.trim.headerEdit.fastq  --fastq2 Cont_B${i}_r2.trim.headerEdit.fastq --bam1 mapped/Cont_B${i}_r1.trim.headerEdit.fastq.sam --bam2 mapped/Cont_B${i}_r2.trim.headerEdit.fastq.sam --sort --output mapped/Cont_B${i}.sort.bam
java -Xmx4g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_O${i}_r1.trim.headerEdit.fastq  --fastq2 Sel_O${i}_r2.trim.headerEdit.fastq --bam1 mapped/Sel_O${i}_r1.trim.headerEdit.fastq.sam --bam2 mapped/Sel_O${i}_r2.trim.headerEdit.fastq.sam --sort --output mapped/Sel_O${i}.sort.bam
done

#5) Generate the ppileup
mkdir Carnes/popTE2
cd Carnes/SRA/reads/gzip_fastq/unzipped/trimmed/mapped/
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar ppileup --bam Cont_B1.sort.bam --bam Cont_B2.sort.bam --bam Cont_B3.sort.bam --bam Cont_B4.sort.bam --bam Cont_B5.sort.bam --bam Sel_O1.sort.bam --bam Sel_O2.sort.bam --bam Sel_O3.sort.bam --bam Sel_O4.sort.bam --bam Sel_O5.sort.bam --map-qual 15 --hier reference/te-hierarchy_edit.txt --output Carnes/popTE2/Cont_Sel.ppileup.gz

#6) Joint analysis to get TE positions and frequencies
mkdir joint_analysis
cd Carnes/popTE2

#7) Get signatures - min count 5 (10 populations)
cd Carnes/popTE2/

java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar identifySignatures --ppileup Cont_Sel.ppileup.gz --mode joint --output joint_analysis/Cont_Sel.ppileup.signatures --min-count 5 --signature-window fix500 --min-valley fix150 --detailed-log
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar frequency --ppileup Cont_Sel.ppileup.gz --signature joint_analysis/Cont_Sel.ppileup.signatures --output joint_analysis/Cont_Sel.freqsig
java -Xmx7g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature joint_analysis/Cont_Sel.freqsig --detailed-log --output-detail medium --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --min-distance -200 --max-distance 300 --output joint_analysis/Cont_Sel.teinsertions

#8) Filter low quality TE insertions
cd Carnes/popTE2/joint_analysis
mkdir filter
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar filterSignatures --input Cont_Sel.freqsig --output filter/Cont_Sel.filter.freqsig --max-otherte-count 2 --max-structvar-count 2
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature filter/Cont_Sel.filter.freqsig --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --detailed-log --output-detail medium --min-distance -200 --max-distance 300 --output filter/Cont_Sel.filter.teinsertions


##################################################################################################################
#FABIAN2018
##################################################################################################################
#IMPORTANT: We worked with the original sequences (i.e. not downloaded from repository)
#If you download sequences from the repository some adjustments might be necessary
#This might include changing options of quality encoding AND changing headers of reads in fastq files (pair 1 end with /1 and pair 2 with /2, otherwise error in SE2PE below)

#1) Trim reads in paired end mode
cd Fabian/reads/unzipped/
mkdir trimmed

cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Cont_Ra1.trim.fastq -p Cont_Ra2.trim.fastq Cont_Ra1.fastq Cont_Ra2.fastq
cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Cont_Rb1.trim.fastq -p Cont_Rb2.trim.fastq Cont_Rb1.fastq Cont_Rb2.fastq
cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Sel_La1.trim.fastq -p Sel_La2.trim.fastq Sel_La1.fastq Sel_La2.fastq
cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Sel_Lb1.trim.fastq -p Sel_Lb2.trim.fastq Sel_Lb1.fastq Sel_Lb2.fastq
cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Sel_2La1.trim.fastq -p Sel_2La2.trim.fastq Sel_2La1.fastq Sel_2La2.fastq
cutadapt -j 15 --quality-base=64 --minimum-length 50 -q 18 -o Sel_2Lb1.trim.fastq -p Sel_2Lb2.trim.fastq Sel_2Lb1.fastq Sel_2Lb2.fastq

#2) Map using bwa
cd Fabian/reads/unzipped/trimmed
mkdir mapped
for file in *.fastq.trim
do
bsub -o mapped/std_map.txt -e mapped/err_map.txt -R "rusage[mem=9000]" -M 10000 -n 18 "bwa bwasw -t 18 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta $file > mapped/$file.sam"
done

#3) Restore PE info
cd Fabian/reads/unzipped/trimmed/mapped
mkdir restorePE

#Create tmp folders, otherwise possibly error indicating disk space is too low
mkdir Fabian/reads/unzipped/trimmed/tmp
mkdir Fabian/reads/unzipped/trimmed/tmp1
mkdir Fabian/reads/unzipped/trimmed/tmp2
mkdir Fabian/reads/unzipped/trimmed/tmp3
mkdir Fabian/reads/unzipped/trimmed/tmp4
mkdir Fabian/reads/unzipped/trimmed/tmp5

cd Fabian/reads/unzipped/trimmed/
java -Xmx5g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Cont_Ra1.trim.fastq --fastq2 Cont_Ra2.trim.fastq --bam1 mapped/Cont_Ra1.trim.fastq.sam --bam2 mapped/Cont_Ra2.trim.fastq.sam --sort --output mapped/restorePE/Cont_Ra.sort.bam
java -Xmx5g -Djava.io.tmpdir=`pwd`/tmp1 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Cont_Rb1.trim.fastq --fastq2 Cont_Rb2.trim.fastq --bam1 mapped/Cont_Rb1.trim.fastq.sam --bam2 mapped/Cont_Rb2.trim.fastq.sam --sort --output mapped/restorePE/Cont_Rb.sort.bam
java -Xmx5g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_La1.trim.fastq --fastq2 Sel_La2.trim.fastq --bam1 mapped/Sel_La1.trim.fastq.sam --bam2 mapped/Sel_La2.trim.fastq.sam --sort --output mapped/restorePE/Sel_La.sort.bam
java -Xmx5g -Djava.io.tmpdir=`pwd`/tmp3 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_Lb1.trim.fastq --fastq2 Sel_Lb2.trim.fastq --bam1 mapped/Sel_Lb1.trim.fastq.sam --bam2 mapped/Sel_Lb2.trim.fastq.sam --sort --output mapped/restorePE/Sel_Lb.sort.bam
java -Xmx5g -Djava.io.tmpdir=`pwd`/tmp4 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_2La1.trim.fastq --fastq2 Sel_2La2.trim.fastq --bam1 mapped/Sel_2La1.trim.fastq.sam --bam2 mapped/Sel_2La2.trim.fastq.sam --sort --output mapped/restorePE/Sel_2La.sort.bam
java -Xmx8g -Djava.io.tmpdir=`pwd`/tmp5 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_2Lb1.trim.fastq --fastq2 Sel_2Lb2.trim.fastq --bam1 mapped/Sel_2Lb1.trim.fastq.sam --bam2 mapped/Sel_2Lb2.trim.fastq.sam --sort --output mapped/restorePE/Sel_2Lb.sort.bam

#5) Generate the ppileup
mkdir Fabian/popTE2
java -Xmx8g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar ppileup --bam Cont_Ra.sort.bam --bam Cont_Rb.sort.bam --bam Sel_La.sort.bam --bam Sel_Lb.sort.bam --bam Sel_2La.sort.bam --bam Sel_2Lb.sort.bam --map-qual 15 --hier reference/te-hierarchy_edit.txt --output Fabian/popTE2/Cont_Sel_Fabian.ppileup.gz

#6) Joint analysis to get TE positions and frequencies
cd Fabian/popTE2
mkdir joint_analysis

#7) Get signatures - min count 3 (6 populations)
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar identifySignatures --ppileup Cont_Sel_Fabian.ppileup.gz --mode joint --output joint_analysis/Cont_Sel_Fabian.ppileup.signatures --min-count 3 --signature-window fix500 --min-valley fix150 --detailed-log
java -Xmx5g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar frequency --ppileup Cont_Sel_Fabian.ppileup.gz --signature joint_analysis/Cont_Sel_Fabian.ppileup.signatures --output joint_analysis/Cont_Sel_Fabian.freqsig
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature joint_analysis/Cont_Sel_Fabian.freqsig --detailed-log --output-detail medium --ref-genome reference/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --min-distance -200 --max-distance 300 --output joint_analysis/less_cons/Cont_Sel_Fabian.teinsertions

#8) Filter low quality TE insertions
cd Fabian/popTE2/joint_analysis
mkdir filter
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar filterSignatures --input Cont_Sel_Fabian.freqsig --output filter/Cont_Sel_Fabian.filter.freqsig --max-otherte-count 2 --max-structvar-count 2
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature filter/Cont_Sel_Fabian.filter.freqsig --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --detailed-log --output-detail medium --min-distance -200 --max-distance 300 --output filter/Cont_Sel_Fabian.filter.teinsertions


##################################################################################################################
#REMOLINA2012
##################################################################################################################

#1) Trim reads in paired end mode
cd Remolina/SRA/reads/gzip_fastq/unzipped
mkdir trimmed

for i in {1..6}
do
for j in {1..3}
do
cutadapt -j 15 --minimum-length 50 -q 18 -o Cont_${j}_${i}_r1.trim.fastq -p Cont_${j}_${i}_r2.trim.fastq Cont_${j}_${i}_r1.fastq Cont_${j}_${i}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o Sel_${j}_${i}_r1.trim.fastq -p Sel_${j}_${i}_r2.trim.fastq Sel_${j}_${i}_r1.fastq Sel_${j}_${i}_r2.fastq
done
done

mv Remolina/SRA/reads/gzip_fastq/unzipped/*.trim.fastq Remolina/SRA/reads/gzip_fastq/unzipped/trimmed

#2) Change headers: read pair 1 end with /1 and pair 2 with /2 (otherwise error in SE2PE below)
for i in {1..6}
do
for j in {1..3}
do
sed 's/ length=//g' Cont_${j}_${i}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > Cont_${j}_${i}_r1.trim.headerEdit.fastq
sed 's/ length=//g' Cont_${j}_${i}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > Cont_${j}_${i}_r2.trim.headerEdit.fastq
sed 's/ length=//g' Sel_${j}_${i}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > Sel_${j}_${i}_r1.trim.headerEdit.fastq
sed 's/ length=//g' Sel_${j}_${i}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > Sel_${j}_${i}_r2.trim.headerEdit.fastq
done
done

#3) Map using bwa
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir mapped
for file in *.headerEdit.fastq
do
bwa bwasw -t 18 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta $file > mapped/$file.sam
done

#4) Restore PE info
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed
for i in {1..6}
do
for j in {1..3}
do
java -Xmx5g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Cont_${j}_${i}_r1.trim.headerEdit.fastq --fastq2 Cont_${j}_${i}_r2.trim.headerEdit.fastq --bam1 mapped/Cont_${j}_${i}_r1.trim.headerEdit.fastq.sam --bam2 mapped/Cont_${j}_${i}_r2.trim.headerEdit.fastq.sam --sort --output mapped/Cont_${j}_${i}.sort.bam
java -Xmx5g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 Sel_${j}_${i}_r1.trim.headerEdit.fastq --fastq2 Sel_${j}_${i}_r2.trim.headerEdit.fastq --bam1 mapped/Sel_${j}_${i}_r1.trim.headerEdit.fastq.sam --bam2 mapped/Sel_${j}_${i}_r2.trim.headerEdit.fastq.sam --sort --output mapped/Sel_${j}_${i}.sort.bam
done
done

#5) Merge multiple lanes into one file with samtools
cd Remolina/SRA/reads/gzip_fastq/unzipped/trimmed/mapped

samtools merge Cont_1_merged.bam Cont_1_1.sort.bam Cont_1_2.sort.bam Cont_1_3.sort.bam Cont_1_4.sort.bam Cont_1_5.sort.bam Cont_1_6.sort.bam
samtools merge Cont_2_merged.bam Cont_2_1.sort.bam Cont_2_2.sort.bam Cont_2_3.sort.bam Cont_2_4.sort.bam Cont_2_5.sort.bam Cont_2_6.sort.bam
samtools merge Cont_3_merged.bam Cont_3_1.sort.bam Cont_3_2.sort.bam Cont_3_3.sort.bam Cont_3_4.sort.bam Cont_3_5.sort.bam Cont_3_6.sort.bam
samtools merge Sel_1_merged.bam Sel_1_1.sort.bam Sel_1_2.sort.bam Sel_1_3.sort.bam Sel_1_4.sort.bam Sel_1_5.sort.bam Sel_1_6.sort.bam
samtools merge Sel_2_merged.bam Sel_2_1.sort.bam Sel_2_2.sort.bam Sel_2_3.sort.bam Sel_2_4.sort.bam Sel_2_5.sort.bam Sel_2_6.sort.bam
samtools merge Sel_3_merged.bam Sel_3_1.sort.bam Sel_3_2.sort.bam Sel_3_3.sort.bam Sel_3_4.sort.bam Sel_3_5.sort.bam Sel_3_6.sort.bam

#6) Generate the ppileup
mkdir Remolina/popTE2/
cd Remolina/reads/unzipped/trimmed/mapped/

java -Xmx8g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar ppileup --bam Cont_1_merged.bam --bam Cont_2_merged.bam --bam Cont_3_merged.bam --bam Sel_1_merged.bam --bam Sel_2_merged.bam --bam Sel_3_merged.bam --map-qual 15 --hier reference/te-hierarchy_edit.txt --output Remolina/popTE2/Cont_Sel_Remo.ppileup.gz

#7) Joint analysis to get TE positions and frequencies
cd Remolina/popTE2/
mkdir joint_analysis
java -Xmx4g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar identifySignatures --ppileup Cont_Sel_Remo.ppileup.gz --mode joint --output joint_analysis/Cont_Sel_Remo.ppileup_TEST.signatures --min-count 2 --detailed-log

#7) Get signatures - min count 3 (6 populations)
java -Xmx8g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar identifySignatures --ppileup Cont_Sel_Remo.ppileup.gz --mode joint --output joint_analysis/Cont_Sel_Remo.ppileup.signatures --min-count 3 --signature-window fix500 --min-valley fix150 --detailed-log
java -Xmx5g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar frequency --ppileup Cont_Sel_Remo.ppileup.gz --signature joint_analysis/Cont_Sel_Remo.ppileup.signatures --output joint_analysis/Cont_Sel_Remo.freqsig
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature joint_analysis/Cont_Sel_Remo.freqsig --detailed-log --output-detail medium --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --min-distance -200 --max-distance 300 --output joint_analysis/Cont_Sel_Remo.teinsertions

#8) Filter low quality TE insertions
cd Remolina/popTE2/joint_analysis/
mkdir filter
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar filterSignatures --input Cont_Sel_Remo.freqsig --output filter/Cont_Sel_Remo.filter.freqsig --max-otherte-count 2 --max-structvar-count 2
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature filter/Cont_Sel_Remo.filter.freqsig --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --detailed-log --output-detail medium --min-distance -200 --max-distance 300 --output filter/Cont_Sel_Remo.filter.teinsertions

##################################################################################################################
#Hoedjes2019
##################################################################################################################

#1) Trim reads in paired end mode
cd Hoedjes/SRA/reads/gzip_fastq/unzipped
mkdir trimmed

for j in {1..4}
do
cutadapt -j 15 --minimum-length 50 -q 18 -o CE${j}_r1.trim.fastq -p CE${j}_r2.trim.fastq CE${j}_r1.fastq CE${j}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o CP${j}_r1.trim.fastq -p CP${j}_r2.trim.fastq CP${j}_r1.fastq CP${j}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o LE${j}_r1.trim.fastq -p LE${j}_r2.trim.fastq LE${j}_r1.fastq LE${j}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o LP${j}_r1.trim.fastq -p LP${j}_r2.trim.fastq LP${j}_r1.fastq LP${j}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o HE${j}_r1.trim.fastq -p HE${j}_r2.trim.fastq HE${j}_r1.fastq HE${j}_r2.fastq
cutadapt -j 15 --minimum-length 50 -q 18 -o HP${j}_r1.trim.fastq -p HP${j}_r2.trim.fastq HP${j}_r1.fastq HP${j}_r2.fastq
done

mv Hoedjes/SRA/reads/gzip_fastq/unzipped/*.trim.fastq Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed

#2) Change headers: read pair 1 end with /1 and pair 2 with /2 (otherwise error in SE2PE below)
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed
for j in {1..4}
do
sed 's/ length=//g' CE${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > CE${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' CE${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > CE${j}_r2.trim.headerEdit.fastq
sed 's/ length=//g' CP${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > CP${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' CP${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > CP${j}_r2.trim.headerEdit.fastq
sed 's/ length=//g' HE${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > HE${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' HE${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > HE${j}_r2.trim.headerEdit.fastq
sed 's/ length=//g' HP${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > HP${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' HP${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > HP${j}_r2.trim.headerEdit.fastq
sed 's/ length=//g' LE${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > LE${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' LE${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/1/g' > LE${j}_r2.trim.headerEdit.fastq
sed 's/ length=//g' LP${j}_r1.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > LP${j}_r1.trim.headerEdit.fastq
sed 's/ length=//g' LP${j}_r2.trim.fastq | sed 's/\s.*$//' | sed -e 's/\(^@SRR.*\)\./\1_/g' | sed -e 's/\(^@SRR.*\)/\1\/2/g' > LP${j}_r2.trim.headerEdit.fastq
done

#3) Map using bwa
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir mapped
for file in *.headerEdit.fastq
do
bwa bwasw -t 22 reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta $file > mapped/$file.sam
done

#4) Restore PE info
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed
mkdir Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/tmp
mkdir Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/tmp1
mkdir Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/tmp2

for j in {1..4}
do
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 CE${j}_r1.trim.headerEdit.fastq --fastq2 CE${j}_r2.trim.headerEdit.fastq --bam1 mapped/CE${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/CE${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/CE${j}.sort.bam
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 CP${j}_r1.trim.headerEdit.fastq --fastq2 CP${j}_r2.trim.headerEdit.fastq --bam1 mapped/CP${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/CP${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/CP${j}.sort.bam
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp1 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 LE${j}_r1.trim.headerEdit.fastq --fastq2 LE${j}_r2.trim.headerEdit.fastq --bam1 mapped/LE${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/LE${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/LE${j}.sort.bam
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp1 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 LP${j}_r1.trim.headerEdit.fastq --fastq2 LP${j}_r2.trim.headerEdit.fastq --bam1 mapped/LP${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/LP${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/LP${j}.sort.bam
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 HE${j}_r1.trim.headerEdit.fastq --fastq2 HE${j}_r2.trim.headerEdit.fastq --bam1 mapped/HE${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/HE${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/HE${j}.sort.bam
java -Xmx15g -Djava.io.tmpdir=`pwd`/tmp2 -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar se2pe --fastq1 HP${j}_r1.trim.headerEdit.fastq --fastq2 HP${j}_r2.trim.headerEdit.fastq --bam1 mapped/HP${j}_r1.trim.headerEdit.fastq.sam --bam2 mapped/HP${j}_r2.trim.headerEdit.fastq.sam --sort --output mapped/HP${j}.sort.bam
done

#5) Generate the ppileup
mkdir Hoedjes/popTE2/
cd Hoedjes/SRA/reads/gzip_fastq/unzipped/trimmed/mapped/
mkdir tmp

java -Djava.io.tmpdir=`pwd`/tmp -Xmx10g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar ppileup --bam CE1.sort.bam --bam CE2.sort.bam --bam CE3.sort.bam --bam CE4.sort.bam --bam CP1.sort.bam --bam CP2.sort.bam --bam CP3.sort.bam --bam CP4.sort.bam --bam LE1.sort.bam --bam LE2.sort.bam --bam LE3.sort.bam --bam LE4.sort.bam --bam LP1.sort.bam --bam LP2.sort.bam --bam LP3.sort.bam --bam LP4.sort.bam --bam HE1.sort.bam --bam HE2.sort.bam --bam HE3.sort.bam --bam HE4.sort.bam --bam HP1.sort.bam --bam HP2.sort.bam --bam HP3.sort.bam --bam HP4.sort.bam --map-qual 15 --hier reference/te-hierarchy_edit.txt --output Hoedjes/popTE2/Cont_Sel_Hoedjes.ppileup.gz

#6) Subset to Dmel only (to reduce load on memory)
cd Hoedjes/popTE2/
zcat Hoedjes/popTE2/Cont_Sel_Hoedjes.ppileup.gz | grep -wE '^@..|^2L|^2R|^3L|^3R|^X|^4' > Hoedjes/popTE2/Cont_Sel_Hoedjes.Dmel.ppileup
gzip Hoedjes/popTE2/Cont_Sel_Hoedjes.Dmel.ppileup

#7) Joint analysis to get TE positions and frequencies
mkdir Hoedjes/popTE2/joint_analysis

#8) Get signatures - min count 12 (present in at least 12 pop out of 24)
cd Hoedjes/popTE2/
mkdir tmp
java -Xmx8g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar identifySignatures --ppileup Cont_Sel_Hoedjes.Dmel.ppileup.gz --mode joint --output joint_analysis/Cont_Sel_Hoedjes.ppileup.signatures --min-count 12 --signature-window fix500 --min-valley fix150 --detailed-log
java -Xmx55g -Djava.io.tmpdir=`pwd`/tmp -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar frequency --ppileup Cont_Sel_Hoedjes.Dmel.ppileup.gz --signature joint_analysis/Cont_Sel_Hoedjes.ppileup.signatures --output joint_analysis/Cont_Sel_Hoedjes.freqsig
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature joint_analysis/Cont_Sel_Hoedjes.freqsig --detailed-log --output-detail medium --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --min-distance -200 --max-distance 300 --output joint_analysis/Cont_Sel_Hoedjes.teinsertions

#9) Filter low quality TE insertions
cd Hoedjes/popTE2/joint_analysis
mkdir filter
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar filterSignatures --input Cont_Sel_Hoedjes.freqsig --output filter/Cont_Sel_Hoedjes.filter.freqsig --max-otherte-count 2 --max-structvar-count 2
java -Xmx6g -jar /PATH_TO_POPTE2/popte2-v1.10.03.jar pairupSignatures --signature filter/Cont_Sel_Hoedjes.filter.freqsig --ref-genome reference/rep_masked/dmel-all-chromosome-r6.27.fasta.masked_with_TE.fasta --hier reference/te-hierarchy_edit.txt --detailed-log --output-detail medium --min-distance -200 --max-distance 300 --output filter/Cont_Sel_Hoedjes.filter.teinsertions



