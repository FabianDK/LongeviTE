#DOWNLOAD OF SRA FILES
#Here we explain how we downloaded raw sequence reads from the SRA

#If necessary, change folder names and edit commands accordingly.
#Most commands were run on a computer cluster with multiple cores and large available memory

#Custom scripts used here (see required_scripts folder):
SRA_to_gzip_fastq.sh #Requires fastq-dump from SRA tool kits, see SRA_to_gzip_fastq.sh script for more details.
fastq_check.sh
unzip_many_files.sh

########################################
#CARNES ET AL (2015) - RAW READS
########################################
mkdir Carnes
mkdir Carnes/reads

#1) Get run info from BioProject
cd Carnes
wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA286855' -O - | tee SraRunInfo.csv

#2) Loop through lines of file, extract download link and use wget to download it to folder "Carnes/reads"
for runinfo in SraRunInfo.csv
do
x=$(sed 's/,/\t/g' $runinfo | cut -f10 | grep "http")
wget $x -P Carnes/reads
echo "$x"
done

#3) Transform .sra files into gzip fastq files using fastq-dump
sh required_scripts/SRA_to_gzip_fastq.sh Carnes/reads

#4) Use FASTQC to perform quality checks on fastq files
sh required_scripts/fastq_check.sh Carnes/reads/gzip_fastq

#5) Unzip .gz files and move to unzipped folder
sh required_scripts/unzip_many_files.sh Carnes/reads/gzip_fastq

#6) Rename files and remove .gz.decomp extension
for filename in Carnes/reads/gzip_fastq/unzipped/*.decomp
do 
[ -f "$filename" ] || continue
mv "$filename" "${filename%.gz.decomp}"
done

#7) Created table with two columns (SRA accession number, and name of population) based onn SraRunInfo.csv: carnes_rename_tab.txt. 
#Then rename file names so that they are labelled according to selection regime and replicate
cd Carnes/reads/gzip_fastq/unzipped/
sed 's/^/mv -vi "/;s/\t/_1.fastq" "/;s/$/_r1.fastq";/' < extra_files/carnes_rename_tab.txt | bash -
sed 's/^/mv -vi "/;s/\t/_2.fastq" "/;s/$/_r2.fastq";/' < extra_files/carnes_rename_tab.txt | bash -
sed 's/^/mv -vi "/;s/\t/.fastq" "/;s/$/_noPair.fastq";/' < extra_files/carnes_rename_tab.txt | bash -

########################################
#REMOLINA ET AL (2012) - RAW READS
########################################
mkdir Remolina
mkdir Remolina/reads

#1) Get run info from BioProject
cd Remolina
wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA185744' -O - | tee SraRunInfo.csv

#2) Loop through lines of file, extract download link and use wget to download it to folder "Remolina/reads"
for runinfo in SraRunInfo.csv
do
x=$(sed 's/,/\t/g' $runinfo | cut -f10 | grep "http")
wget $x -P Remolina/reads
done

#3) Transform .sra files into gzip fastq files using fastq-dump
sh required_scripts/SRA_to_gzip_fastq.sh Remolina/reads

#4) Use FASTQC to perform quality checks on fastq files
sh required_scripts/fastq_check.sh Remolina/reads/gzip_fastq

#5) Unzip .gz files and move to unzipped folder
sh required_scripts/unzip_many_files.sh Remolina/reads/gzip_fastq

#6) Rename files and remove .gz.decomp extension
for filename in Remolina/reads/gzip_fastq/unzipped/*.decomp
do 
[ -f "$filename" ] || continue
mv "$filename" "${filename%.gz.decomp}"
done

#7) Created table with two columns (SRA accession number, and name of population) based onn SraRunInfo.csv: remolina_rename_tab.txt. 
#Then rename file names so that they are labelled according to selection regime and replicate
cd Remolina/reads/gzip_fastq/unzipped/
sed 's/^/mv -vi "/;s/\t/_1.fastq" "/;s/$/_r1.fastq";/' < extra_files/remolina_rename_tab.txt | bash -
sed 's/^/mv -vi "/;s/\t/_2.fastq" "/;s/$/_r2.fastq";/' < extra_files/remolina_rename_tab.txt | bash -


########################################
#HOEDJES ET AL (2019) - RAW DNA READS
########################################
mkdir Hoedjes
mkdir Hoedjes/reads

#1) Get run info from BioProject
cd Hoedjes
wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA564570' -O - | tee SraRunInfo.csv

#2) Loop through lines of file, extract download link and use wget to download it to folder "Hoedjes/reads"
for runinfo in SraRunInfo.csv
do
x=$(sed 's/,/\t/g' $runinfo | cut -f10 | grep "http")
wget $x -P Hoedjes/reads
done

#3) Transform .sra files into gzip fastq files using fastq-dump
sh required_scripts/SRA_to_gzip_fastq.sh Hoedjes/reads

#4) Use FASTQC to perform quality checks on fastq files
sh required_scripts/fastq_check.sh Hoedjes/reads/gzip_fastq

#5) Unzip .gz files and move to unzipped folder
sh required_scripts/unzip_many_files.sh Hoedjes/reads/gzip_fastq

#6) Rename files and remove .gz.decomp extension
for filename in Hoedjes/reads/gzip_fastq/unzipped/*.decomp
do 
[ -f "$filename" ] || continue
mv "$filename" "${filename%.gz.decomp}"
done

#7) Created table with two columns (SRA accession number, and name of population) based onn SraRunInfo.csv: hoedjes_rename_tab.txt. 
#Then rename file names so that they are labelled according to selection regime and replicate
cd Hoedjes/reads/gzip_fastq/unzipped/
sed 's/^/mv -vi "/;s/\t/_1.fastq" "/;s/$/_r1.fastq";/' < extra_files/hoedjes_rename_tab.txt | bash -
sed 's/^/mv -vi "/;s/\t/_2.fastq" "/;s/$/_r2.fastq";/' < extra_files/hoedjes_rename_tab.txt | bash -


########################################
#FABIAN ET AL (2018) - RAW READS
########################################
#We have not downloaded read files from Fabian et al. from the repository, but instead used our own local copies.
#Raw fastq reads used in this study have been deposited to the ENA under the study accession (PRJEB28048): https://www.ebi.ac.uk/ena/data/view/PRJEB28048
#SRA files could be downloaded and extracted similar as for the studies above (points 1 to 3). The BioProject repository includes bam files from an earlier study and raw fastq reads. For the analysis here, we only used the raw fastq reads.

mkdir Fabian
mkdir Fabian/reads

#1) Example of download from ENA, fastqc check and decompressing files: 
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764144/Cont_Ra1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764144/Cont_Ra2.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764145/Cont_Rb1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764145/Cont_Rb2.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764146/Sel_La1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764146/Sel_La2.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764147/Sel_Lb1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764147/Sel_Lb2.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764148/Sel_2La1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764148/Sel_2La2.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764149/Sel_2Lb1.fastq.gz -P Fabian/reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR376/ERR3764149/Sel_2Lb2.fastq.gz -P Fabian/reads

sh required_scripts/fastq_check.sh Fabian/reads/gzip_fastq
sh required_scripts/unzip_many_files.sh Fabian/reads/gzip_fastq

#2) Rename files and remove .gz.decomp extension
for filename in Fabian/reads/unzipped/*.decomp
do 
[ -f "$filename" ] || continue
mv "$filename" "${filename%.gz.decomp}"
done
