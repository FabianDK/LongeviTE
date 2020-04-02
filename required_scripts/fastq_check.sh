#FASTQC ANALYSIS (on gzip compressed fastq files): 
#We downloaded and extracted fastqc_v0.11.8.zip from: http://www.bioinformatics.babraham.ac.uk/projects/download.html

#Please edit code below so that the full path to fastqc is given (i.e. replace /EDIT_PATH/ with correct path).

mkdir $1/fastqc_output

#Loop through *.gz files of a folder to create FASTQC summaries
for file in $1/*.gz
do
/EDIT_PATH/fastqc $file --outdir=$1/fastqc_output
done