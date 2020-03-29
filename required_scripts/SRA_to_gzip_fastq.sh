#TRANSFORM SRA FILES TO GZIP FASTQ
#We downloaded and extracted sratoolkit.2.9.4-centos_linux64.tar.gz from: ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.4/ 

#Please edit code below so that the full path to sratoolkit.2.9.4-centos_linux64/bin/fastq-dump is given (i.e. replace /EDIT_PATH/ with correct path).

#SRA file to fastq.gz:
mkdir $1/gzip_fastq
for file in $1/*
do
/EDIT_PATH/sratoolkit.2.9.4-centos_linux64/bin/fastq-dump $file --gzip --split-3 --outdir $1/gzip_fastq
done
echo "Writing output to $1/gzip_fastq"