#FASTQC ANALYSIS: 
#We downloaded and extracted fastqc_v0.11.8.zip from: http://www.bioinformatics.babraham.ac.uk/projects/download.html

#Please edit code below so that the full path to fastqc is given (i.e. replace /EDIT_PATH/ with correct path).

mkdir $1/fastqc_output

for file in $1/*
do
#echo "$file" >> $1/fastqc_output/test.txt
/EDIT_PATH/fastqc $file --outdir=$1/fastqc_output
done
echo "Output will be saved in $1/fastqc_output"