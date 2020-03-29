#rename fastq headers of reads, so that they are numbered from 1 to N.
#$1 corresponds to a fastq file as input.

awk '{print (NR%4 == 1) ? "@Read_" ++i : $0}' $1 > $1.header_edit.fastq
