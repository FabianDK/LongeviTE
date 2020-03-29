#Example on how to run SNPeff

#1) Install SNPEff: http://snpeff.sourceforge.net/download.html
#We used: snpEff version SnpEff 4.0e

#2) Create database for genome annotation (please also see http://snpeff.sourceforge.net/SnpEff_manual.html#databases)
mkdir /PATH/snpEff/dm6.27
cd /PATH/snpEff/dm6.27
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.27_FB2019_02/gtf/dmel-all-r6.27.gtf.gz
mv dmel-all-r6.27.gtf.gz genes.gtf.gz

cd /PATH/snpEff/data/genomes
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.27.fasta.gz
mv dmel-all-chromosome-r6.27.fasta.gz dm6.27.fa.gz

#Add the new genome to the config file (see section "How to add a new genome to the configuration file" for details)
cd /PATH/snpEff/
java -Xmx15g -jar snpEff.jar build -gff3 -v dm6.27

#3) Edit frequency files so that they are vcf-like:
# Chrom Pos     V4      TEfam   Class   C1      C2      Eff     C3      S1      S2      S3      AvFreq  AvFreqC AvFreqS dFreq   Pval    FDR
2R      4194623 .       GYPSY7  LTR     0.887   0.899           0.91    0.928   0.916   0.879   0.903166667     0.898666667     0.907666667     0.009   0.574970169     0.996913484
2R      653171  .       GYPSY7  LTR     0.9     0.775           0.834   0.9     0.875   0.691   0.829166667     0.836333333     0.822   -0.014333333    0.90534451      0.997040011
2R      657831  .       GYPSY7  LTR     0.123   0.135           0.1     0.05    0.202   0.119   0.1215  0.119333333     0.123666667     0.004333333     0.956341079     0.998789389
2R      1485659 .       GYPSY6  LTR     0.942   0.832           0.99    0.95    0.932   0.871   0.9195  0.921333333     0.917666667     -0.003666667    0.781756729     0.997040011
#Header row needs a hashtag.
#Column 8 needs to be empty. This is where annotations are pasted.

#4) Run SNPeff with -ud 1000
java -Xmx15g -jar /PATH/snpEff/snpEff.jar dm6.27 edited_frequency_file.txt -classic -ud 1000 -s edited_frequency_file_1kb_summary.html > edited_frequency_file_1kb.snpeff

#5) Go to R and use our custom SNPeff function to edit .snpeff output file (see file "R_functions.R" for all custom functions)
snpeff.edit(edited_frequency_file_1kb.snpeff)
