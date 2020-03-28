#ESTIMATION OF TE ABUNDACNE WITH DEVIATE
#TE abundance was estimated using the tool deviaTE (Weilguny and Kofler, 2019, Mol.Ecol.Res.)
#Also see: https://github.com/W-L/deviaTE 

#If necessary, change folder names and edit commands accordingly.
#Most commands were run on a computer cluster with multiple cores.

#Custom scripts used here (see required_scripts folder):
rename_fastq_headers.sh

#Raw outputs from DeviaTE and edited files for analysis in R are available at dryad: https://doi.org/10.5061/dryad.s7h44j13r

#########################
#INSTALLATION OF DEVIATE
#########################

#1) Install Anaconda3: here we used (Anaconda3-5.0.1-Linux-x86_64.sh)
bash Anaconda3-5.0.1-Linux-x86_64.sh

#2) Install conda environment with deviaTE
conda create deviaTE -c r -c defaults -c conda-forge -c bioconda -c w-l -n deviaTE_env

#3) Add anaconda path and library to bash profile, for instance with nano
nano ~/.bashrc
export PATH="/homes/user/anaconda3/bin:$PATH"
export LD_LIBRARY_PATH=/homes/user/anaconda3/lib:/homes/user/anaconda3/lib64:$LD_LIBRARY_PATH


###########################################################################
#RUN DEVIATE TO OBTAIN TE ABUNDANCE FOR CARNES ET AL 2015
###########################################################################
#1) deviaTE is not requiring paired-end information. We therefore merged fastq file with forward/reverse reads for each sample:
mkdir Carnes/reads/gzip_fastq/unzipped/merged/

for i in {1..5}
do
cat Carnes/reads/gzip_fastq/unzipped/Sel_O${i}_*.fastq > Carnes/reads/gzip_fastq/unzipped/merged/Sel_O${i}_merge.fastq
cat Carnes/reads/gzip_fastq/unzipped/Cont_B${i}_*.fastq > Carnes/reads/gzip_fastq/unzipped/merged/Cont_B${i}_merge.fastq
done

#2) Rename headers (to ensure deviaTE runs without raising an error) and move edited files to edit folder. Then rename files to remove unnecessary parts in file name.
cd Carnes/reads/gzip_fastq/unzipped/merged/
mkdir edit

for file in Carnes/reads/gzip_fastq/unzipped/merged/*.fastq
do
sh required_scripts/rename_fastq_headers.sh $file
done

mv *.fastq.* Carnes/reads/gzip_fastq/unzipped/merged/edit

cd Carnes/reads/gzip_fastq/unzipped/merged/edit
rename .fastq. . *.fastq

#3) Run deviaTE for Carnes et al 2015:
cd Carnes/reads/gzip_fastq/unzipped/merged/edit
source activate deviaTE_env
for file in *.fastq
do
deviaTE --input_fq $file --read_type phred+33 --min_read_len 50 --quality_threshold 18 --hq_threshold 20 --min_alignment_len 30 --threads 20 --families ALL --single_copy_genes Dmel_rpl32,Dmel_RpII140,Dmel_Act5C,Dmel_piwi,Dmel_p53
done
source deactivate

#4) Move output files to different folders to have a clearer structure
mkdir Carnes/TE_maps
mkdir Carnes/TE_maps/figures_out
mkdir Carnes/TE_maps/res_files
mkdir Carnes/TE_maps/res_files/raw
mkdir Carnes/TE_maps/res_files/corrected
mkdir Carnes/TE_maps/res_files/corrected/edit

mv Carnes/reads/gzip_fastq/unzipped/merged/edit/*.fused.sort.* Carnes/TE_maps
mv Carnes/reads/gzip_fastq/unzipped/merged/edit/*.pdf Carnes/TE_maps/figures_out
mv Carnes/reads/gzip_fastq/unzipped/merged/edit/*.raw Carnes/TE_maps/res_files/raw
mv Carnes/reads/gzip_fastq/unzipped/merged/edit/*.fastq.* Carnes/TE_maps/res_files/corrected

#5) Edit output file: subset to interesting columns, change sample names, and remove single copy genes from table to avoid fitting them accidentaly into downstream analysis models
for file in Carnes/TE_maps/res_files/corrected/*.fastq.*
do
grep -v '^#' $file | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$9,$10,$11}' | sed 's/_merge.header_edit.fastq//g' | grep -v "^Dmel" >> Carnes/TE_maps/res_files/corrected/edit/carnes_all_TE_edited.txt
done

###########################################################################
#RUN DEVIATE TO OBTAIN TE ABUNDANCE FOR FABIAN ET AL 2018
###########################################################################

#1) deviaTE is not requiring paired-end information. We therefore merged fastq files from different lanes and forward/reverse reads for each sample:
mkdir Fabian/reads/unzipped/merge/
cat Fabian/reads/unzipped/Cont_Ra1.fastq Fabian/reads/unzipped/Cont_Ra2.fastq > Fabian/reads/unzipped/merge/Cont_Ra_merge.fastq
cat Fabian/reads/unzipped/Cont_Rb1.fastq Fabian/reads/unzipped/Cont_Rb2.fastq > Fabian/reads/unzipped/merge/Cont_Rb_merge.fastq
cat Fabian/reads/unzipped/Sel_La1.fastq Fabian/reads/unzipped/Sel_La2.fastq > Fabian/reads/unzipped/merge/Sel_La_merge.fastq
cat Fabian/reads/unzipped/Sel_Lb1.fastq Fabian/reads/unzipped/Sel_Lb2.fastq > Fabian/reads/unzipped/merge/Sel_Lb_merge.fastq
cat Fabian/reads/unzipped/Sel_2La1.fastq Fabian/reads/unzipped/Sel_2La2.fastq > Fabian/reads/unzipped/merge/Sel_2La_merge.fastq
cat Fabian/reads/unzipped/Sel_2Lb1.fastq Fabian/reads/unzipped/Sel_2Lb2.fastq > Fabian/reads/unzipped/merge/Sel_2Lb_merge.fastq

#2) Rename headers (to ensure deviaTE runs without raising an error) and move edited files to edit folder. Then rename files to remove unnecessary parts in file name.
cd Fabian/reads/unzipped/merge/
mkdir edit

for file in *.fastq
do
sh required_scripts/rename_fastq_headers.sh $file
done

#3) Run deviaTE for Fabian et al 2018:
#Please check quality encoding here. If downloaded from SRA, it was likely transformed to Sanger encoding (phred+64). Hence, the option --read-type possibly needs to be changed.
source activate deviaTE_env
for file in Fabian/reads/unzipped/merge/edit/*.fastq
do
deviaTE --input_fq $file --read_type phred+64 --min_read_len 50 --quality_threshold 18 --hq_threshold 20 --min_alignment_len 30 --threads 20 --families ALL --single_copy_genes Dmel_rpl32,Dmel_RpII140,Dmel_Act5C,Dmel_piwi,Dmel_p53
done
source activate deviaTE_env

#4) Move output files to different folders to have a clearer structure
#move files
mkdir Fabian/TE_maps
mkdir Fabian/TE_maps/figures_out
mkdir Fabian/TE_maps/res_files
mkdir Fabian/TE_maps/res_files/raw
mkdir Fabian/TE_maps/res_files/corrected
mkdir Fabian/TE_maps/res_files/corrected/edit

mv Fabian/reads/unzipped/merge/edit/*.fused.sort.* Fabian/TE_maps
mv Fabian/reads/unzipped/merge/edit/*.pdf Fabian/TE_maps/figures_out
mv Fabian/reads/unzipped/merge/edit/*.raw Fabian/TE_maps/res_files/raw
mv Fabian/reads/unzipped/merge/edit/*.fastq.* Fabian/TE_maps/res_files/corrected
mv Fabian/reads/unzipped/merge/edit/*.txt Fabian/TE_maps/res_files/

#5) Edit output file: subset to interesting columns, change sample names, and remove single copy genes from table to avoid fitting them accidentaly into downstream analysis models
for file in Fabian/TE_maps/res_files/corrected/*.fastq.*
do
grep -v '^#' $file | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$9,$10,$11}' | sed 's/Fabian\/unzipped\/merge\/edit\///g' | sed 's/_merge.header_edit.fastq//g' | grep -v "^Dmel" >> Fabian/TE_maps/res_files/corrected/edit/fabian_all_TE_edited.txt
done

###########################################################################
#RUN DEVIATE TO OBTAIN TE ABUNDANCE FOR REMOLINA ET AL 2012
###########################################################################

#1) deviaTE is not requiring paired-end information. We therefore merged fastq files from different lanes and forward/reverse reads for each sample:
mkdir Remolina/reads/gzip_fastq/unzipped/merged/

for i in {1..3}
do
cat Remolina/reads/gzip_fastq/unzipped/Cont_${i}_*.fastq > Remolina/reads/gzip_fastq/unzipped/merged/Cont_${i}_merge.fastq
cat Remolina/reads/gzip_fastq/unzipped/Sel_${i}_*.fastq > Remolina/reads/gzip_fastq/unzipped/merged/Sel_${i}_merge.fastq
done

#2) Rename headers (to ensure deviaTE runs without raising an error) and move edited files to edit folder. Then rename files to remove unnecessary parts in file name.
cd Remolina/reads/gzip_fastq/unzipped/merged/
mkdir edit

for file in Remolina/reads/gzip_fastq/unzipped/merged/*.fastq
do
sh required_scripts/rename_fastq_headers.sh $file
done

mv *.fastq.* Remolina/reads/gzip_fastq/unzipped/merged/edit
cd Remolina/reads/gzip_fastq/unzipped/merged/edit
rename .fastq. . *.fastq

#3) Run deviaTE for Remolina et al 2012:
source activate deviaTE_env
for file inRemolina/reads/gzip_fastq/unzipped/merged/edit/*.fastq
do
	deviaTE --input_fq $file --read_type phred+33 --min_read_len 50 --quality_threshold 18 --hq_threshold 20 --min_alignment_len 30 --threads 20 --families ALL --single_copy_genes Dmel_rpl32,Dmel_RpII140,Dmel_Act5C,Dmel_piwi,Dmel_p53
done
source deactivate

#4) Move output files to different folders to have a clearer structure
mkdir Remolina/TE_maps
mkdir Remolina/TE_maps/figures_out
mkdir Remolina/TE_maps/res_files
mkdir Remolina/TE_maps/res_files/raw
mkdir Remolina/TE_maps/res_files/corrected
mkdir Remolina/TE_maps/res_files/corrected/edit

mv Remolina/reads/gzip_fastq/unzipped/merged/edit/*.fused.sort.* Remolina/TE_maps
mv Remolina/reads/gzip_fastq/unzipped/merged/edit/*.pdf Remolina/TE_maps/figures_out
mv Remolina/reads/gzip_fastq/unzipped/merged/edit/*.raw Remolina/TE_maps/res_files/raw
mv Remolina/reads/gzip_fastq/unzipped/merged/edit/*.fastq.* Remolina/TE_maps/res_files/corrected

#5) Edit output file: subset to interesting columns, change sample names, and remove single copy genes from table to avoid fitting them accidentaly into downstream analysis models
for file in Remolina/TE_maps/res_files/corrected/*.fastq.*
do
grep -v '^#' $file | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$9,$10,$11}' | sed 's/Remolina\/reads\/gzip_fastq\/unzipped\/merged\/edit\///g' | sed 's/_merge.header_edit.fastq//g' | grep -v "^Dmel" >> Remolina/TE_maps/res_files/corrected/edit/remolina_all_TE_edited.txt
done

###########################################################################
#RUN DEVIATE TO OBTAIN TE ABUNDANCE FOR FABIAN ET AL 2018
###########################################################################

#1) deviaTE is not requiring paired-end information. We therefore merged fastq files from different lanes and forward/reverse reads for each sample:
mkdirHoedjes/reads/gzip_fastq/unzipped/merged/

for i in {1..4}
do
cat Hoedjes/reads/gzip_fastq/unzipped/CE${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/CE${i}_merge.fastq
cat Hoedjes/reads/gzip_fastq/unzipped/CP${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/CP${i}_merge.fastq
cat Hoedjes/reads/gzip_fastq/unzipped/HE${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/HE${i}_merge.fastq
cat Hoedjes/reads/gzip_fastq/unzipped/HP${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/HP${i}_merge.fastq
cat Hoedjes/reads/gzip_fastq/unzipped/LE${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/LE${i}_merge.fastq
cat Hoedjes/reads/gzip_fastq/unzipped/LP${i}_*.fastq > Hoedjes/reads/gzip_fastq/unzipped/merged/LP${i}_merge.fastq
don

#2) Rename headers (to ensure deviaTE runs without raising an error) and move edited files to edit folder. Then rename files to remove unnecessary parts in file name.
cd Hoedjes/reads/gzip_fastq/unzipped/merged/
mkdir edit

for file in Hoedjes/reads/gzip_fastq/unzipped/merged/*.fastq
do
sh required_scripts/rename_fastq_headers.sh $file
done

mv *.fastq.* Hoedjes/reads/gzip_fastq/unzipped/merged/edit
cd Hoedjes/reads/gzip_fastq/unzipped/merged/edit
rename .fastq. . *.fastq

#3) Run deviaTE for Hoedjes et al 2019:
source activate deviaTE_env
for file in Hoedjes/reads/gzip_fastq/unzipped/merged/edit/*.fastq
do
deviaTE --input_fq $file --read_type phred+33 --min_read_len 50 --quality_threshold 18 --hq_threshold 20 --min_alignment_len 30 --threads 25 --families ALL --single_copy_genes Dmel_rpl32,Dmel_RpII140,Dmel_Act5C,Dmel_piwi,Dmel_p53
done
source deactivate

#4) Move output files to different folders to have a clearer structure
#move files
mkdir Hoedjes/TE_maps
mkdir Hoedjes/TE_maps/figures_out
mkdir Hoedjes/TE_maps/res_files
mkdir Hoedjes/TE_maps/res_files/raw
mkdir Hoedjes/TE_maps/res_files/corrected
mkdir Hoedjes/TE_maps/res_files/corrected/edit

#move bam files
mv Hoedjes/reads/gzip_fastq/unzipped/merged/edit/*.fused.sort.* Hoedjes/TE_maps
mv Hoedjes/reads/gzip_fastq/unzipped/merged/edit/*.pdf Hoedjes/TE_maps/figures_out
mv Hoedjes/reads/gzip_fastq/unzipped/merged/edit/*.raw Hoedjes/TE_maps/res_files/raw
mv Hoedjes/reads/gzip_fastq/unzipped/merged/edit/*.fastq.* Hoedjes/TE_maps/res_files/corrected

#5) Edit output file: subset to interesting columns, change sample names, and remove single copy genes from table to avoid fitting them accidentaly into downstream analysis models
for file in Hoedjes/TE_maps/res_files/corrected*.fastq.*
do
grep -v '^#' $file | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$9,$10,$11}' | sed 's/Hoedjes\/reads\/gzip_fastq\/unzipped\/merged\/edit\///g' | sed 's/_merge.header_edit.fastq//g' | grep -v "^Dmel" >> Hoedjes/TE_maps/res_files/corrected/edit/hoedjes_all_TE_edited.txt
done
