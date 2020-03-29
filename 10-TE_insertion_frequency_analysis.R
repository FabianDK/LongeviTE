#TE frequency analysis using output from PoPoolationTE2

#For Figure 2, Figure 3BC, Table S9 to S11, Figure S7 to S9

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(mltools)
library(pBrackets)
library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/comp_all/PopTE2/")

dir.create("/Users/danfab/comp_all/PopTE2",recursive = T)
output_path = "/Users/danfab/comp_all/PopTE2/"

#Script does:
#1) Filtering of PopoolationTE2 output
#2) Calculation of average frequencies for each TE family
#3) Calculation of sum of frequencies over all detected insertions for each TE family as measure of abundance (reviewer suggestion)
#4) Correlation between TE abundance (from DeviaTE) and average frequency from PoPTE2
#5) Correlation with average frequencies from Kofler et al 2015 (PLoS Genetics) with plot
#6) Comparison of average frequencies of S>C and C>S using SA population as proxy (with plots)
#7) Candidate TE insertions: stats with manhattan plots
#8) Check if TE family frequency varies between S and C populations, and comparison to TE abundance differences

# Kofler et al 2015, PLoS Genetics: TableS3 with average frequencies (table was edited before analysis)
robert.2015 = read.table("/Users/danfab/extra_files/Robert_Plos2015_AvFreq_edit.txt",header=T)
head(robert.2015)
# Name    TEfam    av_freq Numeric Color1 Color2     Type plotpch
# 1 1360     1360 0.25996450      68   blue  inter Germline      21
# 2 17.6  DMIS176 0.07891631      44   blue   blue Germline      21
# 3 1731 DMTN1731 0.23047342      66   blue  inter Germline      21
# 4  297  DMIS297 0.08123330      47   blue   blue Germline      21
# 5 3S18  DM23420 0.07026803      40   blue   blue Germline      21
# 6  412      412 0.04022767       9   blue   blue     Soma      21

min(robert.2015$av_freq, na.rm=T) #0.02894331
max(robert.2015$av_freq, na.rm=T) #1

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)

#Tables with fitlered for min_coverage and only Bonferroni significant
remo.stat.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE",]
carnes.stat.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE",]
fabian.stat.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE",]
hoed.stat.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == "TRUE",]

#Group into S>C and C>S TEs
carnes.SC = carnes.stat.sign[carnes.stat.sign$Diff_SelCont > 0,]
carnes.CS = carnes.stat.sign[carnes.stat.sign$Diff_SelCont < 0,]
fabian.SC = fabian.stat.sign[fabian.stat.sign$Diff_SelCont > 0,]
fabian.CS = fabian.stat.sign[fabian.stat.sign$Diff_SelCont < 0,]
remo.SC = remo.stat.sign[remo.stat.sign$Diff_SelCont > 0,]
remo.CS = remo.stat.sign[remo.stat.sign$Diff_SelCont < 0,]
hoed.SC = hoed.stat.sign[hoed.stat.sign$Diff_PostEarly > 0,]
hoed.CS = hoed.stat.sign[hoed.stat.sign$Diff_PostEarly < 0,]

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey.dmel = ensembl_casey[ensembl_casey$species == "Dmel",]
ensembl_casey.dmel$flybase_name = gsub("Dmel/","",ensembl_casey.dmel$flybase_name,fixed = T)

#PoPoolationTE2 output
carnes.te2 = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/Cont_Sel.filter.teinsertions")
fabian.te2 = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/Cont_Sel_Fabian.filter.teinsertions")
remo.te2 = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/Cont_Sel_Remo.filter.teinsertions")
hoed.te2 = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/Cont_Sel_Hoedjes.Dmel.filter.teinsertions")

#Add column names, same order as in pileup
names(carnes.te2)[c(2:3,5:7,9:18)] = c("Chrom","Pos","TEfam","Class","Strand",
                                   "C1","C2","C3","C4","C5",
                                   "S1","S2","S3","S4","S5")
names(fabian.te2)[c(2:3,5:7,9:14)] = c("Chrom","Pos","TEfam","Class","Strand",
                                       "Ra","Rb","La","Lb","2La","2Lb")
names(remo.te2)[c(2:3,5:7,9:14)] = c("Chrom","Pos","TEfam","Class","Strand",
                                       "C1","C2","C3","S1","S2","S3")
names(hoed.te2)[c(2:3,5:7,9:32)] = c("Chrom","Pos","TEfam","Class","Strand",
                                                 "CE1","CE2","CE3","CE4",
                                                 "CP1","CP2","CP3","CP4",
                                                 "LE1","LE2","LE3","LE4",
                                                 "LP1","LP2","LP3","LP4",
                                                 "HE1","HE2","HE3","HE4",
                                                 "HP1","HP2","HP3","HP4")

#Create selection vectors for Hoedjes
hoed.cont.vect = c("CE1","CE2","CE3","CE4","LE1","LE2","LE3","LE4","HE1","HE2","HE3","HE4")
hoed.sel.vect = c("CP1","CP2","CP3","CP4","LP1","LP2","LP3","LP4","HP1","HP2","HP3","HP4")

hoed.dietC.vect = c("CE1","CE2","CE3","CE4", "CP1","CP2","CP3","CP4")
hoed.dietH.vect = c("HE1","HE2","HE3","HE4", "HP1","HP2","HP3","HP4")
hoed.dietL.vect = c("LE1","LE2","LE3","LE4","LP1","LP2","LP3","LP4")

hoed.dietCE.vect = c("CE1","CE2","CE3","CE4")
hoed.dietCP.vect = c("CP1","CP2","CP3","CP4")

hoed.dietHE.vect = c("HE1","HE2","HE3","HE4")
hoed.dietHP.vect = c("HP1","HP2","HP3","HP4")

hoed.dietLE.vect = c("LE1","LE2","LE3","LE4")
hoed.dietLP.vect = c("LP1","LP2","LP3","LP4")

nrow(carnes.te2) #783
nrow(fabian.te2) #8764
nrow(remo.te2) #4610
nrow(hoed.te2) #13397

#filter by major chromosomal arms
carnes.te2 = carnes.te2[carnes.te2$Chrom %in% c("2L","2R","3L","3R","X","4"),]
fabian.te2 = fabian.te2[fabian.te2$Chrom %in% c("2L","2R","3L","3R","X","4"),]
remo.te2 = remo.te2[remo.te2$Chrom %in% c("2L","2R","3L","3R","X","4"),]
hoed.te2 = hoed.te2[hoed.te2$Chrom %in% c("2L","2R","3L","3R","X","4"),]

nrow(carnes.te2) #742
nrow(fabian.te2) #8673
nrow(remo.te2) #4573
nrow(hoed.te2) #13397

#############
#AVERAGE FEQUENCY OF TEs
#############

#average frequency of each insertion across all pops
carnes.te2$AvFreq = rowMeans(carnes.te2[,c(9:18)],na.rm=T)
fabian.te2$AvFreq = rowMeans(fabian.te2[,c(9:14)],na.rm=T)
remo.te2$AvFreq = rowMeans(remo.te2[,c(9:14)],na.rm=T)
hoed.te2$AvFreq = rowMeans(hoed.te2[,c(9:32)],na.rm=T)

#average frequency of each insertion in controls and selected
carnes.te2$AvFreqC = rowMeans(carnes.te2[,c(9:13)],na.rm=T)
carnes.te2$AvFreqS = rowMeans(carnes.te2[,c(14:18)],na.rm=T)

fabian.te2$AvFreqC = rowMeans(fabian.te2[,c(9:10)],na.rm=T)
fabian.te2$AvFreqS = rowMeans(fabian.te2[,c(11:14)],na.rm=T)

remo.te2$AvFreqC = rowMeans(remo.te2[,c(9:11)],na.rm=T)
remo.te2$AvFreqS = rowMeans(remo.te2[,c(12:14)],na.rm=T)

hoed.te2$AvFreqC = rowMeans(hoed.te2[,hoed.cont.vect],na.rm=T)
hoed.te2$AvFreqS = rowMeans(hoed.te2[,hoed.sel.vect],na.rm=T)
hoed.te2$AvFreqCdiet = rowMeans(hoed.te2[,hoed.dietC.vect],na.rm=T)
hoed.te2$AvFreqLdiet = rowMeans(hoed.te2[,hoed.dietL.vect],na.rm=T)
hoed.te2$AvFreqHdiet = rowMeans(hoed.te2[,hoed.dietH.vect],na.rm=T)
cor.test(hoed.te2$AvFreqC,hoed.te2$AvFreqS) #highly correlated

#Delta Allele Frequency
carnes.te2$dFreq = carnes.te2$AvFreqS - carnes.te2$AvFreqC
fabian.te2$dFreq = fabian.te2$AvFreqS - fabian.te2$AvFreqC
remo.te2$dFreq = remo.te2$AvFreqS - remo.te2$AvFreqC
hoed.te2$dFreq = hoed.te2$AvFreqS - hoed.te2$AvFreqC

#count NaN per row - only TE frequencies taken were estimate for all populations
table(rowSums(carnes.te2 == "NaN")) #570 with 0, 110 with 1 NaN
table(rowSums(fabian.te2 == "NaN")) #8549
table(rowSums(remo.te2 == "NaN")) #4567
table(rowSums(hoed.te2 == "NaN")) #13215

#remove all rows with NaN
carnes.te2.filt = carnes.te2[rowSums(carnes.te2 == "NaN") == 0,]
nrow(carnes.te2.filt) #570
fabian.te2.filt = fabian.te2[rowSums(fabian.te2 == "NaN") == 0,]
nrow(fabian.te2.filt) #8549
remo.te2.filt = remo.te2[rowSums(remo.te2 == "NaN") == 0,]
nrow(remo.te2.filt) #4567
hoed.te2.filt = hoed.te2[rowSums(hoed.te2 == "NaN") == 0,]
nrow(hoed.te2.filt) #13215

luniq(droplevels(carnes.te2.filt$TEfam)) #93
luniq(carnes.te2$TEfam) #98
luniq(droplevels(fabian.te2.filt$TEfam)) #139
luniq(fabian.te2.filt$TEfam) #139
luniq(droplevels(remo.te2.filt$TEfam)) #130
luniq(remo.te2.filt$TEfam) #120
luniq(droplevels(hoed.te2.filt$TEfam)) #137
luniq(hoed.te2.filt$TEfam) #137

#Remove TEs not thought to occur in Dmel
carnes.te2.filt = carnes.te2.filt[!carnes.te2.filt$TEfam %in% ensembl_casey[!ensembl_casey$species == "Dmel",]$TEfam,]
fabian.te2.filt = fabian.te2.filt[!fabian.te2.filt$TEfam %in% ensembl_casey[!ensembl_casey$species == "Dmel",]$TEfam,]
remo.te2.filt = remo.te2.filt[!remo.te2.filt$TEfam %in% ensembl_casey[!ensembl_casey$species == "Dmel",]$TEfam,]
hoed.te2.filt = hoed.te2.filt[!hoed.te2.filt$TEfam %in% ensembl_casey[!ensembl_casey$species == "Dmel",]$TEfam,]

#Count of TE insertions reported in paper
nrow(carnes.te2.filt) #567
nrow(fabian.te2.filt) #8402
nrow(remo.te2.filt) #4502
nrow(hoed.te2.filt) #13018

#Group TE frequency file into significant C>S and S>C
carnes.te2.filt.CS = carnes.te2.filt[carnes.te2.filt$TEfam %in% carnes.CS$TEfam,]
fabian.te2.filt.CS = fabian.te2.filt[fabian.te2.filt$TEfam %in% fabian.CS$TEfam,]
remo.te2.filt.CS = remo.te2.filt[remo.te2.filt$TEfam %in% remo.CS$TEfam,]
hoed.te2.filt.CS = hoed.te2.filt[hoed.te2.filt$TEfam %in% hoed.CS$TEfam,]

carnes.te2.filt.SC = carnes.te2.filt[carnes.te2.filt$TEfam %in% carnes.SC$TEfam,]
fabian.te2.filt.SC = fabian.te2.filt[fabian.te2.filt$TEfam %in% fabian.SC$TEfam,]
remo.te2.filt.SC = remo.te2.filt[remo.te2.filt$TEfam %in% remo.SC$TEfam,]
hoed.te2.filt.SC = hoed.te2.filt[hoed.te2.filt$TEfam %in% hoed.SC$TEfam,]

carnes.te2.filt.CS$Type = "C>S"
carnes.te2.filt.SC$Type = "S>C"
carnes.te2.filt.CS.SC = rbind(carnes.te2.filt.CS,carnes.te2.filt.SC)

fabian.te2.filt.CS$Type = "C>S"
fabian.te2.filt.SC$Type = "S>C"
fabian.te2.filt.CS.SC = rbind(fabian.te2.filt.CS,fabian.te2.filt.SC)

remo.te2.filt.CS$Type = "C>S"
remo.te2.filt.SC$Type = "S>C"
remo.te2.filt.CS.SC = rbind(remo.te2.filt.CS,remo.te2.filt.SC)

hoed.te2.filt.CS$Type = "C>S"
hoed.te2.filt.SC$Type = "S>C"
hoed.te2.filt.CS.SC = rbind(hoed.te2.filt.CS,hoed.te2.filt.SC)

write.table(carnes.te2.filt,"TE_insertions_Carnes2015_filtered.txt",quote=F,row.names = F,sep="\t")
write.table(remo.te2.filt,"TE_insertions_Remolina2012_filtered.txt",quote=F,row.names = F,sep="\t")
write.table(fabian.te2.filt,"TE_insertions_Fabian2018_filtered.txt",quote=F,row.names = F,sep="\t")
write.table(hoed.te2.filt,"TE_insertions_Hoedjes2019_filtered.txt",quote=F,row.names = F,sep="\t")

#Average frequency per TE family per Regime with number of insertion counts
carnes.te2.stat = as.data.frame(cbind( count = table(carnes.te2.filt$TEfam),
                         FreqC = tapply(carnes.te2.filt$AvFreqC,carnes.te2.filt$TEfam,mean),
                         FreqS =tapply(carnes.te2.filt$AvFreqS,carnes.te2.filt$TEfam,mean),
                         FreqAv = tapply(carnes.te2.filt$AvFreq,carnes.te2.filt$TEfam,mean)))
carnes.te2.stat$dFreq = carnes.te2.stat$FreqS - carnes.te2.stat$FreqC

#looks like mainly large frequency TEs! 
minmax(na.omit(carnes.te2.stat$FreqAv)) #0.1932 0.8958

fabian.te2.stat = as.data.frame(cbind( count = table(fabian.te2.filt$TEfam),
                         FreqC = tapply(fabian.te2.filt$AvFreqC,fabian.te2.filt$TEfam,mean),
                         FreqS =tapply(fabian.te2.filt$AvFreqS,fabian.te2.filt$TEfam,mean),
                         FreqAv = tapply(fabian.te2.filt$AvFreq,fabian.te2.filt$TEfam,mean)))
fabian.te2.stat$dFreq = fabian.te2.stat$FreqS - fabian.te2.stat$FreqC
fabian.te2.stat
minmax(na.omit(fabian.te2.stat$FreqAv)) #0.023 1.000 | excluding nonDmel: 0.023 1.000 

remo.te2.stat = as.data.frame(cbind( count = table(remo.te2.filt$TEfam),
                       FreqC = tapply(remo.te2.filt$AvFreqC,remo.te2.filt$TEfam,mean),
                       FreqS =tapply(remo.te2.filt$AvFreqS,remo.te2.filt$TEfam,mean),
                       FreqAv = tapply(remo.te2.filt$AvFreq,remo.te2.filt$TEfam,mean)))
remo.te2.stat$dFreq = remo.te2.stat$FreqS - remo.te2.stat$FreqC
remo.te2.stat
minmax(na.omit(remo.te2.stat$FreqAv)) #0.03733333 0.91033333

hoed.te2.stat = as.data.frame(cbind( count = table(hoed.te2.filt$TEfam),
                                       FreqC = tapply(hoed.te2.filt$AvFreqC,hoed.te2.filt$TEfam,mean),
                                       FreqS =tapply(hoed.te2.filt$AvFreqS,hoed.te2.filt$TEfam,mean),
                                       FreqAv = tapply(hoed.te2.filt$AvFreq,hoed.te2.filt$TEfam,mean)))
hoed.te2.stat$dFreq = hoed.te2.stat$FreqS - hoed.te2.stat$FreqC
hoed.te2.stat
minmax(na.omit(hoed.te2.stat$FreqAv)) #0.009458333 0.894680556


####################
#COUNTS AND AVERAGE TE FAMILY FREQUENCY PER POPULATION
####################

#Carnes
carnes.te2.stat.pop = as.data.frame(cbind( count = table(carnes.te2.filt$TEfam),
                                           CountC1 = tapply(carnes.te2.filt$C1,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountC2 = tapply(carnes.te2.filt$C2,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountC3 = tapply(carnes.te2.filt$C3,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountC4 = tapply(carnes.te2.filt$C4,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountC5 = tapply(carnes.te2.filt$C5,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountS1 = tapply(carnes.te2.filt$S1,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountS2 = tapply(carnes.te2.filt$S2,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountS3 = tapply(carnes.te2.filt$S3,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountS4 = tapply(carnes.te2.filt$S4,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           CountS5 = tapply(carnes.te2.filt$S5,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                           FreqC1 = tapply(carnes.te2.filt$C1,carnes.te2.filt$TEfam,mean),
                                           FreqC2 = tapply(carnes.te2.filt$C2,carnes.te2.filt$TEfam,mean),
                                           FreqC3 = tapply(carnes.te2.filt$C3,carnes.te2.filt$TEfam,mean),
                                           FreqC4 = tapply(carnes.te2.filt$C4,carnes.te2.filt$TEfam,mean),
                                           FreqC5 = tapply(carnes.te2.filt$C5,carnes.te2.filt$TEfam,mean),
                                           FreqS1 = tapply(carnes.te2.filt$S1,carnes.te2.filt$TEfam,mean),
                                           FreqS2 = tapply(carnes.te2.filt$S2,carnes.te2.filt$TEfam,mean),
                                           FreqS3 = tapply(carnes.te2.filt$S3,carnes.te2.filt$TEfam,mean),
                                           FreqS4 = tapply(carnes.te2.filt$S4,carnes.te2.filt$TEfam,mean),
                                           FreqS5 = tapply(carnes.te2.filt$S5,carnes.te2.filt$TEfam,mean),
                                           FreqC =tapply(carnes.te2.filt$AvFreqC,carnes.te2.filt$TEfam,mean),
                                           FreqS =tapply(carnes.te2.filt$AvFreqS,carnes.te2.filt$TEfam,mean),
                                           FreqAv = tapply(carnes.te2.filt$AvFreq,carnes.te2.filt$TEfam,mean)))

carnes.te2.stat.pop$CountC = rowMeans(carnes.te2.stat.pop[,c(2:6)])
carnes.te2.stat.pop$CountS = rowMeans(carnes.te2.stat.pop[,c(7:11)])

carnes.te2.stat.pop$TEfam = row.names(carnes.te2.stat.pop)
carnes.te2.stat.pop$dFreq = carnes.te2.stat.pop$FreqS - carnes.te2.stat.pop$FreqC
head(carnes.te2.stat.pop)
carnes.te2.stat.pop = carnes.te2.stat.pop[carnes.te2.stat.pop$count > 0,]
# write.table(carnes.te2.stat.pop,"/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/carnes_average_freq_counts.txt",quote=F,row.names = F,sep="\t")

#Fabian
fabian.te2.stat.pop = as.data.frame(cbind( count = table(fabian.te2.filt$TEfam),
                                       CountRa = tapply(fabian.te2.filt$Ra,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       CountRb = tapply(fabian.te2.filt$Rb,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       CountLa = tapply(fabian.te2.filt$La,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       CountLb = tapply(fabian.te2.filt$Lb,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       Count2La = tapply(fabian.te2.filt$`2La`,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       Count2Lb = tapply(fabian.te2.filt$`2Lb`,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                       FreqRa = tapply(fabian.te2.filt$Ra,fabian.te2.filt$TEfam,mean),
                                       FreqRb = tapply(fabian.te2.filt$Rb,fabian.te2.filt$TEfam,mean),
                                       FreqLa = tapply(fabian.te2.filt$La,fabian.te2.filt$TEfam,mean),
                                       FreqLb = tapply(fabian.te2.filt$Lb,fabian.te2.filt$TEfam,mean),
                                       Freq2La = tapply(fabian.te2.filt$`2La`,fabian.te2.filt$TEfam,mean),
                                       Freq2Lb = tapply(fabian.te2.filt$`2Lb`,fabian.te2.filt$TEfam,mean),
                                       FreqC = tapply(fabian.te2.filt$AvFreqC,fabian.te2.filt$TEfam,mean),
                                       FreqS =tapply(fabian.te2.filt$AvFreqS,fabian.te2.filt$TEfam,mean),
                                       FreqAv = tapply(fabian.te2.filt$AvFreq,fabian.te2.filt$TEfam,mean)))

fabian.te2.stat.pop$CountC = rowMeans(fabian.te2.stat.pop[,c(2:3)])
fabian.te2.stat.pop$CountS = rowMeans(fabian.te2.stat.pop[,c(4:7)])

fabian.te2.stat.pop$TEfam = row.names(fabian.te2.stat.pop)
fabian.te2.stat.pop$dFreq = fabian.te2.stat.pop$FreqS - fabian.te2.stat.pop$FreqC
head(fabian.te2.stat.pop)
fabian.te2.stat.pop = fabian.te2.stat.pop[fabian.te2.stat.pop$count > 0,]
# write.table(fabian.te2.stat.pop,"/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/fabian_average_freq_counts.txt",quote=F,row.names = F,sep="\t")

#Remolina
remo.te2.stat.pop = as.data.frame(cbind( count = table(remo.te2.filt$TEfam),
                                     CountC1 = tapply(remo.te2.filt$C1,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     CountC2 = tapply(remo.te2.filt$C2,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     CountC3 = tapply(remo.te2.filt$C3,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     CountS1 = tapply(remo.te2.filt$S1,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     CountS2 = tapply(remo.te2.filt$S2,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     CountS3 = tapply(remo.te2.filt$S3,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                     FreqC1 = tapply(remo.te2.filt$C1,remo.te2.filt$TEfam,mean),
                                     FreqC2 = tapply(remo.te2.filt$C2,remo.te2.filt$TEfam,mean), 
                                     FreqC3 = tapply(remo.te2.filt$C3,remo.te2.filt$TEfam,mean), 
                                     FreqS1 = tapply(remo.te2.filt$S1,remo.te2.filt$TEfam,mean),
                                     FreqS2 = tapply(remo.te2.filt$S2,remo.te2.filt$TEfam,mean), 
                                     FreqS3 = tapply(remo.te2.filt$S3,remo.te2.filt$TEfam,mean), 
                                     FreqC = tapply(remo.te2.filt$AvFreqC,remo.te2.filt$TEfam,mean),
                                     FreqS =tapply(remo.te2.filt$AvFreqS,remo.te2.filt$TEfam,mean),
                                     FreqAv = tapply(remo.te2.filt$AvFreq,remo.te2.filt$TEfam,mean)))

remo.te2.stat.pop$CountC = rowMeans(remo.te2.stat.pop[,c(2:4)])
remo.te2.stat.pop$CountS = rowMeans(remo.te2.stat.pop[,c(5:7)])

remo.te2.stat.pop$TEfam = row.names(remo.te2.stat.pop)
remo.te2.stat.pop$dFreq = remo.te2.stat.pop$FreqS - remo.te2.stat.pop$FreqC
head(remo.te2.stat.pop,25)
remo.te2.stat.pop = remo.te2.stat.pop[remo.te2.stat.pop$count > 0,]
# write.table(remo.te2.stat.pop,"/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/remo_average_freq_counts.txt",quote=F,row.names = F,sep="\t")

#Hoedjes
CountC = NULL
CountS = NULL
FreqCpop = NULL
FreqSpop = NULL
for (i in 1:length(hoed.cont.vect)){
  cont.name = hoed.cont.vect[i]
  sel.name = hoed.sel.vect[i]
  
  C.int = tapply(hoed.te2.filt[,cont.name],hoed.te2.filt$TEfam,function(x) length(x[x>0]))
  S.int = tapply(hoed.te2.filt[,sel.name],hoed.te2.filt$TEfam,function(x) length(x[x>0]))
  
  C.freq.int = tapply(hoed.te2.filt[,cont.name],hoed.te2.filt$TEfam,mean)
  S.freq.int = tapply(hoed.te2.filt[,sel.name],hoed.te2.filt$TEfam,mean)
  
  CountC = cbind(CountC,C.int)
  CountS = cbind(CountS,S.int)
  
  FreqCpop = cbind(FreqCpop,C.freq.int)
  FreqSpop = cbind(FreqSpop,S.freq.int)
}
colnames(CountC) = paste("Count",hoed.cont.vect,sep="")
colnames(CountS) = paste("Count",hoed.sel.vect,sep="")
colnames(FreqCpop) = paste("Freq",hoed.cont.vect,sep="")
colnames(FreqSpop) = paste("Freq",hoed.sel.vect,sep="")

hoed.te2.stat.pop = as.data.frame(cbind( count = table(hoed.te2.filt$TEfam),
                                         CountC,
                                         CountS,
                                         FreqCpop,
                                         FreqSpop,
                                         FreqC = tapply(hoed.te2.filt$AvFreqC,hoed.te2.filt$TEfam,mean),
                                         FreqS =tapply(hoed.te2.filt$AvFreqS,hoed.te2.filt$TEfam,mean),
                                         FreqAv = tapply(hoed.te2.filt$AvFreq,hoed.te2.filt$TEfam,mean)))
                                         
hoed.te2.stat.pop$CountC = rowMeans(hoed.te2.stat.pop[,paste("Count",hoed.cont.vect,sep="")])
hoed.te2.stat.pop$CountS = rowMeans(hoed.te2.stat.pop[,paste("Count",hoed.sel.vect,sep="")])

hoed.te2.stat.pop$TEfam = row.names(hoed.te2.stat.pop)
hoed.te2.stat.pop$dFreq = hoed.te2.stat.pop$FreqS - hoed.te2.stat.pop$FreqC
head(hoed.te2.stat.pop,25)
hoed.te2.stat.pop = hoed.te2.stat.pop[hoed.te2.stat.pop$count > 0,]
#write.table(hoed.te2.stat.pop,"/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/PopTE2/hoed_average_freq_counts.txt",quote=F,row.names = F,sep="\t")

########################################################
#USE SUM OF FREQUENCIES FOR EACHT TE FAMILY TO GET ABUNDANCE (REVIEWER SUGGESTION)
########################################################
#Carnes - SUM OF FREQUENCIES
carnes.te2.stat.pop.sum = as.data.frame(cbind( count = table(carnes.te2.filt$TEfam),
                                               CountC1 = tapply(carnes.te2.filt$C1,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountC2 = tapply(carnes.te2.filt$C2,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountC3 = tapply(carnes.te2.filt$C3,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountC4 = tapply(carnes.te2.filt$C4,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountC5 = tapply(carnes.te2.filt$C5,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountS1 = tapply(carnes.te2.filt$S1,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountS2 = tapply(carnes.te2.filt$S2,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountS3 = tapply(carnes.te2.filt$S3,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountS4 = tapply(carnes.te2.filt$S4,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountS5 = tapply(carnes.te2.filt$S5,carnes.te2.filt$TEfam,function(x) length(x[x>0])),
                                               SumC1 = tapply(carnes.te2.filt$C1,carnes.te2.filt$TEfam,sum),
                                               SumC2 = tapply(carnes.te2.filt$C2,carnes.te2.filt$TEfam,sum),
                                               SumC3 = tapply(carnes.te2.filt$C3,carnes.te2.filt$TEfam,sum),
                                               SumC4 = tapply(carnes.te2.filt$C4,carnes.te2.filt$TEfam,sum),
                                               SumC5 = tapply(carnes.te2.filt$C5,carnes.te2.filt$TEfam,sum),
                                               SumS1 = tapply(carnes.te2.filt$S1,carnes.te2.filt$TEfam,sum),
                                               SumS2 = tapply(carnes.te2.filt$S2,carnes.te2.filt$TEfam,sum),
                                               SumS3 = tapply(carnes.te2.filt$S3,carnes.te2.filt$TEfam,sum),
                                               SumS4 = tapply(carnes.te2.filt$S4,carnes.te2.filt$TEfam,sum),
                                               SumS5 = tapply(carnes.te2.filt$S5,carnes.te2.filt$TEfam,sum)))

carnes.te2.stat.pop.sum$SumC = rowMeans(carnes.te2.stat.pop.sum[,c(12:16)])
carnes.te2.stat.pop.sum$SumS = rowMeans(carnes.te2.stat.pop.sum[,c(17:21)])
carnes.te2.stat.pop.sum$CountC = rowMeans(carnes.te2.stat.pop.sum[,c(2:6)])
carnes.te2.stat.pop.sum$CountS = rowMeans(carnes.te2.stat.pop.sum[,c(7:11)])
carnes.te2.stat.pop.sum$TEfam = row.names(carnes.te2.stat.pop.sum)

carnes.te2.stat.pop.sum$dSum = carnes.te2.stat.pop.sum$SumS - carnes.te2.stat.pop.sum$SumC
carnes.te2.stat.pop.sum = carnes.te2.stat.pop.sum[!is.na(carnes.te2.stat.pop.sum$dSum),]
head(carnes.te2.stat.pop.sum)
table(carnes.te2.stat.pop.sum$dSum>0) #88 more sum in selected, 4 less sum

table(rowMin(carnes.te2.stat.pop.sum,colstart = 12,colend = 16) > rowMax(carnes.te2.stat.pop.sum,colstart = 17,colend = 21))
#0 consistently higher in C
table(rowMax(carnes.te2.stat.pop.sum,colstart = 12,colend = 16) < rowMin(carnes.te2.stat.pop.sum,colstart = 17,colend = 21))
#7 consistently higher frequency in S

#stats
pval.carnes = NULL
for (i in 1:nrow(carnes.te2.stat.pop.sum)){
  pval = t.test(carnes.te2.stat.pop.sum[i,c(12:16)], carnes.te2.stat.pop.sum[i,c(17:21)])$p.value
  pval.carnes = c(pval.carnes,pval)
}
carnes.te2.stat.pop.sum$Pval = pval.carnes
nrow(carnes.te2.stat.pop.sum[carnes.te2.stat.pop.sum$Pval < 0.05,]) #24 more in selected, 0 less

#Check TEs also detected by DeviaTE
table(carnes.te2.stat.pop.sum[row.names(carnes.te2.stat.pop.sum) %in% carnes.stat.filt$TEfam,]$dSum>0) #4 smaller and 77 larger
table(rowMin(carnes.te2.stat.pop.sum[row.names(carnes.te2.stat.pop.sum) %in% carnes.stat.filt$TEfam,],colstart = 12,colend = 16) > rowMax(carnes.te2.stat.pop.sum[row.names(carnes.te2.stat.pop.sum) %in% carnes.stat.filt$TEfam,],colstart = 17,colend = 21))
#0 consistently higher in C
table(rowMax(carnes.te2.stat.pop.sum[row.names(carnes.te2.stat.pop.sum) %in% carnes.stat.filt$TEfam,],colstart = 12,colend = 16) < rowMin(carnes.te2.stat.pop.sum[row.names(carnes.te2.stat.pop.sum) %in% carnes.stat.filt$TEfam,],colstart = 17,colend = 21))
#6 consistently higher frequency in S

#Fabian - SUM OF FREQUENCIES
fabian.te2.stat.pop.sum = as.data.frame(cbind( count = table(fabian.te2.filt$TEfam),
                                               CountRa = tapply(fabian.te2.filt$Ra,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountRb = tapply(fabian.te2.filt$Rb,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountLa = tapply(fabian.te2.filt$La,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               CountLb = tapply(fabian.te2.filt$Lb,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               Count2La = tapply(fabian.te2.filt$`2La`,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               Count2Lb = tapply(fabian.te2.filt$`2Lb`,fabian.te2.filt$TEfam,function(x) length(x[x>0])),
                                               SumRa = tapply(fabian.te2.filt$Ra,fabian.te2.filt$TEfam,sum),
                                               SumRb = tapply(fabian.te2.filt$Rb,fabian.te2.filt$TEfam,sum),
                                               SumLa = tapply(fabian.te2.filt$La,fabian.te2.filt$TEfam,sum),
                                               SumLb = tapply(fabian.te2.filt$Lb,fabian.te2.filt$TEfam,sum),
                                               Sum2La = tapply(fabian.te2.filt$`2La`,fabian.te2.filt$TEfam,sum),
                                               Sum2Lb = tapply(fabian.te2.filt$`2Lb`,fabian.te2.filt$TEfam,sum)))

fabian.te2.stat.pop.sum$SumC = rowMeans(fabian.te2.stat.pop.sum[,c(8:9)])
fabian.te2.stat.pop.sum$SumS = rowMeans(fabian.te2.stat.pop.sum[,c(10:13)])
fabian.te2.stat.pop.sum$CountC = rowMeans(fabian.te2.stat.pop.sum[,c(2:3)])
fabian.te2.stat.pop.sum$CountS = rowMeans(fabian.te2.stat.pop.sum[,c(4:7)])
fabian.te2.stat.pop.sum$TEfam = row.names(fabian.te2.stat.pop.sum)

fabian.te2.stat.pop.sum$dSum = fabian.te2.stat.pop.sum$SumS - fabian.te2.stat.pop.sum$SumC
fabian.te2.stat.pop.sum = fabian.te2.stat.pop.sum[!is.na(fabian.te2.stat.pop.sum$dSum),]
head(fabian.te2.stat.pop.sum)
table(fabian.te2.stat.pop.sum$dSum>0) #80 more sum in selected, 43 less sum

table(rowMin(fabian.te2.stat.pop.sum,colstart = 8,colend = 9) > rowMax(fabian.te2.stat.pop.sum,colstart = 10,colend = 13))
#10 consistently higher in C
table(rowMax(fabian.te2.stat.pop.sum,colstart = 8,colend = 9) < rowMin(fabian.te2.stat.pop.sum,colstart = 10,colend = 13))
#20 consistently higher frequency in S

#stats
pval.fabian = NULL
for (i in 1:nrow(fabian.te2.stat.pop.sum)){
  if (row.names(fabian.te2.stat.pop.sum[i,]) == "STALKER3") {
    pval = NA
  } else {
    pval = t.test(fabian.te2.stat.pop.sum[i,c(8:9)], fabian.te2.stat.pop.sum[i,c(10:13)])$p.value
  }
  pval.fabian = c(pval.fabian,pval)
}
fabian.te2.stat.pop.sum$Pval = pval.fabian
fabian.te2.stat.pop.sum[fabian.te2.stat.pop.sum$Pval < 0.05,]
table(fabian.te2.stat.pop.sum[fabian.te2.stat.pop.sum$Pval < 0.05,]$dSum>0) #9 more in selected, 2 less

#Check TEs also detected by DeviaTE
table(fabian.te2.stat.pop.sum[row.names(fabian.te2.stat.pop.sum) %in% fabian.stat.filt$TEfam,]$dSum>0) #37 smaller and 72 larger
table(rowMin(fabian.te2.stat.pop.sum[row.names(fabian.te2.stat.pop.sum) %in% fabian.stat.filt$TEfam,],colstart = 8,colend = 9) > rowMax(fabian.te2.stat.pop.sum[row.names(fabian.te2.stat.pop.sum) %in% fabian.stat.filt$TEfam,],colstart = 10,colend = 13))
#8 consistently higher in C
table(rowMax(fabian.te2.stat.pop.sum[row.names(fabian.te2.stat.pop.sum) %in% fabian.stat.filt$TEfam,],colstart = 8,colend = 9) < rowMin(fabian.te2.stat.pop.sum[row.names(fabian.te2.stat.pop.sum) %in% fabian.stat.filt$TEfam,],colstart = 10,colend = 13))
#19 consistently higher frequency in S

#Remolina - SUM OF FREQUENCIES
remo.te2.stat.pop.sum = as.data.frame(cbind( count = table(remo.te2.filt$TEfam),
                                             CountC1 = tapply(remo.te2.filt$C1,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             CountC2 = tapply(remo.te2.filt$C2,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             CountC3 = tapply(remo.te2.filt$C3,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             CountS1 = tapply(remo.te2.filt$S1,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             CountS2 = tapply(remo.te2.filt$S2,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             CountS3 = tapply(remo.te2.filt$S3,remo.te2.filt$TEfam,function(x) length(x[x>0])),
                                             SumC1 = tapply(remo.te2.filt$C1,remo.te2.filt$TEfam,sum),
                                             SumC2 = tapply(remo.te2.filt$C2,remo.te2.filt$TEfam,sum), 
                                             SumC3 = tapply(remo.te2.filt$C3,remo.te2.filt$TEfam,sum), 
                                             SumS1 = tapply(remo.te2.filt$S1,remo.te2.filt$TEfam,sum),
                                             SumS2 = tapply(remo.te2.filt$S2,remo.te2.filt$TEfam,sum), 
                                             SumS3 = tapply(remo.te2.filt$S3,remo.te2.filt$TEfam,sum)))

remo.te2.stat.pop.sum$SumC = rowMeans(remo.te2.stat.pop.sum[,c(8:10)])
remo.te2.stat.pop.sum$SumS = rowMeans(remo.te2.stat.pop.sum[,c(11:13)])
remo.te2.stat.pop.sum$CountC = rowMeans(remo.te2.stat.pop.sum[,c(2:4)])
remo.te2.stat.pop.sum$CountS = rowMeans(remo.te2.stat.pop.sum[,c(5:7)])
remo.te2.stat.pop.sum$TEfam = row.names(remo.te2.stat.pop.sum)

remo.te2.stat.pop.sum$dSum = remo.te2.stat.pop.sum$SumS - remo.te2.stat.pop.sum$SumC
remo.te2.stat.pop.sum = remo.te2.stat.pop.sum[!is.na(remo.te2.stat.pop.sum$dSum),]
head(remo.te2.stat.pop.sum,25)
table(remo.te2.stat.pop.sum$dSum>0) #60 more, 60 less

table(rowMin(remo.te2.stat.pop.sum,colstart = 8,colend = 10) > rowMax(remo.te2.stat.pop.sum,colstart = 11,colend = 13))
#5 consistently higher in C
table(rowMax(remo.te2.stat.pop.sum,colstart = 8,colend = 10) < rowMin(remo.te2.stat.pop.sum,colstart = 11,colend = 13))
#4 consistently higher frequency in S

#stats
pval.remo = NULL
for (i in 1:nrow(remo.te2.stat.pop.sum)){
  pval = t.test(remo.te2.stat.pop.sum[i,c(8:10)], remo.te2.stat.pop.sum[i,c(11:13)])$p.value
  pval.remo = c(pval.remo,pval)
}
remo.te2.stat.pop.sum$Pval = pval.remo
remo.te2.stat.pop.sum[remo.te2.stat.pop.sum$Pval < 0.05,] #3 more in selected, 1 less (G2)

#Check TEs also detected by DeviaTE
table(remo.te2.stat.pop.sum[row.names(remo.te2.stat.pop.sum) %in% remo.stat.filt$TEfam,]$dSum>0) #53 smaller and 55 larger
table(rowMin(remo.te2.stat.pop.sum[row.names(remo.te2.stat.pop.sum) %in% remo.stat.filt$TEfam,],colstart = 8,colend = 10) > rowMax(remo.te2.stat.pop.sum[row.names(remo.te2.stat.pop.sum) %in% remo.stat.filt$TEfam,],colstart = 11,colend = 13))
#5 consistently higher in C
table(rowMax(remo.te2.stat.pop.sum[row.names(remo.te2.stat.pop.sum) %in% remo.stat.filt$TEfam,],colstart = 8,colend = 10) < rowMin(remo.te2.stat.pop.sum[row.names(remo.te2.stat.pop.sum) %in% remo.stat.filt$TEfam,],colstart = 11,colend = 13))
#4 consistently higher frequency in S

#Hoedjes - SUM OF FREQUENCIES
CountC = NULL
CountS = NULL
SumCpop = NULL
SumSpop = NULL
for (i in 1:length(hoed.cont.vect)){
  cont.name = hoed.cont.vect[i]
  sel.name = hoed.sel.vect[i]
  
  C.int = tapply(hoed.te2.filt[,cont.name],hoed.te2.filt$TEfam,function(x) length(x[x>0]))
  S.int = tapply(hoed.te2.filt[,sel.name],hoed.te2.filt$TEfam,function(x) length(x[x>0]))
  
  C.sum.int = tapply(hoed.te2.filt[,cont.name],hoed.te2.filt$TEfam,sum)
  S.sum.int = tapply(hoed.te2.filt[,sel.name],hoed.te2.filt$TEfam,sum)
  
  CountC = cbind(CountC,C.int)
  CountS = cbind(CountS,S.int)
  
  SumCpop = cbind(SumCpop,C.sum.int)
  SumSpop = cbind(SumSpop,S.sum.int)
}
colnames(CountC) = paste("Count",hoed.cont.vect,sep="")
colnames(CountS) = paste("Count",hoed.sel.vect,sep="")
colnames(SumCpop) = paste("Sum",hoed.cont.vect,sep="")
colnames(SumSpop) = paste("Sum",hoed.sel.vect,sep="")

hoed.te2.stat.pop.sum = as.data.frame(cbind( count = table(hoed.te2.filt$TEfam),
                                             CountC,
                                             CountS,
                                             SumCpop,
                                             SumSpop))

hoed.te2.stat.pop.sum$SumC = rowMeans(hoed.te2.stat.pop.sum[,c('SumCE1' , 'SumCE2' , 'SumCE3' , 'SumCE4' , 'SumLE1' , 'SumLE2' , 'SumLE3' , 'SumLE4' , 'SumHE1' , 'SumHE2' , 'SumHE3' , 'SumHE4')])
hoed.te2.stat.pop.sum$SumS = rowMeans(hoed.te2.stat.pop.sum[,c('SumCP1' , 'SumCP2' , 'SumCP3' , 'SumCP4' , 'SumLP1' , 'SumLP2' , 'SumLP3' , 'SumLP4' , 'SumHP1' , 'SumHP2' , 'SumHP3' , 'SumHP4')])
hoed.te2.stat.pop.sum$CountC = rowMeans(hoed.te2.stat.pop.sum[,paste("Count",hoed.cont.vect,sep="")])
hoed.te2.stat.pop.sum$CountS = rowMeans(hoed.te2.stat.pop.sum[,paste("Count",hoed.sel.vect,sep="")])
hoed.te2.stat.pop.sum$TEfam = row.names(hoed.te2.stat.pop.sum)

hoed.te2.stat.pop.sum$dSum = hoed.te2.stat.pop.sum$SumS - hoed.te2.stat.pop.sum$SumC
hoed.te2.stat.pop.sum = hoed.te2.stat.pop.sum[hoed.te2.stat.pop.sum$count > 0,]
head(hoed.te2.stat.pop.sum,25)
table(hoed.te2.stat.pop.sum$dSum>0) #99 more, 22 less

#Consistency across all C and all S
table(rowMin(hoed.te2.stat.pop.sum,colstart = 26,colend = 37) > rowMax(hoed.te2.stat.pop.sum,colstart = 38,colend = 49))
#0 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum,colstart = 26,colend = 37) < rowMin(hoed.te2.stat.pop.sum,colstart = 38,colend = 49))
#0 consistently higher frequency in S

#Consistency in control diet 
table(rowMin(hoed.te2.stat.pop.sum,colstart = 26,colend = 29) > rowMax(hoed.te2.stat.pop.sum,colstart = 38,colend = 41))
#0 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum,colstart = 26,colend = 29) < rowMin(hoed.te2.stat.pop.sum,colstart = 38,colend = 41))
#22 consistently higher frequency in S

#Consistency in low diet 
table(rowMin(hoed.te2.stat.pop.sum,colstart = 30,colend = 33) > rowMax(hoed.te2.stat.pop.sum,colstart = 42,colend = 45))
#2 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum,colstart = 30,colend = 33) < rowMin(hoed.te2.stat.pop.sum,colstart = 42,colend = 45))
#22 consistently higher frequency in S

#Consistency in high diet 
table(rowMin(hoed.te2.stat.pop.sum,colstart = 34,colend = 37) > rowMax(hoed.te2.stat.pop.sum,colstart = 46,colend = 49))
#1 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum,colstart = 34,colend = 37) < rowMin(hoed.te2.stat.pop.sum,colstart = 46,colend = 49))
#8 consistently higher frequency in S

#stats
pval.hoed = NULL
for (i in 1:nrow(hoed.te2.stat.pop.sum)){
  int = as.data.frame(t(hoed.te2.stat.pop.sum[i,c(26:49)]))
  diet = mgsub(pattern = c("SumC[A-Z][1-9]","SumL[A-Z][1-9]","SumH[A-Z][1-9]"),replacement = c("C","L","H"),row.names(int))
  regime = mgsub(pattern = c("Sum[A-Z]E[1-9]","Sum[A-Z]P[1-9]"),replacement = c("Cont","Sel"),row.names(int))
  int$diet = diet
  int$regime= regime
  model = aov(int[,1] ~ diet + regime + diet*regime,data=int)
  pval = summary(model)[[1]][['Pr(>F)']][2]
  pval.hoed = c(pval.hoed,pval)
}
hoed.te2.stat.pop.sum$Pval = pval.hoed
nrow(hoed.te2.stat.pop.sum[hoed.te2.stat.pop.sum$Pval < 0.05,]) #37 more in selected, 1 less

#Check TEs also detected by DeviaTE
table(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,]$dSum>0) #19 smaller and 94 larger
#Consistency across all C and all S
table(rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 26,colend = 37) > rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 38,colend = 49))
#0 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 26,colend = 37) < rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 38,colend = 49))
#0 consistently higher frequency in S

#Consistency in control diet 
table(rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 26,colend = 29) > rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 38,colend = 41))
#0 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 26,colend = 29) < rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 38,colend = 41))
#22 consistently higher frequency in S

#Consistency in low diet 
table(rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 30,colend = 33) > rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 42,colend = 45))
#2 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 30,colend = 33) < rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 42,colend = 45))
#22 consistently higher frequency in S

#Consistency in high diet 
table(rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 34,colend = 37) > rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 46,colend = 49))
#1 consistently higher in C
table(rowMax(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 34,colend = 37) < rowMin(hoed.te2.stat.pop.sum[row.names(hoed.te2.stat.pop.sum) %in% hoed.stat.filt$TEfam,],colstart = 46,colend = 49))
#8 consistently higher frequency in S

#Merge PoPoolationTE2 and DeviaTE tab
carnes.te2.stat.pop.combi = merge(carnes.te2.stat.pop, carnes.stat.filt, by ="TEfam",all=T)
fabian.te2.stat.pop.combi = merge(fabian.te2.stat.pop, fabian.stat.filt, by ="TEfam",all=T)
remo.te2.stat.pop.combi = merge(remo.te2.stat.pop, remo.stat.filt, by ="TEfam",all=T)
hoed.te2.stat.pop.combi = merge(hoed.te2.stat.pop, hoed.stat.filt, by ="TEfam",all=T) 

#Subset to TE families that were detected in PoPoolationTE2 and DeviaTE
carnes.te2.stat.pop.combi = carnes.te2.stat.pop.combi[!is.na(carnes.te2.stat.pop.combi$Mean_Cont),] 
fabian.te2.stat.pop.combi = fabian.te2.stat.pop.combi[!is.na(fabian.te2.stat.pop.combi$Mean_Cont),] 
remo.te2.stat.pop.combi = remo.te2.stat.pop.combi[!is.na(remo.te2.stat.pop.combi$Mean_Cont),] 
hoed.te2.stat.pop.combi = hoed.te2.stat.pop.combi[!is.na(hoed.te2.stat.pop.combi$Mean_Early),] 

nrow(carnes.te2.stat.pop.combi) #112
nrow(fabian.te2.stat.pop.combi) #110
nrow(remo.te2.stat.pop.combi) #110
nrow(hoed.te2.stat.pop.combi) #115


#####################################################
#Correlation between DeviaTE abundance and frequency from PoPTE2
#####################################################
cor.test(carnes.te2.stat.pop.combi$Mean_Cont, carnes.te2.stat.pop.combi$FreqC, method = "spearman") #0.1293586 , ns 
cor.test(carnes.te2.stat.pop.combi$Mean_Sel, carnes.te2.stat.pop.combi$FreqS, method = "spearman") #-0.2362014 , p 0.033
cor.test(carnes.te2.stat.pop.combi$Cont_B1, carnes.te2.stat.pop.combi$FreqC1, method = "spearman") #0.09435464, ns
cor.test(carnes.te2.stat.pop.combi$Cont_B2, carnes.te2.stat.pop.combi$FreqC2, method = "spearman") #-0.1569575, ns 
cor.test(carnes.te2.stat.pop.combi$Cont_B3, carnes.te2.stat.pop.combi$FreqC3, method = "spearman") #-0.01225207, ns
cor.test(carnes.te2.stat.pop.combi$Cont_B4, carnes.te2.stat.pop.combi$FreqC4, method = "spearman") #0.122473, ns
cor.test(carnes.te2.stat.pop.combi$Cont_B5, carnes.te2.stat.pop.combi$FreqC5, method = "spearman") #0.2812264 , p 0.011
cor.test(carnes.te2.stat.pop.combi$Sel_O1, carnes.te2.stat.pop.combi$FreqS1, method = "spearman") #-0.1418343 , ns
cor.test(carnes.te2.stat.pop.combi$Sel_O2, carnes.te2.stat.pop.combi$FreqS2, method = "spearman") #-0.1592117 , ns 
cor.test(carnes.te2.stat.pop.combi$Sel_O3, carnes.te2.stat.pop.combi$FreqS3, method = "spearman") #-0.09327115 , ns
cor.test(carnes.te2.stat.pop.combi$Sel_O4, carnes.te2.stat.pop.combi$FreqS4, method = "spearman") #-0.1365595, ns
cor.test(carnes.te2.stat.pop.combi$Sel_O5, carnes.te2.stat.pop.combi$FreqS5, method = "spearman") #-0.2853156 p = 0.009828
mean(c(0.09435464,-0.1569575,-0.01225207,0.122473,0.2812264,
  -0.1418343,-0.1592117,-0.09327115,-0.1365595,-0.2853156))
minmax(c(0.09435464,-0.1569575,-0.01225207,0.122473,0.2812264,
       -0.1418343,-0.1592117,-0.09327115,-0.1365595,-0.2853156))
#mean = -0.05; minmax = -0.29  0.28

cor.test(fabian.te2.stat.pop.combi$Mean_Cont, fabian.te2.stat.pop.combi$FreqC, method = "spearman") #-0.3697618 
cor.test(fabian.te2.stat.pop.combi$Mean_Sel, fabian.te2.stat.pop.combi$FreqS, method = "spearman") #-0.3355397
cor.test(fabian.te2.stat.pop.combi$Cont_Ra, fabian.te2.stat.pop.combi$FreqRa,method = "spearman") #-0.383 sign
cor.test(fabian.te2.stat.pop.combi$Cont_Rb, fabian.te2.stat.pop.combi$FreqRb,method = "spearman") #-0.343 sign
cor.test(fabian.te2.stat.pop.combi$Sel_La, fabian.te2.stat.pop.combi$FreqLa,method = "spearman") #-0.353489 
cor.test(fabian.te2.stat.pop.combi$Sel_Lb, fabian.te2.stat.pop.combi$FreqLb,method = "spearman") #-0.3309425 
cor.test(fabian.te2.stat.pop.combi$Sel_2La, fabian.te2.stat.pop.combi$Freq2La,method = "spearman") #-0.2772786  
cor.test(fabian.te2.stat.pop.combi$Sel_2Lb, fabian.te2.stat.pop.combi$Freq2Lb,method = "spearman") #-0.31752 
mean(-c(0.383,0.343,0.353,0.331,0.278,0.318))
minmax(-c(0.383,0.343,0.353,0.331,0.278,0.318))

cor.test(remo.te2.stat.pop.combi$Mean_Cont, remo.te2.stat.pop.combi$FreqC,method = "spearman") #-0.4041842 
cor.test(remo.te2.stat.pop.combi$Mean_Sel, remo.te2.stat.pop.combi$FreqS,method = "spearman") #-0.4017453 
cor.test(remo.te2.stat.pop.combi$Cont_1, remo.te2.stat.pop.combi$FreqC1,method = "spearman") #-0.4142731 sign - 
cor.test(remo.te2.stat.pop.combi$Cont_2, remo.te2.stat.pop.combi$FreqC2,method = "spearman") #-0.3945049 sign - 
cor.test(remo.te2.stat.pop.combi$Cont_3, remo.te2.stat.pop.combi$FreqC3,method = "spearman") #-0.3949155 sign - 
cor.test(remo.te2.stat.pop.combi$Sel_1, remo.te2.stat.pop.combi$FreqS1,method = "spearman") #-0.3966199 sign - 
cor.test(remo.te2.stat.pop.combi$Sel_2, remo.te2.stat.pop.combi$FreqS2,method = "spearman") #-0.3906752 sign - 
cor.test(remo.te2.stat.pop.combi$Sel_3, remo.te2.stat.pop.combi$FreqS3,method = "spearman") #sign - 
mean(c(-0.4142731, -0.3945049, -0.3949155, -0.3966199, -0.3966199, -0.3906752, -0.414797))
minmax(c(-0.4142731, -0.3945049, -0.3949155, -0.3966199, -0.3966199, -0.3906752, -0.414797))

cor.test(hoed.te2.stat.pop.combi$Mean_Early, hoed.te2.stat.pop.combi$FreqC,method = "spearman") #-0.3746839 
cor.test(hoed.te2.stat.pop.combi$Mean_Posponed, hoed.te2.stat.pop.combi$FreqS,method = "spearman") #-0.3909026
mean(as.numeric(c(cor.test(hoed.te2.stat.pop.combi$CE1, hoed.te2.stat.pop.combi$FreqCE1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CE2, hoed.te2.stat.pop.combi$FreqCE2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CE3, hoed.te2.stat.pop.combi$FreqCE3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CE4, hoed.te2.stat.pop.combi$FreqCE4,method = "spearman")$estimate, #sign -, rho = -0.2427728, p= 0.009573 
cor.test(hoed.te2.stat.pop.combi$CP1, hoed.te2.stat.pop.combi$FreqCP1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CP2, hoed.te2.stat.pop.combi$FreqCP2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CP3, hoed.te2.stat.pop.combi$FreqCP3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$CP4, hoed.te2.stat.pop.combi$FreqCP4,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LE1, hoed.te2.stat.pop.combi$FreqLE1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LE2, hoed.te2.stat.pop.combi$FreqLE2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LE3, hoed.te2.stat.pop.combi$FreqLE3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LE4, hoed.te2.stat.pop.combi$FreqLE4,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LP1, hoed.te2.stat.pop.combi$FreqLP1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LP2, hoed.te2.stat.pop.combi$FreqLP2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LP3, hoed.te2.stat.pop.combi$FreqLP3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$LP4, hoed.te2.stat.pop.combi$FreqLP4,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HE1, hoed.te2.stat.pop.combi$FreqHE1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HE2, hoed.te2.stat.pop.combi$FreqHE2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HE3, hoed.te2.stat.pop.combi$FreqHE3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HE4, hoed.te2.stat.pop.combi$FreqHE4,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HP1, hoed.te2.stat.pop.combi$FreqHP1,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HP2, hoed.te2.stat.pop.combi$FreqHP2,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HP3, hoed.te2.stat.pop.combi$FreqHP3,method = "spearman")$estimate, #sign - 
cor.test(hoed.te2.stat.pop.combi$HP4, hoed.te2.stat.pop.combi$FreqHP4,method = "spearman")$estimate))) #sign - 
#mean is -0.3723869 - all significant
#minmax -0.4082441 -0.2427728


#####################################################
#Preparation for check if difference in frequency between C>S and S>C
#####################################################
carnes.te2.stat.pop.combi.ns = carnes.te2.stat.pop.combi[carnes.te2.stat.pop.combi$Bonf == "FALSE",]
carnes.te2.stat.pop.combi.ns$Type = "ns"
fabian.te2.stat.pop.combi.ns = fabian.te2.stat.pop.combi[fabian.te2.stat.pop.combi$Bonf == "FALSE",]
fabian.te2.stat.pop.combi.ns$Type = "ns"
remo.te2.stat.pop.combi.ns = remo.te2.stat.pop.combi[remo.te2.stat.pop.combi$Bonf == "FALSE",]
remo.te2.stat.pop.combi.ns$Type = "ns"
hoed.te2.stat.pop.combi.ns = hoed.te2.stat.pop.combi[hoed.te2.stat.pop.combi$Bonf_Full == "FALSE",]
hoed.te2.stat.pop.combi.ns$Type = "ns"
 
carnes.te2.stat.pop.combi.sign = carnes.te2.stat.pop.combi[carnes.te2.stat.pop.combi$Bonf == "TRUE",]
fabian.te2.stat.pop.combi.sign = fabian.te2.stat.pop.combi[fabian.te2.stat.pop.combi$Bonf == "TRUE",]
remo.te2.stat.pop.combi.sign = remo.te2.stat.pop.combi[remo.te2.stat.pop.combi$Bonf == "TRUE",]
hoed.te2.stat.pop.combi.sign = hoed.te2.stat.pop.combi[hoed.te2.stat.pop.combi$Bonf_Full == "TRUE",]

carnes.te2.stat.pop.combi.sign$Type = mgsub(c("TRUE","FALSE"),c("S>C","C>S"),carnes.te2.stat.pop.combi.sign$Diff_SelCont >0)
fabian.te2.stat.pop.combi.sign$Type = mgsub(c("TRUE","FALSE"),c("S>C","C>S"),fabian.te2.stat.pop.combi.sign$Diff_SelCont >0)
remo.te2.stat.pop.combi.sign$Type = mgsub(c("TRUE","FALSE"),c("S>C","C>S"),remo.te2.stat.pop.combi.sign$Diff_SelCont >0)
hoed.te2.stat.pop.combi.sign$Type = mgsub(c("TRUE","FALSE"),c("S>C","C>S"),hoed.te2.stat.pop.combi.sign$Diff_PostEarly >0)

carnes.te2.stat.pop.combi = rbind(carnes.te2.stat.pop.combi.sign, carnes.te2.stat.pop.combi.ns)
fabian.te2.stat.pop.combi = rbind(fabian.te2.stat.pop.combi.sign, fabian.te2.stat.pop.combi.ns)
remo.te2.stat.pop.combi = rbind(remo.te2.stat.pop.combi.sign, remo.te2.stat.pop.combi.ns)
hoed.te2.stat.pop.combi = rbind(hoed.te2.stat.pop.combi.sign, hoed.te2.stat.pop.combi.ns)

################################################################
#Correlation with average frequencies from Kofler et al 2015 (PLoS Genetics)
################################################################

#Carnes
carnes.te2.stat.pop.combi.rob = merge(carnes.te2.stat.pop.combi,robert.2015,by="TEfam")
carnes.te2.stat.pop.combi.sign.rob = merge(carnes.te2.stat.pop.combi.sign,robert.2015,by="TEfam")
cor.test(carnes.te2.stat.pop.combi.rob$FreqAv, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #r = 0.09535112, p = 0.403

cor.test(carnes.te2.stat.pop.combi.rob$FreqC1, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #neg, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqC2, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqC3, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqC4, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqC5, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #neg, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqS1, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqS2, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqS3, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqS4, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #pos, ns
cor.test(carnes.te2.stat.pop.combi.rob$FreqS5, carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") #r = 0.2662396, significant, 0.01771

#Fabian
fabian.te2.stat.pop.combi.rob = merge(fabian.te2.stat.pop.combi,robert.2015,by="TEfam")
fabian.te2.stat.pop.combi.sign.rob = merge(fabian.te2.stat.pop.combi.sign,robert.2015,by="TEfam") 
cor.test(fabian.te2.stat.pop.combi.rob$FreqAv, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.65 and p = 0.000000000000083623

cor.test(fabian.te2.stat.pop.combi.rob$FreqC, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.65, sign
cor.test(fabian.te2.stat.pop.combi.rob$FreqS, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.65, sign
cor.test(fabian.te2.stat.pop.combi.rob$FreqRa, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.64
cor.test(fabian.te2.stat.pop.combi.rob$FreqRb, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.66
cor.test(fabian.te2.stat.pop.combi.rob$FreqLa, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.65
cor.test(fabian.te2.stat.pop.combi.rob$FreqLb, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.66
cor.test(fabian.te2.stat.pop.combi.rob$Freq2La, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.65
cor.test(fabian.te2.stat.pop.combi.rob$Freq2Lb, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.66

#Remolina
remo.te2.stat.pop.combi.rob = merge(remo.te2.stat.pop.combi,robert.2015,by="TEfam")
remo.te2.stat.pop.combi.sign.rob = merge(remo.te2.stat.pop.combi.sign,robert.2015,by="TEfam") #sign
cor.test(remo.te2.stat.pop.combi.rob$FreqAv, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.5825918  and p = 1.076e-10

cor.test(remo.te2.stat.pop.combi.rob$FreqC, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.5825369   and p = 1.082e-10
cor.test(remo.te2.stat.pop.combi.rob$FreqC1, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.58
cor.test(remo.te2.stat.pop.combi.rob$FreqC2, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.57
cor.test(remo.te2.stat.pop.combi.rob$FreqC3, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.61
cor.test(remo.te2.stat.pop.combi.rob$FreqS1, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.59
cor.test(remo.te2.stat.pop.combi.rob$FreqS2, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.57
cor.test(remo.te2.stat.pop.combi.rob$FreqS3, remo.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.57

#Hoedjes
hoed.te2.stat.pop.combi.rob = merge(hoed.te2.stat.pop.combi,robert.2015,by="TEfam")
hoed.te2.stat.pop.combi.sign.rob = merge(hoed.te2.stat.pop.combi.sign,robert.2015,by="TEfam") #sign
cor.test(hoed.te2.stat.pop.combi.rob$FreqAv, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.612694  and p = 1.827e-12

cor.test(hoed.te2.stat.pop.combi.rob$FreqC, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.6093739    and p = 2.588e-12
cor.test(hoed.te2.stat.pop.combi.rob$FreqS, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.617 , p =  8.448e-16
cor.test(hoed.te2.stat.pop.combi.rob$FreqCE1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.62
cor.test(hoed.te2.stat.pop.combi.rob$FreqCE2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqCE3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.61
cor.test(hoed.te2.stat.pop.combi.rob$FreqCE4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.59
cor.test(hoed.te2.stat.pop.combi.rob$FreqCP1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqCP2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqCP3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqCP4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.63
cor.test(hoed.te2.stat.pop.combi.rob$FreqLE1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLE2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLE3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLE4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLP1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLP2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLP3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqLP4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHE1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHE2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #0.63
cor.test(hoed.te2.stat.pop.combi.rob$FreqHE3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHE4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHP1, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHP2, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHP3, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #
cor.test(hoed.te2.stat.pop.combi.rob$FreqHP4, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman") #

###
#Plot of correlation between frequencies in experimental evolution studies and Kofler2015
pdf(paste(output_path,"FigS7_correlation_with_Kofler2015_new.pdf",sep=""), width=12, height=20)
par(mfrow=c(4,2),font.axis=1,font.lab=2,lwd=1.2,cex=1,mar=c(4.5,4.5,2,3))

plot(carnes.te2.stat.pop.combi.rob$FreqAv, carnes.te2.stat.pop.combi.rob$av_freq, xlim=c(0,1),ylim = c(0,1),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Frequency (Carnes)", ylab="Frequency (SA)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
abline(0,1,lty = 2)

plot(rank(carnes.te2.stat.pop.combi.rob$FreqAv), rank(carnes.te2.stat.pop.combi.rob$av_freq), xlim= c(0,118), ylim = c(0,118),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Frequency Ranks (Carnes)", ylab="Frequency Ranks (SA)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
abline(0,1,lty = 2)
text(rank(carnes.te2.stat.pop.combi.rob$FreqAv), rank(carnes.te2.stat.pop.combi.rob$av_freq)+2, labels = carnes.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(14, 118, 
     labels = expression(paste(italic(rho), ' = 0.1 (ns)')),cex = 1.4)

plot(fabian.te2.stat.pop.combi.rob$FreqAv, fabian.te2.stat.pop.combi.rob$av_freq, 
     xlab="Frequency (Fabian)", ylab="Frequency (SA)", xlim= c(0,1), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
abline(0,1,lty = 2)

plot(rank(fabian.te2.stat.pop.combi.rob$FreqAv), rank(fabian.te2.stat.pop.combi.rob$av_freq), xlim= c(0,118), ylim = c(0,118),
     xlab="Frequency Ranks (Fabian)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
abline(0,1,lty = 2)
text(rank(fabian.te2.stat.pop.combi.rob$FreqAv), rank(fabian.te2.stat.pop.combi.rob$av_freq)+4, labels = fabian.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(14, 118, 
     labels = expression(paste(italic(rho), ' = 0.65***')),cex = 1.4)

plot(remo.te2.stat.pop.combi.rob$FreqAv, remo.te2.stat.pop.combi.rob$av_freq,
     xlab="Frequency (Remolina)", ylab="Frequency (SA)", xlim= c(0,1), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
abline(0,1, lty = 2)

plot(rank(remo.te2.stat.pop.combi.rob$FreqAv), rank(remo.te2.stat.pop.combi.rob$av_freq),xlim= c(0,118), ylim = c(0,118),
     xlab="Frequency Ranks (Remolina)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
abline(0,1, lty = 2)
text(rank(remo.te2.stat.pop.combi.rob$FreqAv), rank(remo.te2.stat.pop.combi.rob$av_freq)+4, 
     labels = remo.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(14, 118, 
     labels = expression(paste(italic(rho), ' = 0.58***')),cex = 1.4)

plot(hoed.te2.stat.pop.combi.rob$FreqAv, hoed.te2.stat.pop.combi.rob$av_freq,
     xlab="Frequency (Hoedjes)", ylab="Frequency (SA)", xlim= c(0,1), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
abline(0,1, lty = 2)

plot(rank(hoed.te2.stat.pop.combi.rob$FreqAv), rank(hoed.te2.stat.pop.combi.rob$av_freq),xlim= c(0,118), ylim = c(0,118),
     xlab="Frequency Ranks (Hoedjes)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
abline(0,1, lty = 2)
text(rank(hoed.te2.stat.pop.combi.rob$FreqAv), rank(hoed.te2.stat.pop.combi.rob$av_freq)+4, 
     labels = hoed.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(14, 118, 
     labels = expression(paste(italic(rho), ' = 0.61***')),cex = 1.4)

dev.off()

###
#Plot of correlation between genomic copy number in experimental evolution studies and Kofler2015 frequencies
carnes.te2.stat.pop.combi.rob$MeanIns = rowMeans(carnes.te2.stat.pop.combi.rob[,c("Mean_Cont","Mean_Sel")])
fabian.te2.stat.pop.combi.rob$MeanIns = rowMeans(fabian.te2.stat.pop.combi.rob[,c("Mean_Cont","Mean_Sel")])
remo.te2.stat.pop.combi.rob$MeanIns = rowMeans(remo.te2.stat.pop.combi.rob[,c("Mean_Cont","Mean_Sel")])
hoed.te2.stat.pop.combi.rob$MeanIns = rowMeans(hoed.te2.stat.pop.combi.rob[,c("Mean_Early","Mean_Postponed")])

pdf(paste(output_path,"FigS8_correlation_RobertFreq_vs_CopyNumb.pdf",sep=""), width=12, height=20)
par(mfrow=c(4,2),font.axis=1,font.lab=2,lwd=1.2,cex=1,mar=c(4.5,4.5,2,3))
cor.test(carnes.te2.stat.pop.combi.rob$MeanIns,carnes.te2.stat.pop.combi.rob$av_freq,method="spearman") 
#r = -0.529, p 5.409e-09
plot(carnes.te2.stat.pop.combi.rob$MeanIns, carnes.te2.stat.pop.combi.rob$av_freq, xlim=c(0,140),ylim = c(0,1),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Insertions (Carnes)", ylab="Frequency (SA)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)

plot(rank(carnes.te2.stat.pop.combi.rob$MeanIns), rank(carnes.te2.stat.pop.combi.rob$av_freq), xlim= c(0,118), ylim = c(0,118),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Insertion Ranks (Carnes)", ylab="Frequency Ranks (SA)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)
text(rank(carnes.te2.stat.pop.combi.rob$MeanIns), rank(carnes.te2.stat.pop.combi.rob$av_freq)+2, labels = carnes.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.53***')),cex = 1.4)

cor.test(fabian.te2.stat.pop.combi.rob$MeanIns, fabian.te2.stat.pop.combi.rob$av_freq,method="spearman")
#r = -0.4634679 , p 6.395e-07
plot(fabian.te2.stat.pop.combi.rob$MeanIns, fabian.te2.stat.pop.combi.rob$av_freq, 
     xlab="Insertions (Fabian)", ylab="Frequency (SA)", xlim= c(0,170), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)

plot(rank(fabian.te2.stat.pop.combi.rob$MeanIns), rank(fabian.te2.stat.pop.combi.rob$av_freq), xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Fabian)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)
text(rank(fabian.te2.stat.pop.combi.rob$MeanIns), rank(fabian.te2.stat.pop.combi.rob$av_freq)+4, labels = fabian.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.46***')),cex = 1.4)

cor.test(remo.te2.stat.pop.combi.rob$MeanIns, remo.te2.stat.pop.combi.rob$av_freq,method="spearman")
#r = -0.5122253 , p 2.331e-08
plot(remo.te2.stat.pop.combi.rob$MeanIns, remo.te2.stat.pop.combi.rob$av_freq,
     xlab="Insertions (Remolina)", ylab="Frequency (SA)", xlim= c(0,140), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)

plot(rank(remo.te2.stat.pop.combi.rob$MeanIns), rank(remo.te2.stat.pop.combi.rob$av_freq),xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Remolina)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)
text(rank(remo.te2.stat.pop.combi.rob$MeanIns), rank(remo.te2.stat.pop.combi.rob$av_freq)+4, 
     labels = remo.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.51***')),cex = 1.4)

cor.test(hoed.te2.stat.pop.combi.rob$MeanIns, hoed.te2.stat.pop.combi.rob$av_freq,method="spearman")
#r = -0.4386907 , p 1.829e-06
plot(hoed.te2.stat.pop.combi.rob$MeanIns, hoed.te2.stat.pop.combi.rob$av_freq,
     xlab="Insertions (Hoedjes)", ylab="Frequency (SA)", xlim= c(0,140), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)

plot(rank(hoed.te2.stat.pop.combi.rob$MeanIns), rank(hoed.te2.stat.pop.combi.rob$av_freq),xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Hoedjes)", ylab="Frequency Ranks (SA)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)
text(rank(hoed.te2.stat.pop.combi.rob$MeanIns), rank(hoed.te2.stat.pop.combi.rob$av_freq)+4, 
     labels = hoed.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.43***')),cex = 1.4)

dev.off()

###
#Plot of correlation between genomic copy number and TE family frequency in experimental evolution studies
pdf(paste(output_path,"FigS9_correlation_OurFreq_vs_CopyNumb.pdf",sep=""), width=12, height=20)
par(mfrow=c(4,2),font.axis=1,font.lab=2,lwd=1.2,cex=1,mar=c(4.5,4.5,2,3))
cor.test(carnes.te2.stat.pop.combi.rob$MeanIns,carnes.te2.stat.pop.combi.rob$FreqAv,method="spearman") 
#p 0.54, r = -0.06874435
plot(carnes.te2.stat.pop.combi.rob$MeanIns, carnes.te2.stat.pop.combi.rob$FreqAv, xlim=c(0,140),ylim = c(0,1),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Insertions (Carnes)", ylab="Frequency (Carnes)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)

plot(rank(carnes.te2.stat.pop.combi.rob$MeanIns), rank(carnes.te2.stat.pop.combi.rob$FreqAv), xlim= c(0,118), ylim = c(0,118),
     cex.lab = 1.8,cex.axis = 1.3,
     xlab="Insertion Ranks (Carnes)", ylab="Frequency Ranks (Carnes)",
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),carnes.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)
text(rank(carnes.te2.stat.pop.combi.rob$MeanIns), rank(carnes.te2.stat.pop.combi.rob$FreqAv)+2, labels = carnes.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.07 (ns)')),cex = 1.4)

cor.test(fabian.te2.stat.pop.combi.rob$MeanIns, fabian.te2.stat.pop.combi.rob$FreqAv,method="spearman")
#-0.36 , p 0.0001511
plot(fabian.te2.stat.pop.combi.rob$MeanIns, fabian.te2.stat.pop.combi.rob$FreqAv, 
     xlab="Insertions (Fabian)", ylab="Frequency (Fabian)", xlim= c(0,170), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)

plot(rank(fabian.te2.stat.pop.combi.rob$MeanIns), rank(fabian.te2.stat.pop.combi.rob$FreqAv), xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Fabian)", ylab="Frequency Ranks (Fabian)", cex.lab = 1.8,cex.axis = 1.3,
     pch = 21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),fabian.te2.stat.pop.combi.rob$Type.x))
#abline(0,1,lty = 2)
text(rank(fabian.te2.stat.pop.combi.rob$MeanIns), rank(fabian.te2.stat.pop.combi.rob$FreqAv)+4, labels = fabian.te2.stat.pop.combi.rob$flybase_name, cex = 0.6, col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.36**')),cex = 1.4)

cor.test(remo.te2.stat.pop.combi.rob$MeanIns, remo.te2.stat.pop.combi.rob$FreqAv,method="spearman")
#-0.407671 , p 1.422e-05
plot(remo.te2.stat.pop.combi.rob$MeanIns, remo.te2.stat.pop.combi.rob$FreqAv,
     xlab="Insertions (Remolina)", ylab="Frequency (Remolina)", xlim= c(0,140), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)

plot(rank(remo.te2.stat.pop.combi.rob$MeanIns), rank(remo.te2.stat.pop.combi.rob$FreqAv),xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Remolina)", ylab="Frequency Ranks (Remolina)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),remo.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)
text(rank(remo.te2.stat.pop.combi.rob$MeanIns), rank(remo.te2.stat.pop.combi.rob$FreqAv)+4, 
     labels = remo.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.41***')),cex = 1.4)

cor.test(hoed.te2.stat.pop.combi.rob$MeanIns, hoed.te2.stat.pop.combi.rob$FreqAv,method="spearman")
#-0.38, p  4.021e-05
plot(hoed.te2.stat.pop.combi.rob$MeanIns, hoed.te2.stat.pop.combi.rob$FreqAv,
     xlab="Insertions (Hoedjes)", ylab="Frequency (Hoedjes)", xlim= c(0,140), ylim = c(0,1), cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)

plot(rank(hoed.te2.stat.pop.combi.rob$MeanIns), rank(hoed.te2.stat.pop.combi.rob$FreqAv),xlim= c(0,118), ylim = c(0,118),
     xlab="Insertion Ranks (Hoedjes)", ylab="Frequency Ranks (Hoedjes)", cex.lab = 1.8,cex.axis = 1.3,
     pch=21, bg = mgsub(c("ns","C>S","S>C"), c("grey","blue","red"),hoed.te2.stat.pop.combi.rob$Type.x))
#abline(0,1, lty = 2)
text(rank(hoed.te2.stat.pop.combi.rob$MeanIns), rank(hoed.te2.stat.pop.combi.rob$FreqAv)+4, 
     labels = hoed.te2.stat.pop.combi.rob$flybase_name, cex = 0.6,col = "lightgrey")
text(100, 118, 
     labels = expression(paste(italic(rho), ' = -0.38***')),cex = 1.4)

dev.off()



#################################################################################
#COMPARISON OF AVERAGE FREQUENCIES BETWEEN S>C AND C>S FROM SA POPULATION
#################################################################################
#Carnes
carnes.te2.stat.pop.combi.sign.rob$plotcol = mgsub(c("TRUE","FALSE"),c("red","blue"),carnes.te2.stat.pop.combi.sign.rob$Diff_SelCont > 0)
carnes.te2.stat.pop.combi.sign.rob$plotpch = 21
carnes.te2.stat.pop.combi.sign.rob = carnes.te2.stat.pop.combi.sign.rob[!is.na(carnes.te2.stat.pop.combi.sign.rob$av_freq),]

pdf(paste(output_path,"Fig2_Carnes_AvFreq_BoxPlots_log2FC.pdf",sep=""), width=8, height=6)
m = matrix(c(1,1,2,2),nrow = 1,ncol = 4,byrow = TRUE)
layout(mat = m)
par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(carnes.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(carnes.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(carnes.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(carnes.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),carnes.te2.stat.pop.combi.sign.rob$Type.x)),factor=0.5), 
       carnes.te2.stat.pop.combi.sign.rob$av_freq, 
       pch = carnes.te2.stat.pop.combi.sign.rob$plotpch, 
       bg= carnes.te2.stat.pop.combi.sign.rob$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)

title(ylab = "Frequency",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
t.test(carnes.te2.stat.pop.combi.sign.rob[carnes.te2.stat.pop.combi.sign.rob$Type.x == "S>C",]$av_freq,
       carnes.te2.stat.pop.combi.sign.rob[carnes.te2.stat.pop.combi.sign.rob$Type.x == "C>S",]$av_freq)$p.value #0.0008277877

text(1.5,1.14,"***",cex=1.8)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)

par(mar=c(4,0,1,5),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
top10.c = headtail(carnes.te2.stat.pop.combi.sign.rob[order(carnes.te2.stat.pop.combi.sign.rob$log2RatSelCont),],10)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),top10.c$Type.x)),factor=0.5), 
       top10.c$av_freq, 
       pch = top10.c$plotpch, 
       bg= top10.c$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
title(xlab = "Top 10", line = 2.2,cex.lab=1.4)
t.test(top10.c$av_freq ~ top10.c$Type.x)$p.value #0.026
text(1.5,1.14,"*",cex=1.8)
brackets(x1 = 1, y1 = 1.05,
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()

#Fabian
fabian.te2.stat.pop.combi.sign.rob$plotcol = mgsub(c("TRUE","FALSE"),c("red","blue"),fabian.te2.stat.pop.combi.sign.rob$Diff_SelCont > 0)
fabian.te2.stat.pop.combi.sign.rob$plotpch = 21
fabian.te2.stat.pop.combi.sign.rob = fabian.te2.stat.pop.combi.sign.rob[!is.na(fabian.te2.stat.pop.combi.sign.rob$av_freq),]

pdf(paste(output_path,"Fig2_Fabian_AvFreq_BoxPlots_log2FC.pdf",sep=""), width=8, height=6)
m = matrix(c(1,1,2,2),nrow = 1,ncol = 4,byrow = TRUE)
layout(mat = m)
par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(fabian.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(fabian.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(fabian.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(fabian.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),fabian.te2.stat.pop.combi.sign.rob$Type.x)),factor=0.5), 
       fabian.te2.stat.pop.combi.sign.rob$av_freq, 
       pch = fabian.te2.stat.pop.combi.sign.rob$plotpch, 
       bg= fabian.te2.stat.pop.combi.sign.rob$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)
# axis(4,seq(0,1,0.2),labels = F, line = 0,mgp=c(3, .5, 0))
title(ylab = "Frequency",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
t.test(fabian.te2.stat.pop.combi.sign.rob[fabian.te2.stat.pop.combi.sign.rob$Type.x == "S>C",]$av_freq,
       fabian.te2.stat.pop.combi.sign.rob[fabian.te2.stat.pop.combi.sign.rob$Type.x == "C>S",]$av_freq)$p.value #0.000120587

text(1.5,1.14,"***",cex=1.8)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)

par(mar=c(4,0,1,5),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
top10.c = headtail(fabian.te2.stat.pop.combi.sign.rob[order(fabian.te2.stat.pop.combi.sign.rob$log2RatSelCont),],10)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),top10.c$Type.x)),factor=0.5), 
       top10.c$av_freq, 
       pch = top10.c$plotpch, 
       bg= top10.c$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
title(xlab = "Top 10", line = 2.2,cex.lab=1.4)
t.test(top10.c$av_freq ~ top10.c$Type.x)$p.value #0.77
text(1.5,1.14,"ns",cex=1.4)
brackets(x1 = 1, y1 = 1.05,
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()

#Remolina
remo.te2.stat.pop.combi.sign.rob$plotcol = mgsub(c("TRUE","FALSE"),c("red","blue"),remo.te2.stat.pop.combi.sign.rob$Diff_SelCont > 0)
remo.te2.stat.pop.combi.sign.rob$plotpch = 21
remo.te2.stat.pop.combi.sign.rob = remo.te2.stat.pop.combi.sign.rob[!is.na(remo.te2.stat.pop.combi.sign.rob$av_freq),]

pdf(paste(output_path,"Fig2_Remolina_AvFreq_BoxPlots_log2FC.pdf",sep=""), width=8, height=6)
m = matrix(c(1,1,2,2),nrow = 1,ncol = 4,byrow = TRUE)
layout(mat = m)

par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(remo.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(remo.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(remo.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(remo.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),remo.te2.stat.pop.combi.sign.rob$Type.x)),factor=0.5), 
       remo.te2.stat.pop.combi.sign.rob$av_freq, 
       pch = remo.te2.stat.pop.combi.sign.rob$plotpch, 
       bg= remo.te2.stat.pop.combi.sign.rob$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)
# axis(4,seq(0,1,0.2),labels = F, line = 0,mgp=c(3, .5, 0))
title(ylab = "Frequency",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
t.test(remo.te2.stat.pop.combi.sign.rob[remo.te2.stat.pop.combi.sign.rob$Type.x == "S>C",]$av_freq,
       remo.te2.stat.pop.combi.sign.rob[remo.te2.stat.pop.combi.sign.rob$Type.x == "C>S",]$av_freq)$p.value #3.846508e-05

text(1.5,1.14,"***",cex=1.8)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)

par(mar=c(4,0,1,5),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
top10.c = headtail(remo.te2.stat.pop.combi.sign.rob[order(remo.te2.stat.pop.combi.sign.rob$log2RatSelCont),],10)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),top10.c$Type.x)),factor=0.5), 
       top10.c$av_freq, 
       pch = top10.c$plotpch, 
       bg= top10.c$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
title(xlab = "Top 10", line = 2.2,cex.lab=1.4)
t.test(top10.c$av_freq ~ top10.c$Type.x)$p.value #0.23
text(1.5,1.14,"ns",cex=1.4)
brackets(x1 = 1, y1 = 1.05,
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()

#Hoedjes
hoed.te2.stat.pop.combi.sign.rob$plotcol = mgsub(c("TRUE","FALSE"),c("red","blue"),hoed.te2.stat.pop.combi.sign.rob$Diff_PostEarly > 0)
hoed.te2.stat.pop.combi.sign.rob$plotpch = 21
hoed.te2.stat.pop.combi.sign.rob = hoed.te2.stat.pop.combi.sign.rob[!is.na(hoed.te2.stat.pop.combi.sign.rob$av_freq),]

pdf(paste(output_path,"Fig2_Hoedjes_AvFreq_BoxPlots_log2FC.pdf",sep=""), width=8, height=6)
m = matrix(c(1,1,2,2),nrow = 1,ncol = 4,byrow = TRUE)
layout(mat = m)

par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(hoed.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(hoed.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(hoed.te2.stat.pop.combi.sign.rob$av_freq ~ as.factor(hoed.te2.stat.pop.combi.sign.rob$Type.x), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),hoed.te2.stat.pop.combi.sign.rob$Type.x)),factor=0.5), 
       hoed.te2.stat.pop.combi.sign.rob$av_freq, 
       pch = hoed.te2.stat.pop.combi.sign.rob$plotpch, 
       bg= hoed.te2.stat.pop.combi.sign.rob$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)
# axis(4,seq(0,1,0.2),labels = F, line = 0,mgp=c(3, .5, 0))
title(ylab = "Frequency",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
t.test(hoed.te2.stat.pop.combi.sign.rob[hoed.te2.stat.pop.combi.sign.rob$Type.x == "S>C",]$av_freq,
       hoed.te2.stat.pop.combi.sign.rob[hoed.te2.stat.pop.combi.sign.rob$Type.x == "C>S",]$av_freq)$p.value #0.004043846

text(1.5,1.14,"**",cex=1.8)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)

par(mar=c(4,0,1,5),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
top10.c = headtail(hoed.te2.stat.pop.combi.sign.rob[order(hoed.te2.stat.pop.combi.sign.rob$log2Rat_PostEarly),],10)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(top10.c$av_freq ~ as.factor(top10.c$Type.x),
        outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C"),c(1,2),top10.c$Type.x)),factor=0.5), 
       top10.c$av_freq, 
       pch = top10.c$plotpch, 
       bg= top10.c$plotcol,cex=1.4)
axis(1,seq(1,2,1),labels = c("C>S","S>C"),cex.axis=1.2)
title(xlab = "Top 10", line = 2.2,cex.lab=1.4)
t.test(top10.c$av_freq ~ top10.c$Type.x)$p.value #0.52
text(1.5,1.14,"ns",cex=1.4)
brackets(x1 = 1, y1 = 1.05,
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()


################################################################################
#IDENTIFY TE INSERTIONS WITH SIGNIFICANT FREQUENCY DIFFERENCES
################################################################################
#stats all insertions
p.val.vect2 = NULL
for (i in 1:nrow(fabian.te2.filt)){
  int.tab = as.data.frame(t(fabian.te2.filt[i,c(9:14)]))
  int.tab$Reg = c("C","C",rep("S",4))
  model = aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg)
  p.val = summary(model)[[1]][["Pr(>F)"]][1]
  Df = summary(model)[[1]][["Df"]][1]
  Resid = summary(model)[[1]][["Df"]][2]
  F.val = summary(model)[[1]][["F value"]][1]
  p.val.vect1 = c(Df, Resid, F.val,p.val)
  p.val.vect2 = rbind(p.val.vect2,p.val.vect1)
}
p.val.vect2
fabian.te2.filt$Df = p.val.vect2[,1]
fabian.te2.filt$Resid = p.val.vect2[,2]
fabian.te2.filt$Fval = p.val.vect2[,3]
fabian.te2.filt$Pval = p.val.vect2[,4]
#fabian.te2.filt$FDR = p.adjust(fabian.te2.filt$Pval,"fdr")
bonf.005 = 0.05 / length(fabian.te2.filt$Pval)
table(fabian.te2.filt$Pval < bonf.005) #38 significant
fabian.te2.filt$Bonf005 = fabian.te2.filt$Pval < bonf.005
table(fabian.te2.filt$Bonf005)
sum(na.omit(fabian.te2.filt[fabian.te2.filt$Pval < bonf.005,])$dFreq > 0) #8
sum(na.omit(fabian.te2.filt[fabian.te2.filt$Pval < bonf.005,])$dFreq < 0) #30
#fabian.te2.filt.signBonf = na.omit(fabian.te2.filt[fabian.te2.filt$Pval < bonf.005,])

fabian.te2.filt.edit = fabian.te2.filt[,-c(1,4,8,15:20)]
fabian.te2.filt.edit1 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
      fabian.te2.filt.edit,by="TEfam")
      
write.table(fabian.te2.filt.edit1, paste(output_path,"Fabian2018_StatsPerInsertion.txt",sep=""),quote=F,sep="\t",row.names = F)

#manhattan plot Fabian
fabian.te2.filt.test = fabian.te2.filt
fabian.te2.filt.test = fabian.te2.filt.test[order(fabian.te2.filt.test$Pval,decreasing = F),]
fabian.te2.filt.test[1,c("Pval")] = 2e-16
#fabian.te2.filt.test[1,25] = 2e-16

names(fabian.te2.filt.test)[c(2,3,28)] = c("CHR","BP","P")
fabian.te2.filt.test$CHR
fabian.te2.filt.test$CHR = factor(fabian.te2.filt.test$CHR, levels = c("X", "2L", "2R", "3L", "3R","4"))
fabian.te2.filt.test$SNP = paste(fabian.te2.filt.test$CHR,fabian.te2.filt.test$BP,fabian.te2.filt.test$TEfam,sep = ".")
# fabian.te2.filt.test$SNP = fabian.te2.filt.test$TEfam
mysnps = fabian.te2.filt.test[fabian.te2.filt.test$dFreq >0 & fabian.te2.filt.test$P < bonf.005,]$SNP
length(mysnps)
mysnps2 = fabian.te2.filt.test[fabian.te2.filt.test$dFreq <0 & fabian.te2.filt.test$P < bonf.005,]$SNP
length(mysnps2)
mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
mypalette <- c("black", "grey60")
sig = bonf.005 # significant threshold line
sugg = bonf.005 # suggestive threshold line

#Some functions used by the manhattan plot function might be in use already
#Probably need to detach packages or restart R. Examples:
detach("package:biomaRt", unload=TRUE)
detach("package:topGO", unload=TRUE)
detach("package:org.Dm.eg.db", unload=TRUE)

fabian.manh = gg.manhattan(df = fabian.te2.filt.test, 
                           threshold=bonf.005,
                           hlight=mysnps, 
                           hlight2=mysnps2,
                           col=mypalette, 
                           ylims=c(0,17), title="Fabian2018")
fabian.manh

ggsave(fabian.manh, file=paste(output_path,"Fig3_Manhattan_plot_Bonf005_correct1.pdf",sep=""), height=7, width=10)

#HOEDJES
#stats all insertions
p.val.vect2 = NULL
for (i in 1:nrow(hoed.te2.filt)){
  int.tab = as.data.frame(t(hoed.te2.filt[i,c(9:32)]))
  
  int.tab$Type = gsub("[1-9]","",row.names(int.tab))
  int.tab$Reg = mgsub(c("CE","CP","HE","HP","LE","LP"),c("C","S","C","S","C","S"),int.tab$Type)
  int.tab$Diet = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Control","Control","High","High","Low","Low"),int.tab$Type)
  model = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Diet + int.tab$Reg + int.tab$Diet * int.tab$Reg))
 
  Df.val.diet = model[[1]][["Df"]][1]
  Df.val.regime = model[[1]][["Df"]][2]
  Df.val.inter = model[[1]][["Df"]][3]
  resid = model[[1]][["Df"]][4]
  
  F.val.diet = model[[1]][["F value"]][1]
  F.val.regime = model[[1]][["F value"]][2]
  F.val.inter = model[[1]][["F value"]][3]
  
  p.val.diet = model[[1]][["Pr(>F)"]][1]
  p.val.regime = model[[1]][["Pr(>F)"]][2]
  p.val.inter = model[[1]][["Pr(>F)"]][3]
  
  stat.val = c(Df.val.diet, F.val.diet,p.val.diet,
               Df.val.regime, F.val.regime,p.val.regime,
               Df.val.inter,F.val.inter, p.val.inter,
               resid)
  p.val.vect2 = rbind(p.val.vect2,stat.val)
}
p.val.vect2 = as.data.frame(p.val.vect2)
names(p.val.vect2) = c('Df.diet', 'Fval.diet','Pval.diet',
                       'Df.regime', 'Fval.regime','Pval.regime',
                       'Df.inter','Fval.inter','Pval.inter','ResidDf')

hoed.te2.filt.new = cbind(hoed.te2.filt, p.val.vect2)
hoed.te2.filt.new = as.data.frame(hoed.te2.filt.new)

bonf.005 = 0.05 / length(hoed.te2.filt.new$Pval.regime)
hoed.te2.filt.new$Bonf005.regime = hoed.te2.filt.new$Pval.regime < bonf.005
hoed.te2.filt.new$Bonf005.diet = hoed.te2.filt.new$Pval.diet < bonf.005
hoed.te2.filt.new$Bonf005.inter = hoed.te2.filt.new$Pval.inter < bonf.005
table(hoed.te2.filt.new$Bonf005.regime) #100
table(hoed.te2.filt.new$Bonf005.diet) #137
table(hoed.te2.filt.new$Bonf005.inter) #37

sum(hoed.te2.filt.new[hoed.te2.filt.new$Bonf005.regime == TRUE,]$dFreq < 0) #44
sum(hoed.te2.filt.new[hoed.te2.filt.new$Bonf005.regime == TRUE,]$dFreq > 0) #56

hoed.te2.filt.new1 = hoed.te2.filt.new[,-c(1,4,8,33:56)]
hoed.te2.filt.new2 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
      hoed.te2.filt.new1,by="TEfam")

write.table(hoed.te2.filt.new2,paste(output_path,"Hoedjes2019_StatsPerInsertion.txt",sep=""),quote=F,sep="\t",row.names = F)

#manhattan REGIME, HOEDJES2019
hoed.te2.filt.test = hoed.te2.filt.new
hoed.te2.filt.test.reg = hoed.te2.filt.test[,c("Chrom","Pos","Pval.regime","Bonf005.regime","dFreq","TEfam")]

names(hoed.te2.filt.test.reg) = c("CHR","BP","P","Bonf","dFreq","TEfam")
hoed.te2.filt.test.reg$CHR
hoed.te2.filt.test.reg$CHR = factor(hoed.te2.filt.test.reg$CHR, levels = c("X", "2L", "2R", "3L", "3R","4"))
hoed.te2.filt.test.reg$SNP = paste(hoed.te2.filt.test.reg$CHR,hoed.te2.filt.test.reg$BP,hoed.te2.filt.test.reg$TEfam,sep = ".")
# hoed.te2.filt.test$SNP = hoed.te2.filt.test$TEfam
mysnps = hoed.te2.filt.test.reg[hoed.te2.filt.test.reg$dFreq >0 & hoed.te2.filt.test.reg$P < bonf.005,]$SNP
length(mysnps) #56
mysnps.cont = hoed.te2.filt.test.reg[hoed.te2.filt.test.reg$dFreq <0 & hoed.te2.filt.test.reg$P < bonf.005,]$SNP
length(mysnps.cont) #44

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
mypalette <- c("black", "grey60")

sig = bonf.005 # significant threshold line
sugg = bonf.005 # suggestive threshold line

detach("package:topGO", unload=TRUE)
detach("package:org.Dm.eg.db", unload=TRUE)
library(dplyr)
hoed.manh = gg.manhattan(df = hoed.te2.filt.test.reg, 
                         threshold=max(hoed.te2.filt.test.reg$P), 
                         hlight=mysnps,
                         hlight2 = mysnps.cont,
                         col=mypalette, 
                         ylims=c(0,17), title="Hoedjes2019")
hoed.manh

ggsave(hoed.manh, file=paste(output_path,"Fig3_Manhattan_Hoedjes2019_Regime_Bonf005_1.pdf",sep=""), height=7, width=10)


#REMOLINA2012
#stats all insertions
p.val.vect2 = NULL
for (i in 1:nrow(remo.te2.filt)){
  int.tab = as.data.frame(t(remo.te2.filt[i,c(9:14)]))
  int.tab$Reg = c(rep("C",3),rep("S",3))
  Df = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Df"]][1]
  Resid = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Df"]][2]
  F.val = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["F value"]][1]
  p.val = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Pr(>F)"]][1]
  stat.vect = c(Df,Resid,F.val,p.val)
  p.val.vect2 = rbind(p.val.vect2,stat.vect)
}
p.val.vect2

remo.te2.filt$Df= p.val.vect2[,1]
remo.te2.filt$ResidDf = p.val.vect2[,2]
remo.te2.filt$Fval = p.val.vect2[,3]
remo.te2.filt$Pval = p.val.vect2[,4]

table(p.adjust(p.val.vect2[,4],"fdr") < 0.05) #0 significant
table(p.val.vect2[,4] < 0.05/length(p.val.vect2[,4])) #0 significant
remo.te2.filt$Bonf005 = remo.te2.filt$Pval < 0.05/length(remo.te2.filt$Pval)

remo.te2.filt.new = remo.te2.filt[,-c(1,4,8,15:20)]
remo.te2.filt.new2 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                           remo.te2.filt.new,by="TEfam")
write.table(remo.te2.filt.new2,paste(output_path,"Remolina2012_StatsPerInsertion.txt",sep=""),quote=F,sep="\t",row.names = F)

#CARNES2015
#stats all insertions
p.val.vect2 = NULL
for (i in 1:nrow(carnes.te2.filt)){
  int.tab = as.data.frame(t(carnes.te2.filt[i,c(9:18)]))
  int.tab$Reg = c(rep("C",5),rep("S",5))
  Df = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Df"]][1]
  ResidDf = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Df"]][2]
  F.val = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["F value"]][1]
  p.val = summary(aov(asin(sqrt(int.tab[,1])) ~ int.tab$Reg))[[1]][["Pr(>F)"]][1]
  stat.val = c(Df, ResidDf, F.val,p.val)
  p.val.vect2 = rbind(p.val.vect2,stat.val)
}
p.val.vect2

carnes.te2.filt$Df = p.val.vect2[,1]
carnes.te2.filt$ResidDf = p.val.vect2[,2]
carnes.te2.filt$Fval = p.val.vect2[,3]
carnes.te2.filt$Pval = p.val.vect2[,4]

table(p.adjust(p.val.vect2[,4],"fdr") < 0.05) #0 significant
table(p.val.vect2[,4] < 0.05/length(p.val.vect2[,4])) #0 significant
carnes.te2.filt$Bonf005 = carnes.te2.filt$Pval < (0.05/length(carnes.te2.filt$Pval))

carnes.te2.filt.new = carnes.te2.filt[,-c(1,4,8,19:28)]
carnes.te2.filt.new2 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                            carnes.te2.filt.new,by="TEfam")
write.table(carnes.te2.filt.new2,paste(output_path,"Carnes2015_StatsPerInsertion.txt",sep=""),quote=F,sep="\t",row.names = F)

#################################################################
#COMPARE IF TE FAMILY FREQUENCY VARIES BETWEEN S AND C POPULATIONS
#################################################################

#Collapse PopTE2 Table so that each position in own row
hoed.edit = melt(hoed.te2.filt[,c(2,3,5,9:32)],id.vars = c("Chrom","Pos","TEfam"),variable.name = "Pop",value.name = "Freq")
fabian.edit = melt(fabian.te2.filt[,c(2,3,5,9:14)],id.vars = c("Chrom","Pos","TEfam"),variable.name = "Pop",value.name = "Freq")
remo.edit = melt(remo.te2.filt[,c(2,3,5,9:14)],id.vars = c("Chrom","Pos","TEfam"),variable.name = "Pop",value.name = "Freq")
carnes.edit = melt(carnes.te2.filt[,c(2,3,5,9:18)],id.vars = c("Chrom","Pos","TEfam"),variable.name = "Pop",value.name = "Freq")

#Edit tables
fabian.edit$Regime = mgsub(c("Ra","Rb","2La","2Lb","La","Lb"),c("Cont","Cont","Sel","Sel","Sel","Sel"),fabian.edit$Pop)
remo.edit$Regime = mgsub(c("C1","C2","C3","S1","S2","S3"),c("Cont","Cont","Cont","Sel","Sel","Sel"),remo.edit$Pop)
carnes.edit$Regime = mgsub(c("C1","C2","C3","C4","C5","S1","S2","S3","S4","S5"),c("Cont","Cont","Cont","Cont","Cont","Sel","Sel","Sel","Sel","Sel"),carnes.edit$Pop)
Type = gsub("[1-9]","",hoed.edit$Pop)
unique(Type) #"CE" "CP" "HE" "HP" "LE" "LP"
Breeding = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Early","Postponed","Early","Postponed","Early","Postponed"),Type)
unique(Breeding)
Diet = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Control","Control","High","High","Low","Low"),Type)
unique(Diet)
hoed.edit$Type = Type
hoed.edit$Breeding = Breeding
hoed.edit$Diet = Diet
head(hoed.edit)

#Do the TEs that differ signficiantly in abundance also differ significantly in frequency?
#HOEDJES 2019
te.fam.hoed = unique(hoed.edit$TEfam)
stat_output = NULL
for (i in 1:length(te.fam.hoed)) {
  te.subset = as.character(te.fam.hoed[i])
  tab.subset = hoed.edit[hoed.edit$TEfam == te.subset,]
  
  #averages of types (i.e. 6 regimes)
  mean_Breeding = tapply(tab.subset$Freq, INDEX = tab.subset$Breeding, mean)
  diff_PostEarly =  mean_Breeding[2] - mean_Breeding[1]
  diff_PostEarly
  
  log2rat_PostEarly =  log2(mean_Breeding[2] / mean_Breeding[1])
  log2rat_PostEarly
  
  model.lm = lm(asin(sqrt(Freq)) ~ Diet + Breeding + Diet * Breeding, data=tab.subset) 
  anova(model.lm)
  
  Df_Diet = anova(model.lm)[["Df"]][1]
  Df_Breeding = anova(model.lm)[["Df"]][2]
  Df_inter = anova(model.lm)[["Df"]][3]
  Df_Resid = anova(model.lm)[["Df"]][4]
  
  Fval_Diet = anova(model.lm)[["F value"]][1]
  Fval_Breeding = anova(model.lm)[["F value"]][2]
  Fval_inter = anova(model.lm)[["F value"]][3]
  
  pvalFull_Diet = anova(model.lm)[["Pr(>F)"]][1]
  pvalFull_Breeding = anova(model.lm)[["Pr(>F)"]][2]
  pvalFull_inter = anova(model.lm)[["Pr(>F)"]][3]
  
  int.tab = cbind(te.subset, 
                  round(matrix(mean_Breeding,ncol=2),3), #Early Postponed 
                  round(diff_PostEarly,3),
                  round(log2rat_PostEarly,3),
                  Df_Breeding, Fval_Breeding, pvalFull_Breeding,
                  Df_Diet, Fval_Diet, pvalFull_Diet,
                  Df_inter, Fval_inter, pvalFull_inter,
                  Df_Resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","Mean_Early","Mean_Postponed","diff_Freq","log2diff_Freq","Df_Breeding","Fval_Breeding","pvalFull_Breeding",'Df_Diet', 'Fval_Diet','pvalFull_Diet','Df_inter', 'Fval_inter',"pvalFull_inter","Df_Resid")
for (i in 2:ncol(stat_output)){
  stat_output[,i] = numfact(stat_output[,i])
}
stat_output = stat_output[order(stat_output$pvalFull_Breeding, decreasing = F),]
stat_output$fdrFull_Breeding = p.adjust(stat_output$pvalFull_Breeding,method="fdr")
head(stat_output,20 )
table(stat_output$fdrFull_Breeding<0.05) #3

stat_output.sign = stat_output[stat_output$fdrFull_Breeding<0.05,]
stat_output.sign.merge = merge(stat_output.sign,hoed.stat.sign,by="TEfam")
stat_output.sign.merge[,c('flybase_name',"Mean_Postponed.x" ,"Mean_Early.x","Diff_PostEarly" ,'Mean_Postponed.y',"Mean_Early.y","diff_Freq")]
#hoed.stat.filt.sign[hoed.stat.filt.sign$TEfam %in% stat_output.sign$TEfam,][,c('TEfam','Diff_PostEarly','flybase_name')]

#roo, jockey, P-element more abundant and more frequency in selected than in controls!
stat_output1 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                     stat_output,by="TEfam")

write.table(stat_output1,paste(output_path,"TableS11_Hoedjes_Frequency_diff_per_family.txt",sep=""),quote = F, row.names=F, sep="\t")

#FABIAN 2018
te.fam.fabian = unique(fabian.edit$TEfam)
stat_output = NULL
for (i in 1:length(te.fam.fabian)) {
  te.subset = as.character(te.fam.fabian[i])
  tab.subset = fabian.edit[fabian.edit$TEfam == te.subset,]
  #averages of types (i.e. 6 regimes)
  mean_Regime = tapply(tab.subset$Freq, INDEX = tab.subset$Regime, mean)
  diff_SelCont =  mean_Regime[2] - mean_Regime[1]
  diff_SelCont
  
  log2rat_SelCont =  log2(mean_Regime[2] / mean_Regime[1])
  log2rat_SelCont
  
  model.lm = lm(asin(sqrt(Freq)) ~ Regime + Regime/Pop, data=tab.subset) 
  anova(model.lm)
  
  Df = anova(model.lm)[["Df"]][1]
  Fvalue = anova(model.lm)[["F value"]][1]
  pvalue = anova(model.lm)[["Pr(>F)"]][1]
  
  Df.nested = anova(model.lm)[["Df"]][2]
  Fvalue.nested = anova(model.lm)[["F value"]][2]
  pvalue.nested = anova(model.lm)[["Pr(>F)"]][2]
  
  Df_Resid = Df.nested = anova(model.lm)[["Df"]][3]
  
  int.tab = cbind(te.subset, 
                  round(matrix(mean_Regime,ncol=2),3), #Early Postponed 
                  round(diff_SelCont,3),
                  round(log2rat_SelCont,3), 
                  Df, Fvalue, pvalue,
                  Df.nested, Fvalue.nested, pvalue.nested,
                  Df_Resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","Mean_Cont","Mean_Sel","diff_Freq","log2diff_Freq","Df_Regime","Fvalue_Regime","pval_Regime","Df_nested","Fvalue_nested","pval_nested","Df_Resid")
for (i in 2:ncol(stat_output)){
  stat_output[,i] = numfact(stat_output[,i])
}
stat_output = stat_output[order(stat_output$pval_Regime, decreasing = F),]
stat_output$fdr_Regime = p.adjust(stat_output$pval_Regime,method="fdr")
head(stat_output,20 )
table(stat_output$fdr_Regime<0.05) #1

stat_output.sign = stat_output[stat_output$fdr_Regime<0.05,]
stat_output.sign.merge = merge(stat_output.sign,fabian.stat.sign,by="TEfam")
stat_output.sign.merge[,c('flybase_name','Mean_Sel.y' ,'Mean_Cont.y',"Diff_SelCont",'Mean_Sel.x','Mean_Cont.x',"diff_Freq")] 
#P-element more abundance and more frequency in controls

stat_output1 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                     stat_output,by="TEfam")
write.table(stat_output1,paste(output_path,"TableS11_Fabian_Frequency_diff_per_family.txt",sep=""),quote = F, row.names=F, sep="\t")


#CARNES 2015
te.fam.carnes = unique(carnes.edit$TEfam)
stat_output = NULL
for (i in 1:length(te.fam.carnes)) {
  te.subset = as.character(te.fam.carnes[i])
  tab.subset = carnes.edit[carnes.edit$TEfam == te.subset,]
  #print(nrow(tab.subset))
  #averages of types (i.e. 6 regimes)
  mean_Regime = tapply(tab.subset$Freq, INDEX = tab.subset$Regime, mean)
  diff_SelCont =  mean_Regime[2] - mean_Regime[1]
  diff_SelCont
  
  log2rat_SelCont =  log2(mean_Regime[2] / mean_Regime[1])
  log2rat_SelCont
  
  model.lm = lm(asin(sqrt(Freq)) ~ Regime + Regime/Pop, data=tab.subset) 
  anova(model.lm)
  
  Df = anova(model.lm)[["Df"]][1]
  Fvalue = anova(model.lm)[["F value"]][1]
  pvalue = anova(model.lm)[["Pr(>F)"]][1]
  
  Df.nested = anova(model.lm)[["Df"]][2]
  Fvalue.nested = anova(model.lm)[["F value"]][2]
  pvalue.nested = anova(model.lm)[["Pr(>F)"]][2]
  
  Df_Resid = anova(model.lm)[["Df"]][3]
    
  int.tab = cbind(te.subset, 
                  round(matrix(mean_Regime,ncol=2),3), #Early Postponed 
                  round(diff_SelCont,3),
                  round(log2rat_SelCont,3), 
                  Df, Fvalue, pvalue,
                  Df.nested, Fvalue.nested, pvalue.nested,
                  Df_Resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","Mean_Cont","Mean_Sel","diff_Freq","log2diff_Freq","Df_Regime","Fvalue_Regime","pval_Regime","Df_nested","Fvalue_nested","pval_nested","Df_Resid")
for (i in 2:ncol(stat_output)){
  stat_output[,i] = numfact(stat_output[,i])
}
stat_output = stat_output[order(stat_output$pval_Regime, decreasing = F),]
stat_output$fdr_Regime = p.adjust(stat_output$pval_Regime,method="fdr")
head(stat_output,20 )
table(stat_output$fdr_Regime<0.05) #39
stat_output.sign = stat_output[stat_output$fdr_Regime<0.05,]
stat_output.sign
stat_output.sign.merge = merge(stat_output.sign,carnes.stat.sign,by="TEfam")
stat_output.sign.merge[,c('flybase_name','Mean_Sel.y' ,'Mean_Cont.y',"Diff_SelCont",'Mean_Sel.x','Mean_Cont.x',"diff_Freq")]
nrow(stat_output.sign.merge) #34

table(stat_output.sign.merge$diff_Freq >0 & stat_output.sign.merge$Diff_SelCont >0) #27: higher frequency and higher abundance in S populations 
table(stat_output.sign.merge$diff_Freq >0 & stat_output.sign.merge$Diff_SelCont <0) #7 higher frequency, but lower abundance in S populations 
table(stat_output.sign.merge$diff_Freq <0 & stat_output.sign.merge$Diff_SelCont <0) #none with higher frequency and higehr abundacne in C populations

stat_output1 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                     stat_output,by="TEfam")
write.table(stat_output1,paste(output_path,"TableS11_Carnes_Frequency_diff_per_family.txt",sep=""),quote = F, row.names=F, sep="\t")


#REMOLINA 2012
te.fam.remo = unique(remo.edit$TEfam)
stat_output = NULL
for (i in 1:length(te.fam.remo)) {
  te.subset = as.character(te.fam.remo[i])
  tab.subset = remo.edit[remo.edit$TEfam == te.subset,]
  #print(nrow(tab.subset))
  #averages of types (i.e. 6 regimes)
  mean_Regime = tapply(tab.subset$Freq, INDEX = tab.subset$Regime, mean)
  diff_SelCont =  mean_Regime[2] - mean_Regime[1]
  diff_SelCont
  
  log2rat_SelCont =  log2(mean_Regime[2] / mean_Regime[1])
  log2rat_SelCont
  
  model.lm = lm(asin(sqrt(Freq)) ~ Regime + Regime/Pop, data=tab.subset) 
  anova(model.lm)
  Df = anova(model.lm)[["Df"]][1]
  Fvalue = anova(model.lm)[["F value"]][1]
  pvalue = anova(model.lm)[["Pr(>F)"]][1]
  
  Df.nested = anova(model.lm)[["Df"]][2]
  Fvalue.nested = anova(model.lm)[["F value"]][2]
  pvalue.nested = anova(model.lm)[["Pr(>F)"]][2]
  
  Df_Resid = anova(model.lm)[["Df"]][3]
  
  int.tab = cbind(te.subset, 
                  round(matrix(mean_Regime,ncol=2),3), #Early Postponed 
                  round(diff_SelCont,3),
                  round(log2rat_SelCont,3), 
                  Df, Fvalue, pvalue,
                  Df.nested, Fvalue.nested, pvalue.nested,
                  Df_Resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","Mean_Cont","Mean_Sel","diff_Freq","log2diff_Freq","Df_Regime","Fvalue_Regime","pval_Regime","Df_nested","Fvalue_nested","pval_nested", "Df_Resid")
for (i in 2:ncol(stat_output)){
  stat_output[,i] = numfact(stat_output[,i])
}
stat_output = stat_output[order(stat_output$pval_Regime, decreasing = F),]
stat_output$fdr_Regime = p.adjust(stat_output$pval_Regime,method="fdr")
head(stat_output,20 )
table(stat_output$fdr_Regime<0.05) #0 signficant 

stat_output1 = merge(ensembl_casey.dmel[,c("TEfam","flybase_name","fb_id","TEsubclass", "TEclass")],
                     stat_output,by="TEfam")
write.table(stat_output1,paste(output_path,"TableS11_Remolina_Frequency_diff_per_family.txt",sep=""),quote = F, row.names=F, sep="\t")

