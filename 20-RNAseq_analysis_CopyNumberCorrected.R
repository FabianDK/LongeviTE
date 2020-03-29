#Differential expression analysis of filtered RNA-seq read counts from Carnes et al. 2015 (PLoS One)
#Here: normalized by copy number in the genome to get average expression per TE insertion
#For Table S18

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(DESeq2)

setwd("/Users/danfab/Carnes/RNAseq/")

#Load tables
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
TE.dmel = ensembl_casey[ensembl_casey$species == "Dmel",]
TE.dmel$TEfam = droplevels(TE.dmel$TEfam)
TE.dmel$flybase_name = gsub("Dmel/","",TE.dmel[,c(4)],fixed = T) 
nrow(TE.dmel)

#Abundance
TE.detect = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T,check.names = F)

#transpose for DEseq2 and exclude 5 non-Dmel TEs
tab.mean = read.table("/Users/danfab/Carnes/RNAseq/rawCounts_filtered_Chrom_min400_allDmel_mean.txt",header=T,check.names = F)
all.TEs = colnames(tab.mean)[grep(pattern = "FBgn",x = colnames(tab.mean),invert = T)][-c(1:4)]
non.Dmel.TEs = all.TEs[!all.TEs %in% TE.dmel$TEfam] #"DSRN"      "TARTYAK"   "DNTOMRETA" "DDBARI1"   "DTEII"   
tab.mean[,colnames(tab.mean) %in% non.Dmel.TEs]
ncol(tab.mean)
tab.mean = tab.mean[,!colnames(tab.mean) %in% non.Dmel.TEs]
ncol(tab.mean) #13389

tab.mean.t = as.data.frame(t(tab.mean[,5:ncol(tab.mean)]))
class(tab.mean.t$V2)
names(tab.mean.t) = paste(tab.mean$Regime,
                          tab.mean$Pop,
                          tab.mean$Age,
                          tab.mean$Sex,sep="_")
tab.mean.t.round = round(tab.mean.t)
table(row.names(tab.mean.t.round) %in% setdiff(ensembl_casey$TEfam,TE.detect$TEfam))

################################################################################
#################### CORRECTION BY DIVISION, BUT INCLUDING ALL TEs (even without abundance info)! ####################

#Divide TEs by respective copies and genes by 1
bighead(tab.mean,1,10,1,10)
names(tab.mean)

tab.mean.TEs = tab.mean[,colnames(tab.mean) %in% TE.detect$TEfam] #selecte detected  TEs
row.names(tab.mean.TEs) = paste(tab.mean$Regime,tab.mean$Pop,tab.mean$Age,tab.mean$Sex,sep="_")

tab.mean.Genes = tab.mean[,!colnames(tab.mean) %in% TE.detect$TEfam] #select genes this includes TEs without abundance info
row.names(tab.mean.Genes) = paste(tab.mean$Regime,tab.mean$Pop,tab.mean$Age,tab.mean$Sex,sep="_")

names(tab.mean.Genes[,grep(pattern = "FBgn",x = colnames(tab.mean.Genes),invert = T)])[-c(1:4)]
length(names(tab.mean.Genes[,grep(pattern = "FBgn",x = colnames(tab.mean.Genes),invert = T)])[-c(1:4)]) #13

head(TE.detect)
length(TE.detect$TEfam) #112

#select TEs with abundance data also in expression data
TE.detect.new = TE.detect[TE.detect$TEfam %in% names(tab.mean.TEs),] 

nrow(TE.detect.new)
#110
nrow(TE.detect)
#112
ncol(tab.mean.Genes) #13279

#select mean insertion values only
TE.detect.new1 = as.data.frame(t(TE.detect.new[,c(1,3:12)]))
colnames(TE.detect.new1) = TE.detect.new$TEfam
TE.detect.new1 = TE.detect.new1[-1,]
for (i in 1:ncol(TE.detect.new1)){
  TE.detect.new1[,i] = numfact(TE.detect.new1[,i])
}
head(TE.detect.new1)
head(tab.mean.TEs)

#Sort expression tab by same columns as abundance
tab.mean.TEs.sort = tab.mean.TEs[names(TE.detect.new1)]
names(tab.mean.TEs.sort)
names(TE.detect.new1)
nrow(TE.detect.new1) #10

#Create matrix for division with TE abundance estiamtes: rbind 4 times as we have young/old, males/females for which abundance will be used
TE.detect.new1.rep4 = rbind(TE.detect.new1,
                            TE.detect.new1,
                            TE.detect.new1,
                            TE.detect.new1)

#Divide expression counts by copy number estimates
tab.mean.TEs.sort.corrected = tab.mean.TEs.sort / TE.detect.new1.rep4
tab.mean.TEs.sort[,1] / TE.detect.new1.rep4[,1]
tab.mean.TEs.sort.corrected[,1]

tab.mean.TEs.sort[,100] / TE.detect.new1.rep4[,100]
tab.mean.TEs.sort.corrected[,100]

head(tab.mean.TEs.sort)
head(TE.detect.new1)

#merge corrected counts with genes table (genes assumed to be present once)
tab.mean.corrected1 = cbind(tab.mean$Regime,tab.mean$Pop,tab.mean$Age,tab.mean$Sex,
                            tab.mean.TEs.sort.corrected,
                            tab.mean.Genes[,-c(1:4)])
bighead(tab.mean.corrected1,1,10,1,10)
names(tab.mean.corrected1)[1:4] = c("Regime","Pop","Age","Sex")

for (i in 5:ncol(tab.mean.corrected1)){
  tab.mean.corrected1[,i] = numfact(tab.mean.corrected1[,i])
}

ncol(tab.mean.corrected1) #13389
nrow(tab.mean.corrected1) #40
setdiff(colnames(tab.mean),colnames(tab.mean.corrected1))
setdiff(colnames(tab.mean.corrected1), colnames(tab.mean))

########################################
#DESeq2 models on corrected counts (female data only)
########################################

tab.mean.corrected.t = as.data.frame(t(tab.mean.corrected1[,5:ncol(tab.mean.corrected1)]))

names(tab.mean.corrected.t) = paste(tab.mean.corrected1$Regime,
                                    tab.mean.corrected1$Pop,
                                    tab.mean.corrected1$Age,
                                    tab.mean.corrected1$Sex,sep="_")
tab.mean.corrected.t.round = round(tab.mean.corrected.t)

#Subset to females
tab.mean.corrected1.F = tab.mean.corrected1[tab.mean.corrected1$Sex == "F",]
tab.mean.corrected1.F.t = as.data.frame(t(tab.mean.corrected1.F[,5:ncol(tab.mean.corrected1)]))

names(tab.mean.corrected1.F.t) = paste(tab.mean.corrected1.F$Regime,
                                       tab.mean.corrected1.F$Pop,
                                       tab.mean.corrected1.F$Age,
                                       tab.mean.corrected1.F$Sex,sep="_")

ncol(tab.mean.corrected1.F.t)
tab.mean.corrected1.F.t.round = round(tab.mean.corrected1.F.t)

bighead(tab.mean.corrected1.F.t,1,10,1,10)
bighead(tab.mean.corrected1.F.t.round,1,10,1,10)

#DEseq analysis
variables.F = tab.mean.corrected1.F[,1:3]
variables.F$Regime
variables.F$Age = relevel(variables.F$Age, ref = "young")
variables.F$Age

dds.F = DESeqDataSetFromMatrix(countData = tab.mean.corrected1.F.t.round, 
                               colData = variables.F, 
                               design = ~ Age + Regime + Regime*Age)

sizeFactors(dds1.F)
dds1.F = DESeq(dds.F)
# sizeFactors(dds1.F) = rep(1, length(sizeFactors(dds1.F))) #SET THIS, BECAUSE WE ALREADY USED COUNTS NORMALIZED FOR LIBRARY SIZE
# sizeFactors(dds1.F)
res.TE.summary.F = deseq.results(DEseqobj = dds1.F, 
                                 TEvect = ensembl_casey$TEfam,
                                 cutoff = 0.05,
                                 type = "summary")
res.TE.summary.F
# Factor Padj005  ns FC>0 FC<0 Total
# 1   Age_old_vs_young      13 109    7    6   122
# 2 Regime_Sel_vs_Cont      67  56    1   66   123
# 3   Ageold.RegimeSel       4 116    4    0   120

res.TE.list.F = deseq.results(DEseqobj = dds1.F, 
                              TEvect = ensembl_casey$TEfam,
                              cutoff = 0.05, 
                              type="table")

res.TE.list.F.sign = lapply(res.TE.list.F, function(x) na.omit(x[x$padj<0.05,]))


#Genes as well
res.Genes.table.F = deseq.results(DEseqobj = dds1.F, 
                                  TEvect = row.names(tab.mean.corrected1.F.t.round)[!row.names(tab.mean.corrected1.F.t.round) %in% ensembl_casey$TEfam],
                                  cutoff = 0.05,
                                  type = "table")
res.Genes.summary.F= deseq.results(DEseqobj = dds1.F, 
                                   TEvect = row.names(tab.mean.corrected1.F.t.round)[!row.names(tab.mean.corrected1.F.t.round) %in% ensembl_casey$TEfam],
                                   cutoff = 0.05,
                                   type = "summary")




#res.Genes.summary.M
res.Genes.summary.F
# Factor Padj005    ns FC>0 FC<0 Total
# Age_old_vs_young    3335  8287 1679 1656 11622
# Regime_Sel_vs_Cont  2519 10392 1175 1344 12911
# Ageold.RegimeSel    1140  9710  558  582 10850

res.Genes.table.F.sign = lapply(res.Genes.table.F, function(x) na.omit(x[x$padj<0.05,]))

write.table(as.data.frame(res.TE.list.F),file="TE_expr_females_INTERACTION_correctedByDivision_allTEs.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.table.F),file="Genes_expr_females_INTERACTION_correctedByDivision_allTEs.txt",quote=F,sep="\t")

#MAIN MODEL
dds.F.main = DESeqDataSetFromMatrix(countData = tab.mean.corrected1.F.t.round, 
                               colData = variables.F, 
                               design = ~ Age + Regime)

#dds1.M = DESeq(dds.M)
dds1.F.main = DESeq(dds.F.main)

res.TE.summary.F.main = deseq.results(DEseqobj = dds1.F.main, 
                                 TEvect = ensembl_casey$TEfam,
                                 cutoff = 0.05,
                                 type = "summary")
res.TE.summary.F.main
# Factor Padj005  ns FC>0 FC<0 Total
# 1   Age_old_vs_young      11 111    9    2   122
# 2 Regime_Sel_vs_Cont      66  56    3   63   122

res.Genes.summary.F.main = deseq.results(DEseqobj = dds1.F.main, 
                                   TEvect = row.names(tab.mean.corrected1.F.t.round)[!row.names(tab.mean.corrected1.F.t.round) %in% ensembl_casey$TEfam],
                                   cutoff = 0.05,
                                   type = "summary")

res.Genes.summary.F.main
# Factor Padj005   ns FC>0 FC<0 Total
# 1   Age_old_vs_young    2017 9604  929 1088 11621
# 2 Regime_Sel_vs_Cont    2445 9692 1226 1219 12137


res.TE.list.F.main = deseq.results(DEseqobj = dds1.F.main, 
                                      TEvect = ensembl_casey$TEfam,
                                      cutoff = 0.05,
                                      type = "table")
res.Genes.list.F.main = deseq.results(DEseqobj = dds1.F.main, 
                                         TEvect = row.names(tab.mean.corrected1.F.t.round)[!row.names(tab.mean.corrected1.F.t.round) %in% ensembl_casey$TEfam],
                                         cutoff = 0.05,
                                         type = "table")

write.table(as.data.frame(res.TE.list.F.main),file="TE_expr_females_MAIN_FACTORS_correctedByDivision_allTEs.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list.F.main),file="Genes_expr_females_MAIN_FACTORS__correctedByDivision_allTEs.txt",quote=F,sep="\t")
