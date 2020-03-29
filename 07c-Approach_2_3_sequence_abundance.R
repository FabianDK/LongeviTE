#Approach 2 and 3 using sequence abundance (i.e. sum of normalized coverage values across all consensus positions) within a TE family
#For Table S5

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)

#Sum normalized coverage values to get 'sequence abundance'
output_path = "/Users/danfab/comp_all/"

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)
ensembl_casey = ensembl_casey[,c("TEfam","flybase_name","fb_id","TEsubclass","TEclass")]

######################################
#Approach 3 (consistent differences): REMOLINA
######################################

tab=read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/remolina_all_TE_edited_filtered.txt",header=T) #NEEDS TO GO ON DRYAD TOO or change code
names(tab) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab)

pos.sum = as.data.frame(aggregate(tab$hqphys_cov, by=list(tab$TEfam, tab$Rep), sum))
names(pos.sum) = c("TEfam","Rep","Sum")
head(pos.sum)
pos.sum.edit = as.data.frame(dcast(pos.sum, TEfam ~ Rep, value.var = c("Sum")))
head(pos.sum.edit)
pos.sum.edit$Mean_Cont = apply(X = pos.sum.edit[,2:4],1,mean)
pos.sum.edit$Mean_Sel = apply(X = pos.sum.edit[,5:7],1,mean)
head(pos.sum.edit)
table(pos.sum.edit$Mean_Sel > pos.sum.edit$Mean_Cont)

sum.stat = NULL
pos.sum.edit.TEfam = unique(as.character(pos.sum.edit$TEfam))
for (i in 1:length(pos.sum.edit.TEfam)){
  subset.tab = pos.sum.edit[pos.sum.edit$TEfam == pos.sum.edit.TEfam[i],]
  ConsistCS = min(subset.tab[,c('Cont_1','Cont_2' ,'Cont_3')]) > max(subset.tab[,c('Sel_1','Sel_2' ,'Sel_3')])
  ConsistSC = max(subset.tab[,c('Cont_1','Cont_2' ,'Cont_3')]) < min(subset.tab[,c('Sel_1','Sel_2' ,'Sel_3')])
  stat.summary = c(pos.sum.edit.TEfam[i],ConsistSC,ConsistCS)
  sum.stat = rbind(sum.stat, stat.summary)
}
sum.stat = as.data.frame(sum.stat)
head(sum.stat)
names(sum.stat) = c("TEfam","ConsistSC","ConsistCS")
table(sum.stat$ConsistSC) #3 S>C
table(sum.stat$ConsistCS) #0 C>S

remo.sum.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                       list(ensembl_casey,pos.sum.edit,sum.stat))
head(remo.sum.stat)
write.table(remo.sum.stat,paste(output_path,"Remolina_SequenceAbundance_SumNormCov.txt",sep=""),quote =F, row.names=F, sep="\t")

######################################
#Approach 3 (consistent differences): FABIAN
######################################
tab=read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/fabian_all_TE_edited_filtered.txt",header=T)
names(tab) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab)

pos.sum = as.data.frame(aggregate(tab$hqphys_cov, by=list(tab$TEfam, tab$Rep), sum))
names(pos.sum) = c("TEfam","Rep","Sum")
head(pos.sum)
pos.sum.edit = as.data.frame(dcast(pos.sum, TEfam ~ Rep, value.var = c("Sum")))
head(pos.sum.edit)
pos.sum.edit$Mean_Cont = apply(X = pos.sum.edit[,2:3],1,mean)
pos.sum.edit$Mean_Sel = apply(X = pos.sum.edit[,4:7],1,mean)
head(pos.sum.edit)

#Stats
table(pos.sum.edit$Mean_Sel > pos.sum.edit$Mean_Cont)

sum.stat = NULL
pos.sum.edit.TEfam = unique(as.character(pos.sum.edit$TEfam))
for (i in 1:length(pos.sum.edit.TEfam)){
  subset.tab = pos.sum.edit[pos.sum.edit$TEfam == pos.sum.edit.TEfam[i],]
  ConsistCS = min(subset.tab[,c('Cont_Ra','Cont_Rb')]) > max(subset.tab[,c('Sel_2La','Sel_2Lb' ,'Sel_La','Sel_Lb')])
  ConsistSC = max(subset.tab[,c('Cont_Ra','Cont_Rb')]) < min(subset.tab[,c('Sel_2La','Sel_2Lb' ,'Sel_La','Sel_Lb')])
  stat.summary = c(pos.sum.edit.TEfam[i],ConsistSC,ConsistCS)
  sum.stat = rbind(sum.stat, stat.summary)
}
sum.stat = as.data.frame(sum.stat)
head(sum.stat)
names(sum.stat) = c("TEfam","ConsistSC","ConsistCS")

table(sum.stat$ConsistSC) #15
table(sum.stat$ConsistCS) #8

fabian.sum.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                       list(ensembl_casey,pos.sum.edit,sum.stat))
head(fabian.sum.stat)
write.table(fabian.sum.stat,paste(output_path,"Fabian_SequenceAbundance_SumNormCov.txt",sep=""),quote =F, row.names=F, sep="\t")

######################################
#Approach 3 (consistent differences): CARNES
######################################
tab=read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/carnes_all_TE_edited_filtered.txt",header=T)
names(tab) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab)

pos.sum = as.data.frame(aggregate(tab$hqphys_cov, by=list(tab$TEfam, tab$Rep), sum))
names(pos.sum) = c("TEfam","Rep","Sum")
head(pos.sum)
pos.sum.edit = as.data.frame(dcast(pos.sum, TEfam ~ Rep, value.var = c("Sum")))
head(pos.sum.edit)
pos.sum.edit$Mean_Cont = apply(X = pos.sum.edit[,2:6],1,mean)
pos.sum.edit$Mean_Sel = apply(X = pos.sum.edit[,7:11],1,mean)
head(pos.sum.edit)

#Stats
table(pos.sum.edit$Mean_Sel > pos.sum.edit$Mean_Cont)

sum.stat = NULL
pos.sum.edit.TEfam = unique(as.character(pos.sum.edit$TEfam))
for (i in 1:length(pos.sum.edit.TEfam)){
  subset.tab = pos.sum.edit[pos.sum.edit$TEfam == pos.sum.edit.TEfam[i],]
  ConsistCS = min(subset.tab[,c('Cont_B1','Cont_B2','Cont_B3','Cont_B4','Cont_B5')]) > max(subset.tab[,c('Sel_O1','Sel_O2','Sel_O3','Sel_O4','Sel_O5')])
  ConsistSC = max(subset.tab[,c('Cont_B1','Cont_B2','Cont_B3','Cont_B4','Cont_B5')]) < min(subset.tab[,c('Sel_O1','Sel_O2','Sel_O3','Sel_O4','Sel_O5')])
  stat.summary = c(pos.sum.edit.TEfam[i],ConsistSC,ConsistCS)
  sum.stat = rbind(sum.stat, stat.summary)
}
sum.stat = as.data.frame(sum.stat)
head(sum.stat)
names(sum.stat) = c("TEfam","ConsistSC","ConsistCS")

table(sum.stat$ConsistSC) #48
table(sum.stat$ConsistCS) #2

carnes.sum.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                         list(ensembl_casey,pos.sum.edit,sum.stat))
head(carnes.sum.stat)
write.table(carnes.sum.stat,paste(output_path,"Carnes_SequenceAbundance_SumNormCov.txt",sep=""),quote =F, row.names=F, sep="\t")

######################################
#Approach 3 (consistent differences): HOEDJES
######################################
tab=read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/hoedjes_all_TE_edited_filtered.txt",header=T)
names(tab) = c("TEfam","Breeding","Diet","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab)

pos.sum = as.data.frame(aggregate(tab$hqphys_cov, by=list(tab$TEfam, tab$Rep), sum))
names(pos.sum) = c("TEfam","Rep","Sum")
head(pos.sum)
pos.sum.edit = as.data.frame(dcast(pos.sum, TEfam ~ Rep, value.var = c("Sum")))
head(pos.sum.edit)

pos.sum.edit.transp = as.data.frame(t(pos.sum.edit))
colnames(pos.sum.edit.transp) = pos.sum.edit$TEfam
pos.sum.edit.transp = pos.sum.edit.transp[-1,]
head(pos.sum.edit.transp)
Rep = row.names(pos.sum.edit.transp)
Type = gsub("[1-9]","",Rep)
unique(Type) #"CE" "CP" "HE" "HP" "LE" "LP"
Breeding = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Early","Postponed","Early","Postponed","Early","Postponed"),Type)
unique(Breeding)
Diet = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Control","Control","High","High","Low","Low"),Type)
unique(Diet)
pos.sum.edit.transp = cbind(Breeding,Diet,Type,Rep,pos.sum.edit.transp)
for (i in 5:ncol(pos.sum.edit.transp)){
  pos.sum.edit.transp[,i] = numfact(pos.sum.edit.transp[,i])
}
head(pos.sum.edit.transp)

sum.stat = NULL
pos.sum.edit.TEfam = unique(as.character(pos.sum.edit$TEfam))
for (i in 1:length(pos.sum.edit.TEfam)){
  subset.tab = pos.sum.edit[pos.sum.edit$TEfam == pos.sum.edit.TEfam[i],]
  ConsistCS_C = min(subset.tab[,c('CE1','CE2','CE3','CE4')]) > max(subset.tab[,c('CP1','CP2','CP3','CP4')])
  ConsistSC_C = max(subset.tab[,c('CE1','CE2','CE3','CE4')]) < min(subset.tab[,c('CP1','CP2','CP3','CP4')])
  
  ConsistCS_H = min(subset.tab[,c('HE1','HE2','HE3','HE4')]) > max(subset.tab[,c('HP1','HP2','HP3','HP4')])
  ConsistSC_H = max(subset.tab[,c('HE1','HE2','HE3','HE4')]) < min(subset.tab[,c('HP1','HP2','HP3','HP4')])
  
  ConsistCS_L = min(subset.tab[,c('LE1','LE2','LE3','LE4')]) > max(subset.tab[,c('LP1','LP2','LP3','LP4')])
  ConsistSC_L = max(subset.tab[,c('LE1','LE2','LE3','LE4')]) < min(subset.tab[,c('LP1','LP2','LP3','LP4')])

  subset.tab2 = pos.sum.edit.transp[,c('Breeding','Diet',pos.sum.edit.TEfam[i])]
  Mean_Sel = mean(subset.tab2[subset.tab2$Breeding == "Postponed",][,3])
  Mean_Cont = mean(subset.tab2[subset.tab2$Breeding == "Early",][,3])

  stat.summary = c(pos.sum.edit.TEfam[i],Mean_Cont,Mean_Sel,
                   ConsistSC_L,ConsistCS_L,ConsistSC_C,ConsistCS_C,ConsistSC_H,ConsistCS_H) #stat$estimate,stat$p.value,
  sum.stat = rbind(sum.stat, stat.summary)
}
sum.stat = as.data.frame(sum.stat)
head(sum.stat)
for (i in 2:3){
  sum.stat[,i] = numfact(sum.stat[,i])
}
names(sum.stat) = c("TEfam","Mean_Cont","Mean_Sel","ConsistSC_L","ConsistCS_L","ConsistSC_C","ConsistCS_C","ConsistSC_H","ConsistCS_H")

table(sum.stat$ConsistSC_L) #43
table(sum.stat$ConsistCS_L) #2

table(sum.stat$ConsistSC_C) #4
table(sum.stat$ConsistCS_C) #0

table(sum.stat$ConsistSC_H) #3
table(sum.stat$ConsistCS_H) #33

hoedjes.sum.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                         list(ensembl_casey,pos.sum.edit,sum.stat))
head(hoedjes.sum.stat)
write.table(hoedjes.sum.stat,paste(output_path,"Hoedjes_SequenceAbundance_SumNormCov.txt",sep=""),quote =F, row.names=F, sep="\t")


######################################
#Approach 2 (all studied combined) 
######################################

tab.carn=read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/carnes_all_TE_edited_filtered.txt",header=T)
names(tab.carn) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")

tab.fab=read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/fabian_all_TE_edited_filtered.txt",header=T)
names(tab.fab) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")

tab.remo=read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/remolina_all_TE_edited_filtered.txt",header=T)
names(tab.remo) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")

tab.hoed=read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/hoedjes_all_TE_edited_filtered.txt",header=T)
names(tab.hoed) = c("TEfam","Breeding","Diet","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")

pos.sum.remo = as.data.frame(aggregate(tab.remo$hqphys_cov, by=list(tab.remo$TEfam, tab.remo$Rep), sum))
names(pos.sum.remo) = c("TEfam","Rep","Sum")
pos.sum.remo$Study = "Remolina"
sample.split = strsplit(x = as.character(pos.sum.remo$Rep),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
pos.sum.remo$Regime = sample.split.mat[,1]
head(pos.sum.remo)

pos.sum.carn = as.data.frame(aggregate(tab.carn$hqphys_cov, by=list(tab.carn$TEfam, tab.carn$Rep), sum))
names(pos.sum.carn) = c("TEfam","Rep","Sum")
pos.sum.carn$Study = "Carnes"
sample.split = strsplit(x = as.character(pos.sum.carn$Rep),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
pos.sum.carn$Regime = sample.split.mat[,1]
head(pos.sum.carn)

pos.sum.fab = as.data.frame(aggregate(tab.fab$hqphys_cov, by=list(tab.fab$TEfam, tab.fab$Rep), sum))
names(pos.sum.fab) = c("TEfam","Rep","Sum")
pos.sum.fab$Study = "Fabian"
sample.split = strsplit(x = as.character(pos.sum.fab$Rep),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
pos.sum.fab$Regime = sample.split.mat[,1]
head(pos.sum.fab)

pos.sum.hoed = as.data.frame(aggregate(tab.hoed$hqphys_cov, by=list(tab.hoed$TEfam, tab.hoed$Rep), sum))
names(pos.sum.hoed) = c("TEfam","Rep","Sum")
pos.sum.hoed$Study = "Hoedjes"
Type = gsub("[1-9]","",pos.sum.hoed$Rep)
pos.sum.hoed$Regime = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Cont","Sel","Cont","Sel","Cont","Sel"),Type)
head(pos.sum.hoed)

#Combine all
pos.sum.all = rbind(pos.sum.carn,pos.sum.fab,pos.sum.remo,pos.sum.hoed)
pos.total = as.data.frame(t(t(tapply(X = pos.sum.all$Sum,pos.sum.all$Rep,sum))))
pos.total$Rep = row.names(pos.total)
names(pos.total)[1] = "TotalSum"
pos.sum.all1 = merge(pos.sum.all,pos.total,by="Rep")
pos.sum.all1$PropSum = pos.sum.all1$Sum / pos.sum.all1$TotalSum
head(pos.sum.all1)

#Subset to shared TEs
shared.TEfam = Reduce(intersect, list(pos.sum.carn$TEfam,pos.sum.fab$TEfam,pos.sum.remo$TEfam,pos.sum.hoed$TEfam))
length(shared.TEfam) #103

combined.stat = NULL
for (i in 1:length(shared.TEfam)){
  TE.sub = shared.TEfam[i]
  sub.tab = pos.sum.all1[pos.sum.all1$TEfam %in% TE.sub,]
  reg.means = tapply(sub.tab$Sum, sub.tab$Regime, mean)
  dif_SelCont = reg.means[2] - reg.means[1]
  log2_SelCont = log2(reg.means[2] /reg.means[1])
  model = aov(asin(sqrt(PropSum)) ~ Study + Regime + Regime * Study, data = sub.tab)
  summary(model)
  df.val = summary(model)[[1]][['Df']][1:3]
  df.resid = summary(model)[[1]][['Df']][4]
  f.val = summary(model)[[1]][['F value']][1:3]
  p.val = summary(model)[[1]][['Pr(>F)']][1:3]
  int = c(TE.sub, reg.means, 
          dif_SelCont,
          log2_SelCont,
          df.val[1],f.val[1],p.val[1],
          df.val[2],f.val[2],p.val[2],
          df.val[3],f.val[3],p.val[3],
          df.resid)
  combined.stat = rbind(combined.stat,int)
}
combined.stat = as.data.frame(combined.stat)
head(combined.stat)
names(combined.stat) = c("TEfam","Mean_Cont","Mean_Sel","Diff_SelCont","log2_SelCont",
                         "Df_Study","F_value_Study","P_value_Study",
                         "Df_Regime","F_value_Regime","P_value_Regime",
                         "Df_Interaction","F_value_Interaction","P_value_Interaction",
                         "Df_resid")
for (i in 2:ncol(combined.stat)){
  combined.stat[,i] = numfact(combined.stat[,i])
}
combined.stat$FDR_Study = combined.stat$P_value_Study
combined.stat$FDR_Regime = combined.stat$P_value_Regime
combined.stat$FDR_Interaction = combined.stat$P_value_Interaction

head(combined.stat)
nrow(combined.stat) #103
combined.stat.sign = combined.stat[combined.stat$FDR_Regime < 0.05,]
nrow(combined.stat.sign) #51 
nrow(combined.stat[combined.stat$FDR_Regime < 0.05,]) #51
nrow(combined.stat[combined.stat$FDR_Study < 0.05,]) #101
nrow(combined.stat[combined.stat$FDR_Interaction < 0.05,]) #68

table(combined.stat.sign$Diff_SelCont > 0) #40 S>C and 11 C>S
round(40/103,2)
round(11/103,2)
round((103-51)/103,2)

combined.stat.annot = merge(ensembl_casey,combined.stat,by="TEfam")
write.table(combined.stat.annot,paste(output_path,"SequenceAbundance_Approach_#2.txt",sep=""),sep="\t",quote=F,row.names = F)

