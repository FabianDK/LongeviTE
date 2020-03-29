#Plot coverage of 5 single copy genes used per default for normalization of TE family coverage by DeviaTE

#For Figure S2

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

dir.create("/Users/danfab/single_copy_gene_coverage")
setwd("/Users/danfab/single_copy_gene_coverage")

#Load tables with coverage of single copy genes
#Files on dryad
fab = read.table("fabian_all_Dmel_rawCov_edited.txt")
carn = read.table("carnes_all_Dmel_rawCov_edited.txt")
remo = read.table("remolina_all_Dmel_rawCov_edited.txt")
hoed = read.table("hoedjes_all_Dmel_rawCov_edited.txt")

names(fab) = c("ALL","Rep","pos","cov","physcov","hq_cov")
names(carn) = c("ALL","Rep","pos","cov","physcov","hq_cov")
names(remo) = c("ALL","Rep","pos","cov","physcov","hq_cov")
names(hoed) = c("ALL","Rep","pos","cov","physcov","hq_cov")

cov.list = list(fab,carn,remo,hoed)
names(cov.list) = c("fab","carn","remo","hoed")
cov.list = lapply(cov.list, function(x) { x$hq_physcov <- x$hq_cov+x$physcov; return(x)})
cov.list.mean = lapply(cov.list, function(x) tapply(x$hq_physcov,list(x$ALL,x$Rep),mean))

mean(cov.list.mean$fab) #240.6853
mean(cov.list.mean$remo) #47.67849
mean(cov.list.mean$carn) #26.4474
mean(cov.list.mean$hoed) #122.5108

#Plot
cov.list.mean.melt = melt(cov.list.mean,by="row.names")
names(cov.list.mean.melt) = c("Gene","Population","Coverage","Study")
cov.list.mean.melt$Gene = gsub(pattern = "Dmel_","",cov.list.mean.melt$Gene)

#Add regime info
regime1 = unlist(lapply(strsplit(x = as.character(cov.list.mean.melt$Population),split = "_"), function(x) x[1]))
regime2 = gsub("[1-9]","",regime1)
regime3 = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Cont","Sel","Cont","Sel","Cont","Sel"),regime2)
cov.list.mean.melt$Regime = regime3

#Change rpl32
cov.list.mean.melt$Gene = gsub(pattern = "rpl32","RpL32",cov.list.mean.melt$Gene)

#subset to studies
cov.list.mean.melt.fab = cov.list.mean.melt[cov.list.mean.melt$Study == "fab",]
cov.list.mean.melt.remo = cov.list.mean.melt[cov.list.mean.melt$Study == "remo",]
cov.list.mean.melt.carn = cov.list.mean.melt[cov.list.mean.melt$Study == "carn",]
cov.list.mean.melt.hoed = cov.list.mean.melt[cov.list.mean.melt$Study == "hoed",]

p.carn = ggplot(data = cov.list.mean.melt.carn, aes(x=Population, y=Coverage, group = Population, label = Gene)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + THEME_DEF + 
  labs(y = expression(bold("Average Coverage")), x = expression(bold(""))) +
  geom_jitter(width = 0.15,size=2.8,shape=21,aes(fill = Gene),colour = "black") + scale_fill_manual(values = c("blue","red","green","grey40","magenta")) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=15, colour = "black"), axis.text.y = element_text(size=18, colour = "black"))

p.fab = ggplot(data = cov.list.mean.melt.fab, aes(x=Population, y=Coverage, group = Population, label = Gene)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + THEME_DEF + 
  labs(y = expression(bold("")), x = expression(bold(""))) +
  geom_jitter(width = 0.15,size=2.8,shape=21,aes(fill = Gene),colour = "black") + scale_fill_manual(values = c("blue","red","green","grey40","magenta")) +
  theme(axis.text.x = element_text(size=12, colour = "black"), axis.text.y = element_text(size=18, colour = "black"))

p.fab2 = ggplot(data = cov.list.mean.melt.fab, aes(x=Population, y=Coverage, group = Population, label = Gene)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + THEME_DEF + 
  labs(y = expression(bold("Average Coverage")), x = expression(bold(""))) +
  theme(legend.position = "none") +
  geom_jitter(width = 0.15,size=2.8,shape=21,aes(fill = Gene),colour = "black") + scale_fill_manual(values = c("blue","red","green","grey40","magenta")) +
  theme(axis.text.x = element_text(size=14, colour = "black"), axis.text.y = element_text(size=18, colour = "black"))

cov.list.mean.melt.hoed$Population = as.factor(mgsub(c("CE","CP"),c("ME","MP"),cov.list.mean.melt.hoed$Population))
cov.list.mean.melt.hoed$Population = factor(cov.list.mean.melt.hoed$Population, levels = c("LE1", "LE2" ,"LE3" ,"LE4" ,"LP1" ,"LP2" ,"LP3", "LP4",
                                                      "ME1" ,"ME2", "ME3" ,"ME4" ,"MP1", "MP2", "MP3", "MP4",
                                                      "HE1" ,"HE2", "HE3", "HE4", "HP1", "HP2", "HP3", "HP4"))
p.hoed = ggplot(data = cov.list.mean.melt.hoed, aes(x=Population, y=Coverage, group = Population, label = Gene)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + THEME_DEF + 
  labs(y = expression(bold("Average Coverage")), x = expression(bold(""))) +
  geom_jitter(width = 0.15,size=4.8,shape=21,aes(fill = Gene),colour = "black") + scale_fill_manual(values = c("blue","red","green","grey40","magenta")) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=20, colour = "black"), axis.text.y = element_text(size=28, colour = "black"), axis.title.y = element_text(size=31, colour = "black"))

p.remo = ggplot(data = cov.list.mean.melt.remo, aes(x=Population, y=Coverage, group = Population, label = Gene)) + geom_boxplot(outlier.shape = NA) + 
  theme_bw() + THEME_DEF + 
  labs(y = expression(bold("")), x = expression(bold(""))) +
  geom_jitter(width = 0.15,size=2.8,shape=21,aes(fill = Gene),colour = "black") + scale_fill_manual(values = c("blue","red","green","grey40","magenta")) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=18, colour = "black"))

ggsave(filename = "Coverage_Dmel_Carnes.pdf",plot = p.carn,width = 10, height = 4)
ggsave(filename = "Coverage_Dmel_just_for_legend.pdf",plot = p.fab,width = 6, height = 4)
ggsave(filename = "Coverage_Dmel_Fabian.pdf",plot = p.fab2,width = 6, height = 4)
ggsave(filename = "Coverage_Dmel_Hoedjes.pdf",plot = p.hoed,width = 18, height = 7)
ggsave(filename = "Coverage_Dmel_Remolina.pdf",plot = p.remo,width = 6, height = 4)

