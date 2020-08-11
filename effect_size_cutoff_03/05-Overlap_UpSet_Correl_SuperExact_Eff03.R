#Overlap of S>C and C>S TE families significant according to approach 1, and correlation in fold change between studies

#For Figure S3A, Table S8

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(SuperExactTest)
library(FactoMineR)
library("Hmisc")
library(rstatix) 
#library(UpSetR)
#trace("upset",edit=TRUE)
detach("package:UpSetR",unload=TRUE)
#remove.packages("UpSetR")
#install.packages("UpSetR")
library(UpSetR)

dir.create("/Users/danfab/comp_all/")
output_path = "/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output/"
setwd("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output/")

effcut = 0.3

#TE lookup table
ensembl_casey = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/table_key/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name)

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
names(hoed.stat.filt)[38:39] = c("Diff_SelCont","log2RatSelCont")

#Tables with fitlered for min_coverage and only Bonferroni significant
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE" & abs(carnes.stat.filt$Diff_SelCont) > effcut,]
carnes.stat.SC = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont>0,]
carnes.stat.CS = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont<0,]

fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE" & abs(fabian.stat.filt$Diff_SelCont) > effcut,]
fabian.stat.SC = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont>0,]
fabian.stat.CS = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont<0,]

remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
remo.stat.filt.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE" & abs(remo.stat.filt$Diff_SelCont) > effcut,]
remo.stat.SC = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont>0,]
remo.stat.CS = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont<0,]

hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == TRUE & abs(hoed.stat.filt$Diff_PostEarly) > effcut,]
hoed.stat.SC = hoed.stat.sign[hoed.stat.sign$Diff_PostEarly>0,]
hoed.stat.CS = hoed.stat.sign[hoed.stat.sign$Diff_PostEarly<0,]

#Background 
TE_bg = intersect(hoed.stat.filt$TEfam, 
                  intersect(intersect(remo.stat.filt$TEfam, carnes.stat.filt$TEfam),fabian.stat.filt$TEfam) ) #background of TEs
te.number = length(TE_bg) 
te.number #103
te.fam = as.data.frame(TE_bg) #restrict to those shared in all
names(te.fam) = "TEfam"

#select columns to extract
col.sel = c('TEfam','TEpropCov','Mean_Cont', 'Mean_Sel', 'Diff_SelCont',"log2RatSelCont") 
col.sel.hoed = c('TEfam','TEpropCov','Mean_Early', 'Mean_Post', 'Diff_PostEarly',"log2Rat_PostEarly") 

#subset tables
remo.stat.sub = remo.stat.filt[,c(which(names(remo.stat.filt) %in% col.sel))]
carnes.stat.sub = carnes.stat.filt[,c(which(names(carnes.stat.filt) %in% col.sel))]
fabian.stat.sub = fabian.stat.filt[,c(which(names(fabian.stat.filt) %in% col.sel))]
hoed.stat.sub = hoed.stat.filt[,c(which(names(hoed.stat.filt) %in% col.sel.hoed))]

remo.stat.sign.sub = remo.stat.sign[,c(which(names(remo.stat.sign) %in% col.sel))]
carnes.stat.sign.sub = carnes.stat.sign[,c(which(names(carnes.stat.sign) %in% col.sel))]
fabian.stat.sign.sub = fabian.stat.sign[,c(which(names(fabian.stat.sign) %in% col.sel))]
hoed.stat.sign.sub = hoed.stat.sign[,c(which(names(hoed.stat.sign) %in% col.sel.hoed))]

#merge and subset tabs
remo.carnes.fabian.sub.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                    list(te.fam, remo.stat.sub, carnes.stat.sub, fabian.stat.sub))
names(remo.carnes.fabian.sub.stat) = mgsub(c(".x",".y"),c(".remo",".carn"),names(remo.carnes.fabian.sub.stat))

remo.carnes.fabian.hoed.sub.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                           list(te.fam, remo.carnes.fabian.sub.stat, hoed.stat.sub))
names(remo.carnes.fabian.hoed.sub.stat) = mgsub(c(".x",".y"),c(".fab",".hoed"),names(remo.carnes.fabian.hoed.sub.stat),fixed=T)

#all positive diff i.e. more in selected - OVERLAP STATS
pos.list = list(as.character(carnes.stat.SC$TEfam),
                as.character(remo.stat.SC$TEfam),
                as.character(fabian.stat.SC$TEfam),
                as.character(hoed.stat.SC$TEfam))
names(pos.list) = c("Carnes_pos","Remo_pos","Fabian_pos","Hoedjes_pos")
pos.list.bg = lapply(pos.list, function(x) intersect(x,TE_bg))
stat.pos = supertest(pos.list.bg,length(TE_bg))
summary(stat.pos) #
# Remo_pos & Hoedjes_pos 
# Carnes_pos & Remo_pos & Hoedjes_pos 
as.data.frame(summary(stat.pos))

#all negative diff i.e. more in selected - OVERLAP STATS
neg.list = list(as.character(carnes.stat.CS$TEfam),
                           as.character(remo.stat.CS$TEfam),
                           as.character(fabian.stat.CS$TEfam),
                           as.character(hoed.stat.CS$TEfam))
names(neg.list) = c("Carnes_neg","Remo_neg","Fabian_neg","Hoedjes_neg")
neg.list.bg = lapply(neg.list, function(x) intersect(x,TE_bg))
stat.neg = supertest(neg.list.bg,length(TE_bg))
summary(stat.neg) 
# Remo_neg & Hoedjes_neg  
#Remo_neg & Fabian_neg & Hoedjes_neg
#Carnes_neg & Hoedjes_neg 
#Carnes_neg & Remo_neg & Hoedjes_neg
#Carnes_neg & Remo_neg & Fabian_neg & Hoedjes_neg


#Upset plot
pos.m = as.data.frame(list_to_matrix(lapply(pos.list, function(x) as.character(x))))
neg.m = as.data.frame(list_to_matrix(lapply(neg.list, function(x) as.character(x))))
names(pos.m) = c("Carnes2015", "Remolina2012","Fabian2018","Hoedjes2019")
names(neg.m) = c("Carnes2015", "Remolina2012","Fabian2018","Hoedjes2019")
pos.m = pos.m[,c(1,3,4,2)]
neg.m = neg.m[,c(1,3,4,2)]

pos.m[pos.m$Carnes2015 == 1 & pos.m$Fabian2018 == 1 & pos.m$Hoedjes2019 ==1 & pos.m$Remolina2012 ==1,]
#

neg.m[neg.m$Carnes2015 == 1 & neg.m$Fabian2018 == 1 & neg.m$Hoedjes2019 ==1,]
#DMREPG and G2

#Carnes_neg & Remo_neg & Fabian_neg & Hoedjes_neg DONE
#Carnes_neg & Remo_neg & Hoedjes_neg DONE
#Carnes_neg & Hoedjes_neg DONE
#Remo_neg & Fabian_neg & Hoedjes_neg DONE
# Remo_neg & Hoedjes_neg DONE
pdf(file="Upset_CS_ConsDif_EffCut03.pdf",width=7,height=7)
upset(neg.m, order.by = "degree", keep.order = TRUE, empty.intersections = "on", 
      intersections = list(list("Remolina2012","Hoedjes2019","Fabian2018","Carnes2015"), 
                           list("Carnes2015", "Fabian2018","Hoedjes2019"), 
                           list("Carnes2015", "Fabian2018","Remolina2012"), 
                           list("Carnes2015", "Hoedjes2019","Remolina2012"), 
                           list("Fabian2018", "Hoedjes2019","Remolina2012"), 
                           list("Carnes2015","Fabian2018"), 
                           list("Carnes2015", "Hoedjes2019"), 
                           list("Carnes2015", "Remolina2012"),
                           list("Fabian2018", "Hoedjes2019"),
                           list("Fabian2018", "Remolina2012"),
                           list("Remolina2012","Hoedjes2019"),
                           list("Carnes2015"),
                           list("Fabian2018"),
                           list("Hoedjes2019"),
                           list("Remolina2012")),
      queries = list(list(query = intersects, 
                          params = list("Fabian2018","Remolina2012","Carnes2015","Hoedjes2019"), color = "red", active = T, 
                          query.name = "All"), 
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019","Carnes2015"), color = "red", active = T, 
                          query.name = "3"), 
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019","Fabian2018"), color = "red", active = T, 
                          query.name = "3"),
                     list(query = intersects, 
                          params = list("Carnes2015","Hoedjes2019"), color = "red", active = T, 
                          query.name = "2"),
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019"), color = "red", active = T, 
                          query.name = "2")),
      matrix.color="black", 
      main.bar.color = "black",
      sets.bar.color=c("black"),
      point.size=5, 
      text.scale = 2,
      mainbar.y.label = "C>S Intersection Size",
      mainbar.y.max = 35)
dev.off()

#Carnes_pos & Remo_pos & Hoedjes_pos
#Remo_pos & Hoedjes_pos 
pdf(file="Upset_SC_ConsDif_EffCut03.pdf",width=7,height=7)
upset(pos.m, nsets = 4, keep.order = T,order.by = "degree",
      sets = c("Remolina2012","Hoedjes2019","Fabian2018","Carnes2015"),
      intersections = list(list("Remolina2012","Hoedjes2019","Fabian2018","Carnes2015"), 
                           list("Carnes2015", "Fabian2018","Hoedjes2019"), 
                           list("Carnes2015", "Fabian2018","Remolina2012"), 
                           list("Carnes2015", "Hoedjes2019","Remolina2012"), 
                           list("Fabian2018", "Hoedjes2019","Remolina2012"), 
                           list("Carnes2015","Fabian2018"), 
                           list("Carnes2015", "Hoedjes2019"), 
                           list("Carnes2015", "Remolina2012"),
                           list("Fabian2018", "Hoedjes2019"),
                           list("Fabian2018", "Remolina2012"),
                           list("Remolina2012","Hoedjes2019"),
                           list("Carnes2015"),
                           list("Fabian2018"),
                           list("Hoedjes2019"),
                           list("Remolina2012")),
      queries = list(list(query = intersects, 
                          params = list("Carnes2015","Remolina2012","Hoedjes2019"), color = "red", active = T, 
                          query.name = "3"),
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019"), color = "red", active = T, 
                          query.name = "2")),
      matrix.color="black", 
      main.bar.color = "black",
      sets.bar.color=c("black"),
      point.size=5, 
      text.scale = 2,
      empty.intersections = "on",
      mainbar.y.label = "S>C Intersection Size",
      mainbar.y.max = 35)
dev.off()

#Save overlapping TEs table
neg.m$TEfam = row.names(neg.m)
pos.m$TEfam = row.names(pos.m)
neg.m.merge = merge(neg.m,ensembl_casey[,c("TEfam","flybase_name","TEsubclass","TEclass")],by="TEfam")
pos.m.merge = merge(pos.m,ensembl_casey[,c("TEfam","flybase_name","TEsubclass","TEclass")],by="TEfam")
setdiff(pos.m.merge$TEfam,TE_bg)

write.table(x = neg.m.merge, "shared_C>S_table_Effcut03.txt",quote = F,row.names = F,sep="\t")
write.table(x = pos.m.merge, "shared_S>C_table_Effcut03.txt",quote = F,row.names = F,sep="\t")

#CORRELATION ANALYSIS
colsel = c("TEfam","log2RatSelCont","log2RatSelCont.remo","log2RatSelCont.carn","log2Rat_PostEarly","Diff_SelCont","Diff_SelCont.remo","Diff_SelCont.carn","Diff_PostEarly")
rank.tab.filt = remo.carnes.fabian.hoed.sub.stat[,names(remo.carnes.fabian.hoed.sub.stat) %in% colsel]
names(rank.tab.filt)[6:9] = c("Diff_SelCont.fab", "log2RatSelCont.fab","Diff_SelCont.hoed", "log2RatSelCont.hoed")
head(rank.tab.filt)

#Spearman log2
cor.test(rank.tab.filt$log2RatSelCont.remo,
         rank.tab.filt$log2RatSelCont.carn,method="spearman") #0.04900927 cor, p = 0.623
cor.test(rank.tab.filt$log2RatSelCont.remo,
         rank.tab.filt$log2RatSelCont.fab,method="spearman") #-0.03301196, p = 0.7406
cor.test(rank.tab.filt$log2RatSelCont.fab,
         rank.tab.filt$log2RatSelCont.carn,method="spearman") #-0.1616973, p = 0.1027
cor.test(rank.tab.filt$log2RatSelCont.hoed,
         rank.tab.filt$log2RatSelCont.remo,method="spearman") #Hoedjes to Remolina: 0.2817086  , p = 0.003942
cor.test(rank.tab.filt$log2RatSelCont.hoed,
         rank.tab.filt$log2RatSelCont.fab,method="spearman") #0.02316692 , p = 0.8163
cor.test(rank.tab.filt$log2RatSelCont.hoed,
         rank.tab.filt$log2RatSelCont.carn,method="spearman") #0.162795  cor, p = 0.1004

rcorr(as.matrix(cbind(rank.tab.filt$log2RatSelCont.fab,
                      rank.tab.filt$log2RatSelCont.carn,
                      rank.tab.filt$log2RatSelCont.remo,
                      rank.tab.filt$log2RatSelCont.hoed)),type=c("spearman"))


# #Chek if the 14 that are shared have something in common
# names(ensembl_casey)[1] = "TE"
# ensembl_casey.bg = ensembl_casey[ensembl_casey$TE %in% TE_bg,]
# ensembl_casey.bg$TEsubclass = droplevels(ensembl_casey.bg$TEsubclass)
# ensembl_casey.bg$TE = droplevels(ensembl_casey.bg$TE)
# nrow(ensembl_casey.bg) #103
# 
# #length
# Reduce(intersect,pos.list) %in% ensembl_casey.bg$TE
# table(Reduce(intersect,pos.list) %in% ensembl_casey.bg$TE)
# 
# ensembl_casey.bg$sharedUP = ensembl_casey.bg$TE %in% Reduce(intersect,pos.list)
# ensembl_casey.bg$remo_carnUP = ensembl_casey.bg$TE %in% intersect(pos.list$Carnes_pos,pos.list$Remo_pos)
# ensembl_casey.bg$remo_fabUP = ensembl_casey.bg$TE %in% intersect(pos.list$Fabian_pos,pos.list$Remo_pos)
# ensembl_casey.bg$fab_carnUP = ensembl_casey.bg$TE %in% intersect(pos.list$Carnes_pos,pos.list$Fabian_pos)
# 
# ensembl_casey.bg$sharedDOWN = ensembl_casey.bg$TE %in% Reduce(intersect,neg.list)
# ensembl_casey.bg$remo_carnDOWN = ensembl_casey.bg$TE %in% intersect(neg.list$Carnes_pos,pos.list$Remo_pos)
# ensembl_casey.bg$remo_fabDOWN = ensembl_casey.bg$TE %in% intersect(neg.list$Fabian_pos,pos.list$Remo_pos)
# ensembl_casey.bg$fab_carnDOWN = ensembl_casey.bg$TE %in% intersect(neg.list$Carnes_pos,pos.list$Fabian_pos)
# 
# #enrichment of classes
# #all
# for(i in 1:length(unique(ensembl_casey.bg$TEsubclass)) ){
#   shared_up_te = ensembl_casey.bg[ensembl_casey.bg$sharedUP == "TRUE",]
#   type = names(table(shared_up_te$TEsubclass))[i]
#   m = matrix(c(table(shared_up_te$TEsubclass)[i],
#                sum(table(shared_up_te$TEsubclass)) - table(shared_up_te$TEsubclass)[i],
#                table(ensembl_casey.bg$TEsubclass)[i],
#                sum(table(ensembl_casey.bg$TEsubclass)) - table(ensembl_casey.bg$TEsubclass)[i]),
#              ncol=2,nrow=2)
#   pval.2 = fisher.test(m,alternative="two.sided")$p.value
#   pval.greater = fisher.test(m,alternative="greater")$p.value
#   pval.less = fisher.test(m,alternative="less")$p.value
#   print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
# }
# #no enrichment
# 
# #class
# ensembl_casey.bg$TEclass = droplevels(ensembl_casey.bg$TEclass)
# for(i in 1:length(unique(ensembl_casey.bg$TEclass)) ){
#   shared_up_te = ensembl_casey.bg[ensembl_casey.bg$sharedUP == "TRUE",]
#   type = names(table(shared_up_te$TEclass))[i]
#   m = matrix(c(table(shared_up_te$TEclass)[i],
#                sum(table(shared_up_te$TEclass)) - table(shared_up_te$TEclass)[i],
#                table(ensembl_casey.bg$TEclass)[i],
#                sum(table(ensembl_casey.bg$TEclass)) - table(ensembl_casey.bg$TEclass)[i]),
#              ncol=2,nrow=2)
#   pval.2 = fisher.test(m,alternative="two.sided")$p.value
#   pval.greater = fisher.test(m,alternative="greater")$p.value
#   pval.less = fisher.test(m,alternative="less")$p.value
#   print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
# }
# #no enrichment
# 

