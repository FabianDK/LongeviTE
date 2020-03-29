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
output_path = "/Users/danfab/comp_all/"

#TE lookup table
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name)

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
names(hoed.stat.filt)[38:39] = c("Diff_SelCont","log2RatSelCont")

#Tables with fitlered for min_coverage and only Bonferroni significant
remo.stat.sign = remo.stat.filt[remo.stat.filt$Bonf == TRUE,]
carnes.stat.sign = carnes.stat.filt[carnes.stat.filt$Bonf == TRUE,]
fabian.stat.sign = fabian.stat.filt[fabian.stat.filt$Bonf == TRUE,]
hoed.stat.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == TRUE,]

#Background 
TE_bg = intersect(hoed.stat.filt$TEfam, 
                  intersect(intersect(remo.stat.filt$TEfam, carnes.stat.filt$TEfam),fabian.stat.filt$TEfam) ) #background of TEs
te.number = length(TE_bg) 
te.number #103
te.fam = as.data.frame(TE_bg) #restrict to those shared in all
names(te.fam) = "TEfam"

#select columns to extract
col.sel = c('TEfam','TEpropCov','Mean_Cont', 'Mean_Sel', 'Diff_SelCont',"log2RatSelCont") 

#subset tables
remo.stat.sub = remo.stat.filt[,c(which(names(remo.stat.filt) %in% col.sel))]
carnes.stat.sub = carnes.stat.filt[,c(which(names(carnes.stat.filt) %in% col.sel))]
fabian.stat.sub = fabian.stat.filt[,c(which(names(fabian.stat.filt) %in% col.sel))]
hoed.stat.sub = hoed.stat.filt[,c(which(names(hoed.stat.filt) %in% col.sel))]

remo.stat.sign.sub = remo.stat.sign[,c(which(names(remo.stat.sign) %in% col.sel))]
carnes.stat.sign.sub = carnes.stat.sign[,c(which(names(carnes.stat.sign) %in% col.sel))]
fabian.stat.sign.sub = fabian.stat.sign[,c(which(names(fabian.stat.sign) %in% col.sel))]
hoed.stat.sign.sub = hoed.stat.sign[,c(which(names(hoed.stat.sign) %in% col.sel))]

#merge and subset tabs
remo.carnes.fabian.sub.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                     list(te.fam, remo.stat.sub, carnes.stat.sub, fabian.stat.sub))
names(remo.carnes.fabian.sub.stat) = mgsub(c(".x",".y"),c(".remo",".carn"),names(remo.carnes.fabian.sub.stat))

remo.carnes.fabian.hoed.sub.stat = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                          list(te.fam, remo.carnes.fabian.sub.stat, hoed.stat.sub))
names(remo.carnes.fabian.hoed.sub.stat) = mgsub(c(".x",".y"),c(".fab",".hoed"),names(remo.carnes.fabian.hoed.sub.stat))

#all positive diff i.e. more in selected - OVERLAP STATS
pos.list = list(carnes.stat.sign.sub[carnes.stat.sign.sub$Diff_SelCont>0,]$TEfam,
                remo.stat.sign.sub[remo.stat.sign.sub$Diff_SelCont>0,]$TEfam,
                fabian.stat.sign.sub[fabian.stat.sign.sub$Diff_SelCont>0,]$TEfam,
                hoed.stat.sign.sub[hoed.stat.sign.sub$Diff_SelCont>0,]$TEfam)
names(pos.list) = c("Carnes_pos","Remo_pos","Fabian_pos","Hoedjes_pos")
pos.list.bg = lapply(pos.list, function(x) intersect(x,TE_bg))
stat.pos = supertest(pos.list.bg,length(TE_bg))
summary(stat.pos) #Remo/Hoedjes: P ~ 0.02460324
as.data.frame(summary(stat.pos))

#all negative diff i.e. more in selected - OVERLAP STATS
neg.list = list(carnes.stat.sign.sub[carnes.stat.sign.sub$Diff_SelCont<0,]$TEfam,
                remo.stat.sign.sub[remo.stat.sign.sub$Diff_SelCont<0,]$TEfam,
                fabian.stat.sign.sub[fabian.stat.sign.sub$Diff_SelCont<0,]$TEfam,
                hoed.stat.sign.sub[hoed.stat.sign.sub$Diff_SelCont<0,]$TEfam)
names(neg.list) = c("Carnes_neg","Remo_neg","Fabian_neg","Hoedjes_neg")
neg.list.bg = lapply(neg.list, function(x) intersect(x,TE_bg))
stat.neg = supertest(neg.list.bg,length(TE_bg))
summary(stat.neg) # Quadruple: 0.01178729 ; Carn/Remo/Hoed: 0.01294287; Remo/Fab/Hoedj: 0.04491782; Remo/Hoed: 0.00440764

#Upset plot
pos.m = as.data.frame(list_to_matrix(lapply(pos.list, function(x) as.character(x))))
neg.m = as.data.frame(list_to_matrix(lapply(neg.list, function(x) as.character(x))))
names(pos.m) = c("Carnes2015", "Remolina2012","Fabian2018","Hoedjes2019")
names(neg.m) = c("Carnes2015", "Remolina2012","Fabian2018","Hoedjes2019")
pos.m = pos.m[,c(1,3,4,2)]
neg.m = neg.m[,c(1,3,4,2)]

pdf(file=paste(output_path,"Fig2_Upset_CS.pdf",sep=""),width=7,height=7)
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
                          params = list("Remolina2012","Hoedjes2019","Fabian2018","Carnes2015"), color = "red", active = T, 
                          query.name = "All"), 
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019","Carnes2015"), color = "red", active = T, 
                          query.name = "3"), 
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019","Fabian2018"), color = "red", active = T, 
                          query.name = "3"),
                     list(query = intersects, 
                          params = list("Remolina2012","Hoedjes2019"), color = "red", active = T, 
                          query.name = "2")),
      matrix.color="black", 
      main.bar.color = "black",
      sets.bar.color=c("black"),
      point.size=5, 
      text.scale = 2,
      mainbar.y.label = "C>S Intersection Size",
      mainbar.y.max = 20)
dev.off()

pdf(file=paste(output_path,"Fig2_Upset_SC.pdf",sep=""),width=7,height=7)
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
                          params = list("Remolina2012","Hoedjes2019"), color = "red", active = T, 
                          query.name = "2")),
      matrix.color="black", 
      main.bar.color = "black",
      sets.bar.color=c("black"),
      point.size=5, 
      text.scale = 2,
      empty.intersections = "on",
      mainbar.y.label = "S>C Intersection Size",
      mainbar.y.max = 20)
dev.off()


neg.m$TEfam = row.names(neg.m)
pos.m$TEfam = row.names(pos.m)
neg.m.merge = merge(neg.m,ensembl_casey[,c("TEfam","flybase_name","TEsubclass","TEclass")],by="TEfam")
pos.m.merge = merge(pos.m,ensembl_casey[,c("TEfam","flybase_name","TEsubclass","TEclass")],by="TEfam")
setdiff(pos.m.merge$TEfam,TE_bg)

write.table(x = neg.m.merge, paste(output_path,"shared_C>S_table.txt",sep=""),quote = F,row.names = F,sep="\t")
write.table(x = pos.m.merge, paste(output_path,"shared_S>C_table.txt",sep=""),quote = F,row.names = F,sep="\t")

#CORRELATION ANALYSIS
colsel = c("TEfam","log2RatSelCont.fab","log2RatSelCont.remo","log2RatSelCont.carn","log2RatSelCont.hoed","Diff_SelCont.fab","Diff_SelCont.remo","Diff_SelCont.carn","Diff_SelCont.hoed")
rank.tab.filt = remo.carnes.fabian.hoed.sub.stat[,names(remo.carnes.fabian.hoed.sub.stat) %in% colsel]
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
