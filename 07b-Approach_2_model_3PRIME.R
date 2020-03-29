#Approach 2: Analysis of TE abundance by combining all studies and normalizing by total TE content
#Insertions calculated considering only the last 200bp of 3'-ends

#For Table S4

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(SuperExactTest)
library(FactoMineR)
library("Hmisc")
library(rstatix) 

dir.create("/Users/danfab/comp_all/")
dir.create("/Users/danfab/comp_all/interact_TotContNorm/negative/",recursive = T)
dir.create("/Users/danfab/comp_all/interact_TotContNorm/positive/",recursive = T)
dir.create("/Users/danfab/comp_all/interaction_plots/positive/",recursive = T)
dir.create("/Users/danfab/comp_all/interaction_plots/negative/",recursive = T)

output_path = "/Users/danfab/comp_all/"

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)

#Total content
total_cont.all = read.table("/Users/danfab/comp_all/Total_TE_content_3PRIME.txt",header=T)
total_cont.all = total_cont.all[,1:6]

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent_3PRIME_200bp.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent_3PRIME_200bp.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent_3PRIME_200bp.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent_3PRIME_200bp.txt",header = T)
names(hoed.stat.filt)[38:39] = c("Diff_SelCont","log2RatSelCont")

#Background across all studies
TE_bg = intersect(hoed.stat.filt$TEfam, 
                  intersect(intersect(remo.stat.filt$TEfam, carnes.stat.filt$TEfam),fabian.stat.filt$TEfam) ) #background of TEs
te.number = length(TE_bg) 
te.number #60
te.fam = as.data.frame(TE_bg) #restrict to those shared in all
names(te.fam) = "TEfam"

################################################################################################
#PREPARE TABLE: Combine all studies to single tab
################################################################################################

#Create tab with mean ins and all studies combined
remo.carnes.fabian.stat.filt = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                     list(te.fam, remo.stat.filt, carnes.stat.filt, fabian.stat.filt))
names(remo.carnes.fabian.stat.filt) = mgsub(c(".x",".y"),c(".remo",".carn"),names(remo.carnes.fabian.stat.filt))
names(remo.carnes.fabian.stat.filt)[3:8] = c("Cont_1.remo", "Cont_2.remo", "Cont_3.remo", "Sel_1.remo",  "Sel_2.remo",  "Sel_3.remo")

remo.carnes.fabian.hoed.stat.filt = Reduce(function(...) merge(..., all = F, by = "TEfam"),
                                      list(te.fam, remo.carnes.fabian.stat.filt, hoed.stat.filt))
names(remo.carnes.fabian.hoed.stat.filt) = mgsub(c(".x",".y"),c(".fab",".hoed"),names(remo.carnes.fabian.hoed.stat.filt))

head(remo.carnes.fabian.hoed.stat.filt)

col.sel1 = c("TEfam","Cont_1.remo","Cont_2.remo","Cont_3.remo","Sel_1.remo","Sel_2.remo","Sel_3.remo",
             "Cont_B1","Cont_B2","Cont_B3","Cont_B4","Cont_B5",
             "Sel_O1","Sel_O2","Sel_O3","Sel_O4","Sel_O5",
             "Cont_Ra", "Cont_Rb","Sel_2La","Sel_2Lb","Sel_La","Sel_Lb",
             names(hoed.stat.filt)[9:32])

remo.carnes.fabian.hoed.stat.filt.means = remo.carnes.fabian.hoed.stat.filt[,c(which(names(remo.carnes.fabian.hoed.stat.filt) %in% col.sel1))]
remo.carnes.fabian.hoed.stat.filt.means = as.data.frame(t(remo.carnes.fabian.hoed.stat.filt[,c(which(names(remo.carnes.fabian.hoed.stat.filt) %in% col.sel1))]),stringsAsFactors=F)
colnames(remo.carnes.fabian.hoed.stat.filt.means) = c(remo.carnes.fabian.hoed.stat.filt.means[1,])
head(remo.carnes.fabian.hoed.stat.filt.means)

remo.carnes.fabian.hoed.stat.filt.means = remo.carnes.fabian.hoed.stat.filt.means[-1,]
remo.carnes.fabian.hoed.stat.filt.means = cbind(rep=rownames(remo.carnes.fabian.hoed.stat.filt.means),remo.carnes.fabian.hoed.stat.filt.means)

breed = mgsub(c('[A-Z]E','[A-Z]P'),c('Cont_','Sel_'),remo.carnes.fabian.hoed.stat.filt.means$rep)
regime = matrix(unlist(strsplit(x = as.character(breed),split = "_")),ncol=2,byrow=T)[,1]
diet = mgsub(c('C[A-Z][1-9]','L[A-Z][1-9]','H[A-Z][1-9]'),c('Cont','Low','High'),remo.carnes.fabian.hoed.stat.filt.means$rep)
diet[1:22] = NA
diet
study = c(rep("Remolina",6),rep("Carnes",10),rep("Fabian",6),rep("Hoedjes",24))

mean_ins_tab = cbind(study,regime,diet,remo.carnes.fabian.hoed.stat.filt.means)
head(mean_ins_tab)
for (i in 5:ncol(mean_ins_tab)){
  mean_ins_tab[,i] = numfact(mean_ins_tab[,i]) #character columns to numbers
}
nrow(mean_ins_tab) #46
ncol(mean_ins_tab) -5 #60
nrow(te.fam) #60

mean_ins_tab$rep = gsub(pattern = ".remo",replacement = "",mean_ins_tab$rep)
names(total_cont.all)[3] = "rep"
mean_ins_tab = merge(mean_ins_tab,total_cont.all,by="rep")
mean_ins_tab$Regime == mean_ins_tab$regime
mean_ins_tab$study == mean_ins_tab$Study
mean_ins_tab = mean_ins_tab[,-c(65:68)]
ncol(mean_ins_tab) #65
mean_ins_tab = mean_ins_tab[,c(65,1:64)]
#write.table(mean_ins_tab,paste(output_path,"Table_Mean_Insertions_Per_Family_With_TotalCont.txt",sep=""),sep="\t",row.names = F,quote=F)

################
#INTERACTION MODEL- NORMALIZED BY TOTAL COUNTS AND ASIN SQRT TRANSFORMED
################

normintmodel_output = NULL
bonf = 0.05 / te.number
for (i in 6:ncol(mean_ins_tab)) {
  int.tab = mean_ins_tab[,c(1:4,i)]
  int.tab[,5] = int.tab[,5] / int.tab$Ins
  int.tab = int.tab[,-1]
  int.av = matrix(tapply(int.tab[,4],list(int.tab$regime),mean),ncol = 2)
  int.diff = int.av[2] - int.av[1]
  int.log2FC = log2(int.av[2] / int.av[1])
  model = lm(asin(sqrt(int.tab[,4])) ~ study + regime + regime*study,data=int.tab) #positive means more in selected -- NEEDS TO BE ASINSQRT TRANSFORMED ?
  reg.stat = matrix(c(anova(model)[["Df"]][2],
                      round(anova(model)[["F value"]][2],2), 
                      anova(model)[["Pr(>F)"]][2]),ncol=3,byrow=T)
  
  int.stat = matrix(c(anova(model)[["Df"]][3],
                      round(anova(model)[["F value"]][3],2), 
                      anova(model)[["Pr(>F)"]][3]),ncol=3,byrow=T)
  
  study.stat = matrix(c(anova(model)[["Df"]][1],
                        round(anova(model)[["F value"]][1],2), 
                        anova(model)[["Pr(>F)"]][1]),ncol=3,byrow=T)
  
  df.resid = anova(model)[["Df"]][4]
  
  int.tab.stat = cbind(names(int.tab)[4], 
                       int.av, int.diff, int.log2FC,
                       reg.stat,reg.stat[,3] <= bonf,
                       int.stat,int.stat[,3] <= bonf,
                       study.stat,study.stat[,3] <= bonf,
                       summary(model)$adj.r.squared,
                       df.resid)
  normintmodel_output = rbind(normintmodel_output,int.tab.stat)
}
normintmodel_output = as.data.frame(normintmodel_output)
head(normintmodel_output)
for (i in c(2:8,10:12,14:16,18:19)) {
  normintmodel_output[,i] = numfact(normintmodel_output[,i])
}
names(normintmodel_output) = c("TEfam","Cont","Sel","DiffSelCont","log2SelCont",
                               "regime_Df","regime_F","regime_P","regime_Bonf",
                               "inter_Df","inter_F","inter_P","inter_Bonf",
                               "study_Df","study_F","study_P","study_Bonf",
                               "adj_R_squared","Df_resid")
head(normintmodel_output)
normintmodel_output$regime_fdr = p.adjust(as.numeric(normintmodel_output$regime_P),method = "fdr")
normintmodel_output$inter_fdr = p.adjust(as.numeric(normintmodel_output$inter_P),method = "fdr")
normintmodel_output$study_fdr = p.adjust(numfact(normintmodel_output$study_P),method = "fdr")
normintmodel_output

#Create output tble
nrow(normintmodel_output) #60
normintmodel_output1 = merge(normintmodel_output,ensembl_casey,by="TEfam",all=F) #merge to casey
normintmodel_output2 = normintmodel_output1[!is.na(normintmodel_output1$regime_F),]
write.table(normintmodel_output2,paste(output_path,"TotContNorm_interaction_model_ASINSQRT_allStat_3PRIME.txt",sep=""),quote = F,row.names = F,sep="\t")

#At Bonferroni alpha 0.05: 
round(table(normintmodel_output[normintmodel_output$regime_Bonf == "TRUE",]$DiffSelCont>0) / te.number,2) #12% pos, 2&

#At FDR 0.05
table(normintmodel_output$regime_fdr < 0.05) #21
table(normintmodel_output$study_fdr < 0.05) #57
table(normintmodel_output$inter_fdr < 0.05) #34
nrow(normintmodel_output) #60
round(table(normintmodel_output[normintmodel_output$regime_fdr < 0.05,]$DiffSelCont>0) / te.number,2) #27% pos, 8%
table(normintmodel_output[normintmodel_output$regime_fdr < 0.05,]$DiffSelCont>0) #16 TEs sign more in selected, 5 sign less
nrow(normintmodel_output[normintmodel_output$regime_fdr >= 0.05,])/ te.number #65% n.s.
nrow(normintmodel_output[normintmodel_output$regime_fdr < 0.05 & normintmodel_output$inter_fdr > 0.05,]) #2 with regime sign, but not interaction
nrow(normintmodel_output[normintmodel_output$regime_fdr < 0.05 & normintmodel_output$inter_fdr < 0.05,]) #19 with both significant
nrow(normintmodel_output[normintmodel_output$inter_fdr < 0.05,]) #34 with significant interaction
nrow(normintmodel_output[normintmodel_output$regime_fdr > 0.05 & normintmodel_output$inter_fdr < 0.05,]) #15 with only interaction
nrow(normintmodel_output[normintmodel_output$regime_fdr > 0.05 & normintmodel_output$inter_fdr > 0.05,]) #24 n.s. for both

normintmodel_output.sign = normintmodel_output[normintmodel_output$regime_fdr < 0.05,]
nrow(normintmodel_output.sign)

############
#PLOTS OF INTERACTION (not normalized)
mean_ins_tab$study = factor(mean_ins_tab$study, c('Carnes', 'Fabian', 'Hoedjes','Remolina'))
sign.TE.names = as.character(normintmodel_output.sign$TEfam)
for (i in 1:length(sign.TE.names)){
  col.sel = c("study","regime",sign.TE.names[i])
  FB.name = ensembl_casey[ensembl_casey$TEfam == sign.TE.names[i],]$flybase_name
  int.tab = mean_ins_tab[,which(names(mean_ins_tab) %in% col.sel)]
  names(int.tab)[3] = "insertions"

  #PLOT STYLES
  THEME_DEF = theme(axis.title.y = element_text(vjust=0.5,size=18),
                    axis.title.x = element_text(vjust=-0.2,size=18),
                    axis.text.x = element_text(size=18, colour = "black"),
                    axis.text.y = element_text(size=18, colour = "black"),
                    axis.ticks = element_line(size = 1, colour = "black"),
                    axis.ticks.length = unit(0.4, "cm"),
                    axis.line = element_line(size = 1),
                    plot.title = element_text(color="black", size=20, face="bold.italic"),
                    legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
                    legend.key = element_blank(),
                    legend.key.size = unit(1, "cm"),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 20, face = "bold"),
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour="black", size = 1),
                    panel.grid.major = element_line(),
                    panel.grid.minor = element_blank())

  PLOT_LAB = labs(x = expression(bold("Regime")),
                  y = expression(bold("Average Insertions")))

  PLOT_TITLE = ggtitle(FB.name)

  interactplot = ggplot(int.tab, aes(x = regime, color = study, fill=study, group = study, y = int.tab[,3])) + theme_bw() +
    THEME_DEF + PLOT_LAB + PLOT_TITLE +
    stat_summary(fun.y = mean, geom = "point") +
    stat_summary(fun.y = mean, geom = "line",lwd=1.2) +
    scale_colour_manual(name = "Study",values = c("steelblue1", "blue","red","orange"),guide = guide_legend(override.aes = list(shape = c(22,23,21,24)))) +
    scale_fill_manual(values = c("steelblue1", "blue","red","orange")) +
    geom_jitter(position=position_jitter(0.1),alpha=1,size=2.8,
                shape=as.numeric(mgsub(c("Remolina","Carnes","Fabian","Hoedjes"),c(24,22,23,21),mean_ins_tab$study)),color="black") +
    labs(fill = "Study")

  if (intmodel_output.sign[i,4] < 0) {
    plot.save.path = paste(output_path,"interaction_plots/negative/",FB.name,".pdf",sep="")
    ggsave(interactplot, file=plot.save.path, height=6, width=8)
  } else {
    plot.save.path = paste(output_path,"interaction_plots/positive/",FB.name,".pdf",sep="")
    ggsave(interactplot, file=plot.save.path, height=6, width=8)
  }
}

############
#PLOTS OF INTERACTION: PROPORTION OF TOTAL TE CONTENT
sign.TE.names = as.character(normintmodel_output.sign$TEfam)
for (i in 1:length(sign.TE.names)){
  col.sel = c("study","regime",sign.TE.names[i])
  FB.name = ensembl_casey[ensembl_casey$TEfam == TE_bg[i],]$flybase_name
  int.tab = mean_ins_tab[,which(names(mean_ins_tab) %in% col.sel)]
  int.tab[,3] = int.tab[,3]/mean_ins_tab$Ins
  names(int.tab)[3] = "insertions"
  
  #PLOT STYLES
  THEME_DEF = theme(axis.title.y = element_text(vjust=0.5,size=18),
                    axis.title.x = element_text(vjust=-0.2,size=18),
                    axis.text.x = element_text(size=18, colour = "black"),
                    axis.text.y = element_text(size=18, colour = "black"),
                    axis.ticks = element_line(size = 1, colour = "black"),
                    axis.ticks.length = unit(0.4, "cm"),
                    axis.line = element_line(size = 1),
                    plot.title = element_text(color="black", size=20, face="bold.italic"),
                    legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
                    legend.key = element_blank(),
                    legend.key.size = unit(1, "cm"),
                    legend.text = element_text(size = 18),
                    legend.title = element_text(size = 20, face = "bold"),
                    panel.background = element_blank(),
                    panel.border = element_rect(fill = NA, colour="black", size = 1),
                    panel.grid.major = element_line(),
                    panel.grid.minor = element_blank())
  
  PLOT_LAB = labs(x = expression(bold("Regime")), 
                  y = expression(bold("Proportion of Total Insertions")))
  
  PLOT_TITLE = ggtitle(FB.name)
  
  interactplot = ggplot(int.tab, aes(x = regime, color = study, fill=study, group = study, y = int.tab[,3])) + theme_bw() + 
    THEME_DEF + PLOT_LAB + PLOT_TITLE +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.y = mean, geom = "line",lwd=1.2) + 
    scale_colour_manual(name = "Study",values = c("steelblue1", "blue","red","orange"),guide = guide_legend(override.aes = list(shape = c(22,23,21,24)))) +  
    scale_fill_manual(values = c("steelblue1", "blue","red","orange")) + 
    geom_jitter(position=position_jitter(0.1),alpha=1,size=2.8,
                shape=as.numeric(mgsub(c("Remolina","Carnes","Fabian","Hoedjes"),c(24,22,23,21),mean_ins_tab$study)),color="black") + 
    labs(fill = "Study")
  
  if (normintmodel_output.sign[i,4] < 0) {
    plot.save.path = paste(output_path,"interact_TotContNorm/negative/",FB.name,".pdf",sep="")
    ggsave(interactplot, file=plot.save.path, height=6, width=8)
  } else {
    plot.save.path = paste(output_path,"interact_TotContNorm/positive/",FB.name,".pdf",sep="")
    ggsave(interactplot, file=plot.save.path, height=6, width=8)
  }
}
