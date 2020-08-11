#TE abundance analysis for approach 1 (all consensus positions considered) and approach 3 (consistent differences between regimes using mean of all consensus positions)

#For Table 1, and Table S2, and Figure S1, S3

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/")
dir.create("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/sign_TE_plots/effcut03/")
dir.create("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/")

#Settings
output_path = "/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/sign_TE_plots/effcut03/" #for many TE plots
output_path1 = "/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/"

#Settings for filtering
minAvCov = 0.5
propCovCut = 0.8 

#Bonferroni alpha level
alphaBonf = 0.01

#Effectsize cutoff
effcut = 0.3

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)

#Load modified table from DeviaTE
tab = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/carnes_all_TE_edited.txt") #File on dryad: https://doi.org/10.5061/dryad.s7h44j13r
names(tab) = c("TEfam","Rep","pos","cov","physcov","hq_cov")
head(tab)
# TEfam     Rep pos     cov physcov hq_cov
# 1  1360 Cont_B1   0 194.602       0 67.510
# 2  1360 Cont_B1   1 208.434       0 76.594
# 3  1360 Cont_B1   2 221.606       0 86.090
# 4  1360 Cont_B1   3 227.510       0 90.343
# 5  1360 Cont_B1   4 230.648       0 92.614
# 6  1360 Cont_B1   5 236.594       0 96.413
table(sort(tapply(tab$cov,tab$TEfam,max)) == 0) #28 without any reads mapping

#Edit table for analysis
sample.split = strsplit(x = as.character(tab[,2]),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
tab.edit = cbind(tab[,1],as.data.frame(sample.split.mat[,1]),tab[,-1])
tab.edit$hqphys_cov = tab.edit[,6] + tab.edit[,7]
names(tab.edit) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab.edit)
#  TEfam Type     Rep pos     cov physcov hq_cov hqphys_cov
# 1  1360 Cont Cont_B1   0 194.602       0 67.510     67.510
# 2  1360 Cont Cont_B1   1 208.434       0 76.594     76.594
# 3  1360 Cont Cont_B1   2 221.606       0 86.090     86.090
# 4  1360 Cont Cont_B1   3 227.510       0 90.343     90.343
# 5  1360 Cont Cont_B1   4 230.648       0 92.614     92.614
# 6  1360 Cont Cont_B1   5 236.594       0 96.413     96.413

#Subset into TE families in Loop and create output stat table
te.fams = as.vector(unique(tab.edit$TEfam))
stat_output = NULL
for (i in 1:length(te.fams)) { 
  te.subset = te.fams[i]
  tab.subset = tab.edit[tab.edit$TEfam == te.subset,]
  
  #check how much of the TE length span is covered
  mean_Pos = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$pos, mean) #insertions per pos, all populations pooled
  length(mean_Pos) == max(tab.subset$pos)+1
  
  pos.cov = length(mean_Pos[mean_Pos >= minAvCov]) / length(mean_Pos) 
  pos.cov
  
  #average difference between types
  mean_ContSel = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Type, mean)
  diff_SelCont =  mean_ContSel[2] - mean_ContSel[1]
  diff_SelCont
  
  log2rat_SelCont =  log2(mean_ContSel[2] / mean_ContSel[1])
  log2rat_SelCont
  
  #average of replicate populations
  mean_ContSelPop = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Rep, mean)
  
  #stat model
  model.lm = lm(hqphys_cov ~ Type+Type/Rep,data=tab.subset)
  fvalue = round(anova(model.lm)[["F value"]][1],2)
  fvalue.nested = round(anova(model.lm)[["F value"]][2],2)
  pvalue = anova(model.lm)[["Pr(>F)"]][1]
  pvalue.nested = anova(model.lm)[["Pr(>F)"]][1]
  df = anova(model.lm)[["Df"]][1]
  df.nested = anova(model.lm)[["Df"]][2]
  df.resid = anova(model.lm)[["Df"]][3]

  #Create output table
  int.tab = cbind(te.subset, pos.cov,
                  round(matrix(mean_ContSelPop,ncol=length(mean_ContSelPop)),3),
                  round(matrix(mean_ContSel,ncol=2),3), 
                  round(diff_SelCont,3),
                  round(log2rat_SelCont,3), 
                  df, fvalue, pvalue,
                  df.nested,fvalue.nested, pvalue.nested,
                  df.resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","TEpropCov",'Cont_B1','Cont_B2','Cont_B3','Cont_B4','Cont_B5','Sel_O1','Sel_O2','Sel_O3','Sel_O4','Sel_O5' ,
                       "Mean_Cont","Mean_Sel","Diff_SelCont","log2RatSelCont","Df","F_value","P_value","Df_nested","F_value_nested","P_value_nested","Df_resid")

nrow(stat_output) == length(te.fams)

#factor columns to numeric columns
for (i in 2:ncol(stat_output)) {
  stat_output[,i] = numfact(stat_output[,i])
}
head(stat_output)

#Filter table for TEs that have sufficient coverage
stat_output.filter = stat_output[stat_output$TEpropCov >= propCovCut,]
nrow(stat_output.filter) #112

#Plot filter
pdf(file = paste(output_path1,"FigS3_carnes_coverage_filter_minAvCov",minAvCov,"_propCovCut",propCovCut,".pdf", sep =""),width = 8,height = 6)
par(font = 2,font.axis =2,font.lab=2,cex.axis=1.4,cex.lab = 1.4, cex.main = 1.5)
plot(sort(stat_output$TEpropCov), xlab = "TE families sorted by proportion", 
     ylab=paste("Prop. of Positions with NormCov >=",minAvCov,sep=""), 
     main = "Carnes2015", xaxt="n", cex = 1.5) 
axis(1,seq(0,180,20))
abline(h=propCovCut,col="red")
text(130,propCovCut-0.04, paste(nrow(stat_output) - nrow(stat_output.filter), " TE families removed"),col = "red",cex=1.4)
text(130,propCovCut+0.04, paste(nrow(stat_output.filter), "TE families retained"),col = "red",cex=1.4)
dev.off()

#Add Bonferroni column and Filter for Bonferroni cut-off
stat_output.filter$Bonf = stat_output.filter$P_value < alphaBonf/nrow(stat_output.filter)
stat_output.filter = merge(stat_output.filter,ensembl_casey,by="TEfam")
stat_output.filter.sign = stat_output.filter[stat_output.filter$Bonf == "TRUE" & abs(stat_output.filter$Diff_SelCont) > effcut,]
nrow(stat_output.filter.sign) #103

#Save tables
#write.table(stat_output.filter.sign,paste(output_path1,"carnes_TE_stat_covfilter_sign.txt",sep=""),row.names = F,sep = "\t",quote = F)
write.table(stat_output.filter,paste(output_path1,"carnes_TE_stat_covfilter.txt",sep=""),row.names = F,sep = "\t",quote = F)

#Check results without filter 
stat_output.noFilt = stat_output[!stat_output$P_value == "NaN",]
stat_output.noFilt.sign = stat_output.noFilt[stat_output.noFilt$P_value < alphaBonf/nrow(stat_output.noFilt) & abs(stat_output.noFilt$Diff_SelCont) > effcut,]

pdf(file = paste(output_path1,"insertion_diff_Carnes_Filter_vs_NoFilter_EffCut03.pdf", sep =""),width = 7,height = 9)
par(mfrow=c(2,1),font = 2,font.axis =2,font.lab = 2)
plot(sort(stat_output.noFilt$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Carnes2015 (no filters)",ylim = c(-55,25))
abline(h = 0, lty = 2)
diff_count_prop.noFilt.sign = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt.sign) 
diff_count_prop.noFilt.all = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt)
diff_count_prop.noFilt = table(stat_output$Diff_SelCont > 0) / nrow(stat_output)
text(x = 55, y = 24, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.noFilt.all,2)[1],"% of all tested TE families",sep =""),cex=0.9)
text(x = 55, y = 18, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.sign,2)[2],"% of significant TE families",sep =""),cex=0.9)
text(x = 55, y = 12, paste(nrow(stat_output.noFilt.sign),"/",nrow(stat_output.noFilt)," TEs significant",sep =""),cex=0.9)

plot(sort(stat_output.filter$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Carnes2015 (min. 0.5 NormCov at 80% of TE family positions)",ylim = c(-55,25), xlim=c(0,110))
abline(h = 0, lty = 2)
diff_count_prop.sign = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter.sign) 
diff_count_prop.all = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter)
text(x = 40, y = 24, paste("Sel>Cont: ",100*round(diff_count_prop.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.all,2)[1],"% of all tested TE families",sep =""),cex=0.9)
text(x = 40, y = 18, paste("Sel>Cont: ",100*round(diff_count_prop.sign,2)[2],"% of significant TE families",sep =""),cex=0.9)
text(x = 40, y = 12, paste(nrow(stat_output.filter.sign),"/",nrow(stat_output.filter)," TEs significant",sep =""),cex=0.9)
dev.off()

#plot significant ones
dev.off()
for (i in 1:nrow(stat_output.filter.sign)) {
  tab.subset.sign = tab.edit[tab.edit$TEfam %in% stat_output.filter.sign$TEfam[i],]
  tab.subset.sign.TEfam = as.vector(unique(tab.subset.sign$TEfam))
  flybase.name = ensembl_casey[ensembl_casey$TEfam %in% tab.subset.sign.TEfam,]$flybase_name
  
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
  
  THEME_DEF_BOX = theme(axis.title.y = element_text(vjust=0.5,size=18),
                        axis.title.x = element_text(vjust=-0.2,size=18),
                        axis.text.x = element_text(size=18, colour = "black"),
                        axis.text.y = element_text(size=18, colour = "black"),
                        axis.ticks = element_line(size = 1, colour = "black"),
                        axis.ticks.length = unit(0.4, "cm"),
                        axis.line = element_line(size = 1),
                        plot.title = element_text(color="black", size=20, face="bold.italic"),
                        legend.background = element_blank(),
                        legend.key = element_blank(),
                        legend.key.size = unit(0, "cm"),
                        legend.text = element_blank(),
                        legend.title = element_blank(),
                        legend.position="none",
                        panel.background = element_blank(),
                        panel.border = element_rect(fill = NA, colour="black", size = 1),
                        panel.grid.major = element_line(),
                        panel.grid.minor = element_blank())
  
  PLOT_LAB = labs(x = expression(bold("Position (bp)")), 
                  y = expression(bold("NormCov")))
  PLOT_LAB_BOX = labs(x = expression(bold("")), 
                      y = expression(bold("NormCov")))
  PLOT_LAB_BOX2 = labs(x = expression(bold("")), 
                       y = expression(bold("")))
  PLOT_TITLE = ggtitle(paste(flybase.name, "in Carnes et al. (2015)"))
  
  #COVERAGE ACROSS POSITIONS
  covplot = ggplot(data = tab.subset.sign, aes(x=pos, y=hqphys_cov, group = Type, color = Type)) + geom_line(alpha = 0.2) + theme_bw() #geom_density(alpha = .3)
  covplot2 = covplot + PLOT_TITLE + THEME_DEF + PLOT_LAB + scale_colour_manual(values = c("blue", "red") ) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Cont",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('blue')) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Sel",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('red')) +
    guides(colour = guide_legend(override.aes = list(size=3,linetype=1))) + ylim(0, max(tab.subset.sign$hqphys_cov))
  # covplot2
  
  #BOXPLOT FOR REPLICATES
  boxcovplot = ggplot(data = tab.subset.sign, aes(x=Rep, y=hqphys_cov, color = Rep)) + theme_bw() + geom_boxplot()
  boxcovplot2 = boxcovplot + THEME_DEF_BOX + PLOT_LAB_BOX2 + ylim(0, max(tab.subset.sign$hqphys_cov)) + 
    scale_colour_manual(values = c("blue", "blue","blue","blue","blue","red","red","red","red","red")) 
  # boxcovplot2
  
  #BOXPLOT FOR REGIMES
  boxcovplotreg = ggplot(data = tab.subset.sign, aes(x=Type, y=hqphys_cov, color = Type)) + theme_bw() + geom_boxplot()
  boxcovplotreg2 = boxcovplotreg + THEME_DEF_BOX + PLOT_LAB_BOX + ylim(0, max(tab.subset.sign$hqphys_cov)) + 
    scale_colour_manual(values = c("blue","red"))
  # boxcovplotreg2
  
  # Combine plots
  combi.plot = covplot2 / (boxcovplotreg2 | boxcovplot2)
  # combi.plot
  
  # Defining plot names
  combi.plot.path=paste(output_path,flybase.name,"_carnes.pdf", sep="")
  
  # Defining save and export
  ggsave(combi.plot, file=combi.plot.path, height=9, width=14)
}

#Plot mean NormCov boxplots to see if differences are driven by outliers
head(stat_output.filter.sign)
nrow(stat_output.filter.sign) #103
stat_output.filter.Edit = melt(data = stat_output.filter.sign[,c(27,3:12)],value="flybase_name")
stat_output.filter.Edit$Regime = matrix(unlist(strsplit(x = as.character(stat_output.filter.Edit$variable),split = "_")),byrow=T,ncol=2)[,1]
names(stat_output.filter.Edit)[3] = "Insertions"

boxplot.mean = ggplot(data = stat_output.filter.Edit , aes(x=Regime, y=Insertions, color = Regime)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(col = "black",width = 0.1) + 
  scale_colour_manual(values = c("blue","red")) + 
  labs(x = "", y = "Insertions") + theme(axis.title.y = element_text(vjust=0.5,size=18)) + 
  facet_wrap(~flybase_name,scales = "free") 
ggsave(filename = "all_sign_appr1_carnes_effcut03.pdf",plot = boxplot.mean,width = 16,height = 16) #on dryad

#Filter original input tab for sufficiently abundant TEs
tab.edit.filtered = tab.edit[tab.edit$TEfam %in% stat_output.filter$TEfam,]
write.table(tab.edit.filtered, paste(output_path1,"carnes_all_TE_edited_filtered.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)

#Tables with filters and annotation
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter.txt",header = T)
carnes.stat.filt.UP = carnes.stat.filt[carnes.stat.filt$Diff_SelCont>0,]
carnes.stat.filt.DOWN = carnes.stat.filt[carnes.stat.filt$Diff_SelCont<0,]

carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE" & abs(carnes.stat.filt$Diff_SelCont) > effcut,]
carnes.stat.filt.signUP = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont>0,]
carnes.stat.filt.signDOWN = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont<0,]

nrow(carnes.stat.filt) #112
nrow(carnes.stat.filt.UP) #89
nrow(carnes.stat.filt.DOWN) #23

nrow(carnes.stat.filt.sign) #103
nrow(carnes.stat.filt.signUP) #82
nrow(carnes.stat.filt.signDOWN) #21

#Approach #3: Consistent differences
carnes.stat.filt$ConsSC = apply(carnes.stat.filt[,c('Cont_B1','Cont_B2','Cont_B3','Cont_B4','Cont_B5')], 1, FUN=max) < apply(carnes.stat.filt[,c('Sel_O1','Sel_O2','Sel_O3','Sel_O4','Sel_O5')], 1, FUN=min)
carnes.stat.filt$ConsCS = apply(carnes.stat.filt[,c('Cont_B1','Cont_B2','Cont_B3','Cont_B4','Cont_B5')], 1, FUN=min) > apply(carnes.stat.filt[,c('Sel_O1','Sel_O2','Sel_O3','Sel_O4','Sel_O5')], 1, FUN=max)
table(carnes.stat.filt$ConsCS) #2
table(carnes.stat.filt$ConsSC) #48

carnes.consistent.TEs = carnes.stat.filt[carnes.stat.filt$ConsCS == "TRUE" | carnes.stat.filt$ConsSC == "TRUE",]
write.table(carnes.stat.filt, paste(output_path1,"carnes_TE_stat_covfilter_withConsistent.txt",sep=""),quote = F, row.names = F,sep="\t")

#Check correlation between absolute difference in insertions and mean number of insertions
carnes.stat.filt$MeanIns = rowMeans(carnes.stat.filt[,c("Mean_Cont","Mean_Sel")])
cor.test(abs(carnes.stat.filt$Diff_SelCont), carnes.stat.filt$MeanIns)
#t = 13.636, df = 110, p-value < 2.2e-16
#r = 0.7926508

#Calculate total TE content
total_cont = t(t(colSums(carnes.stat.filt[3:12])))
total_cont.tab = as.data.frame(cbind(matrix(unlist(strsplit(row.names(total_cont),"_")),ncol=2,byrow=T)[,1],
                                     row.names(total_cont),
                                     total_cont))
total_cont.tab[,3] = numfact(total_cont.tab[,3])
names(total_cont.tab) = c("Regime","Pop","Ins")
head(total_cont.tab)

#Classes contents
total_cont.tab$LTR = t(t(colSums(carnes.stat.filt[carnes.stat.filt$TEsubclass == "LTR",][3:12])))
total_cont.tab$TIR = t(t(colSums(carnes.stat.filt[carnes.stat.filt$TEsubclass == "TIR",][3:12])))
total_cont.tab$nonLTR = t(t(colSums(carnes.stat.filt[carnes.stat.filt$TEsubclass == "non-LTR",][3:12])))
total_cont.tab$RNA = t(t(colSums(carnes.stat.filt[carnes.stat.filt$TEclass == "RNA",][3:12])))
total_cont.tab$DNA = t(t(colSums(carnes.stat.filt[carnes.stat.filt$TEclass == "DNA",][3:12])))
head(total_cont.tab)
write.table(total_cont.tab,paste(output_path1,"carnes_TE_total_content.txt",sep=""),quote = F, row.names = F,sep="\t")

#Enrichment of significant TEs with fisher test
table(carnes.stat.filt$TEsubclass) / nrow(carnes.stat.filt)
table(carnes.stat.filt.sign$TEsubclass) / nrow(carnes.stat.filt.sign)

#Enrichment for all significant
for(i in 1:length(unique(carnes.stat.filt$TEsubclass)) ){
  type = names(table(carnes.stat.filt.sign$TEsubclass))[i]
  m = matrix(c(table(carnes.stat.filt.sign$TEsubclass)[i],
               sum(table(carnes.stat.filt.sign$TEsubclass)) - table(carnes.stat.filt.sign$TEsubclass)[i],
               table(carnes.stat.filt$TEsubclass)[i],
               sum(table(carnes.stat.filt$TEsubclass)) - table(carnes.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(carnes.stat.filt$TEsubclass)) ){
  type = names(table(carnes.stat.filt.signUP$TEsubclass))[i]
  
  m = matrix(c(table(carnes.stat.filt.signUP$TEsubclass)[i],
               sum(table(carnes.stat.filt.signUP$TEsubclass)) - table(carnes.stat.filt.signUP$TEsubclass)[i],
               table(carnes.stat.filt$TEsubclass)[i],
               sum(table(carnes.stat.filt$TEsubclass)) - table(carnes.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(carnes.stat.filt$TEsubclass)) ){
  type = names(table(carnes.stat.filt.signDOWN$TEsubclass))[i]
  
  m = matrix(c(table(carnes.stat.filt.signDOWN$TEsubclass)[i],
               sum(table(carnes.stat.filt.signDOWN$TEsubclass)) - table(carnes.stat.filt.signDOWN$TEsubclass)[i],
               table(carnes.stat.filt$TEsubclass)[i],
               sum(table(carnes.stat.filt$TEsubclass)) - table(carnes.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}

#TIR significant
# [1] "TIR    twosided: 0.0440487189151329 | greater: 1 | less: 0.0193887622991499"

#C>S have TIR underrepresented
table(carnes.stat.filt.signDOWN$TEsubclass)

##### SAME FOR DNA VS RNA #######
#Enrichment for all significant
for(i in 1:length(unique(carnes.stat.filt$TEclass)) ){
  type = names(table(carnes.stat.filt.sign$TEclass))[i]
  m = matrix(c(table(carnes.stat.filt.sign$TEclass)[i],
               sum(table(carnes.stat.filt.sign$TEclass)) - table(carnes.stat.filt.sign$TEclass)[i],
               table(carnes.stat.filt$TEclass)[i],
               sum(table(carnes.stat.filt$TEclass)) - table(carnes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(carnes.stat.filt$TEclass)) ){
  type = names(table(carnes.stat.filt.signUP$TEclass))[i]
  
  m = matrix(c(table(carnes.stat.filt.signUP$TEclass)[i],
               sum(table(carnes.stat.filt.signUP$TEclass)) - table(carnes.stat.filt.signUP$TEclass)[i],
               table(carnes.stat.filt$TEclass)[i],
               sum(table(carnes.stat.filt$TEclass)) - table(carnes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(carnes.stat.filt$TEclass)) ){
  type = names(table(carnes.stat.filt.signDOWN$TEclass))[i]
  
  m = matrix(c(table(carnes.stat.filt.signDOWN$TEclass)[i],
               sum(table(carnes.stat.filt.signDOWN$TEclass)) - table(carnes.stat.filt.signDOWN$TEclass)[i],
               table(carnes.stat.filt$TEclass)[i],
               sum(table(carnes.stat.filt$TEclass)) - table(carnes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#[1] "DNA    twosided: 0.0239192967512212 | greater: 1 | less: 0.0127730021903184"
#[1] "RNA    twosided: 0.0239192967512212 | greater: 0.0127730021903184 | less: 1"

table(carnes.stat.filt.signDOWN$TEclass) 
#C>S have more than expected RNA TEs, i.e. udnberrepresented DNA TEs
