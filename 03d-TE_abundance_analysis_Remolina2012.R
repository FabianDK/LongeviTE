#TE abundance analysis for approach 1 (all consensus positions considered) and approach 3 (consistent differences between regimes using mean of all consensus positions)

#For Table 1, and Table S2, and Figure S1, S3

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/")
dir.create("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/sign_TE_plots/")
dir.create("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/")

#Settings
output_path = "/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/sign_TE_plots/" 
output_path1 = "/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/"

#Settings for filtering
minAvCov = 0.5
propCovCut = 0.8 

#Bonferroni alpha level
alphaBonf = 0.01

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)

#Load modified table from DeviaTE
tab = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/remolina_all_TE_edited.txt") #File on dryad: https://doi.org/10.5061/dryad.s7h44j13r
names(tab) = c("TEfam","Rep","pos","cov","physcov","hq_cov")
head(tab)
# TEfam    Rep pos     cov physcov hq_cov
# 1  1360 Cont_1   0 139.214       0 69.429
# 2  1360 Cont_1   1 146.950       0 75.424
# 3  1360 Cont_1   2 155.339       0 82.290
# 4  1360 Cont_1   3 160.384       0 86.108
# 5  1360 Cont_1   4 161.868       0 87.157
# 6  1360 Cont_1   5 165.964       0 90.006

#Edit table for analysis
sample.split = strsplit(x = as.character(tab[,2]),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
tab.edit = cbind(tab[,1],as.data.frame(sample.split.mat[,1]),tab[,-1])
tab.edit$hqphys_cov = tab.edit[,6] + tab.edit[,7]
names(tab.edit) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab.edit)
# TEfam Type    Rep pos     cov physcov hq_cov hqphys_cov
# 1  1360 Cont Cont_1   0 139.214       0 69.429     69.429
# 2  1360 Cont Cont_1   1 146.950       0 75.424     75.424
# 3  1360 Cont Cont_1   2 155.339       0 82.290     82.290
# 4  1360 Cont Cont_1   3 160.384       0 86.108     86.108
# 5  1360 Cont Cont_1   4 161.868       0 87.157     87.157
# 6  1360 Cont Cont_1   5 165.964       0 90.006     90.006

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
                  df,fvalue, pvalue,
                  df.nested,fvalue.nested, pvalue.nested,
                  df.resid)
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam","TEpropCov",'Cont_1', 'Cont_2', 'Cont_3', 'Sel_1', 'Sel_2', 'Sel_3',
                       "Mean_Cont","Mean_Sel","Diff_SelCont","log2RatSelCont","Df","F_value","P_value","Df_nested","F_value_nested","P_value_nested","Df_resid")

nrow(stat_output) == length(te.fams)

#factor columns to numeric columns
for (i in 2:ncol(stat_output)) {
  stat_output[,i] = numfact(stat_output[,i])
}
head(stat_output)

#Filter table for TEs that have sufficient coverage
stat_output.filter = stat_output[stat_output$TEpropCov >= propCovCut,]
nrow(stat_output.filter) #110

#Plot filter
pdf(file = paste(output_path1,"FigS3_remolina_coverage_filter.pdf", sep =""),width = 8,height = 6)
par(font = 2,font.axis =2,font.lab=2,cex.axis=1.4,cex.lab = 1.4, cex.main = 1.5)
plot(sort(stat_output$TEpropCov), xlab = "TE families sorted by proportion", 
     ylab=paste("Prop. of Positions with NormCov >=",minAvCov,sep=""), 
     main = "Remolina2012", xaxt="n", cex = 1.5) 
axis(1,seq(0,180,20))
abline(h=propCovCut,col="red")
text(130,propCovCut-0.04, paste(nrow(stat_output) - nrow(stat_output.filter), " TE families removed"),col = "red",cex=1.4)
text(130,propCovCut+0.04, paste(nrow(stat_output.filter), "TE families retained"),col = "red",cex=1.4)
dev.off()

#Add Bonferroni column and annotation, then filter for Bonferroni cut-off
stat_output.filter$Bonf = stat_output.filter$P_value < alphaBonf/nrow(stat_output.filter)
stat_output.filter = merge(stat_output.filter,ensembl_casey,by="TEfam")
stat_output.filter.sign = stat_output.filter[stat_output.filter$Bonf == "TRUE",]
nrow(stat_output.filter.sign) #76 are significant after Bonferroni

#Save tables
write.table(stat_output.filter.sign,paste(output_path1,"remolina_TE_stat_covfilter_sign.txt",sep=""),row.names = F,sep = "\t",quote = F)
write.table(stat_output.filter,paste(output_path1,"remolina_TE_stat_covfilter.txt",sep=""),row.names = F,sep = "\t",quote = F)

#Check results with/without filter 
stat_output.noFilt = stat_output[!stat_output$P_value == "NaN",]
stat_output.noFilt.sign = stat_output.noFilt[stat_output.noFilt$P_value < alphaBonf/nrow(stat_output.noFilt),]

#pdf(file = paste(output_path1,"insertion_diff_Remolina_Filter_vs_NoFilter.pdf", sep =""),width = 7,height = 9)
par(mfrow=c(2,1),font = 2,font.axis =2,font.lab = 2)
plot(sort(stat_output.noFilt$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Remolina2012 (no filters)",ylim = c(-1.5,4))
abline(h = 0, lty = 2)
diff_count_prop.noFilt.sign = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt.sign) 
diff_count_prop.noFilt.all = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt)
diff_count_prop.noFilt = table(stat_output$Diff_SelCont > 0) / nrow(stat_output)
text(x = 55, y = 3.5, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.noFilt.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 55, y = 3, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 55, y = 2.5, paste(nrow(stat_output.noFilt.sign),"/",nrow(stat_output.noFilt)," TEs significant",sep =""),cex=1)

plot(sort(stat_output.filter$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Remolina2012 (min. 0.5 NormCov at 80% of TE family positions)",ylim = c(-1.5,4), xlim=c(0,110))
abline(h = 0, lty = 2)
diff_count_prop.sign = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter.sign) 
diff_count_prop.all = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter)
text(x = 45, y = 3.5, paste("Sel>Cont: ",100*round(diff_count_prop.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 45, y = 3, paste("Sel>Cont: ",100*round(diff_count_prop.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 45, y = 2.5, paste(nrow(stat_output.filter.sign),"/",nrow(stat_output.filter)," TEs significant",sep =""),cex=1)
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
  PLOT_TITLE = ggtitle(paste(flybase.name, "in Remolina et al. (2012)"))

  #COVERAGE ACROSS POSITIONS
  covplot = ggplot(data = tab.subset.sign, aes(x=pos, y=hqphys_cov, group = Type, color = Type)) + geom_line(alpha=0.2) + theme_bw() #geom_density(alpha = .3)
  covplot2 = covplot + PLOT_TITLE + THEME_DEF + PLOT_LAB + scale_colour_manual(values = c("blue", "red")) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Cont",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('blue')) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Sel",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('red')) +
    guides(colour = guide_legend(override.aes = list(size=3,linetype=1))) + ylim(0, max(tab.subset.sign$hqphys_cov))
  # covplot2
  
  #BOXPLOT FOR REPLICATES
  boxcovplot = ggplot(data = tab.subset.sign, aes(x=Rep, y=hqphys_cov, color = Rep)) + theme_bw() + geom_boxplot()
  boxcovplot2 = boxcovplot + THEME_DEF_BOX + PLOT_LAB_BOX2 + ylim(0, max(tab.subset.sign$hqphys_cov)) + 
    scale_colour_manual(values = c("blue", "blue","blue","red","red","red")) 
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
  combi.plot.path=paste(output_path,flybase.name,"_remolina.pdf", sep="")
  
  # Defining save and export
  ggsave(combi.plot, file=combi.plot.path, height=9, width=14)
}

#Plot mean NormCov boxplots to see if differences are driven by outliers
head(stat_output.filter)
names(stat_output.filter)
stat_output.filter.Edit = melt(data = stat_output.filter[,c(23,3:8)],value="flybase_name")
stat_output.filter.Edit$Regime = matrix(unlist(strsplit(x = as.character(stat_output.filter.Edit$variable),split = "_")),byrow=T,ncol=2)[,1]
names(stat_output.filter.Edit)[3] = "Insertions"
head(stat_output.filter.Edit)

boxplot.mean = ggplot(data = stat_output.filter.Edit , aes(x=Regime, y=Insertions, color = Regime)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(col = "black",width = 0.1) + 
  scale_colour_manual(values = c("blue","red")) + 
  labs(x = "", y = "Insertions") + theme(axis.title.y = element_text(vjust=0.5,size=18)) + 
  facet_wrap(~flybase_name,scales = "free") 
ggsave(filename = "all_sign_appr1_remolina.pdf",plot = boxplot.mean,width = 16,height = 16) #on dryad

#Filter original input tab for sufficiently abundant TEs
tab.edit.filtered = tab.edit[tab.edit$TEfam %in% stat_output.filter$TEfam,]
write.table(tab.edit.filtered, paste(output_path1,"remolina_all_TE_edited_filtered.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter.txt",header = T)
remo.stat.filt.UP = remo.stat.filt[remo.stat.filt$Diff_SelCont>0,]
remo.stat.filt.DOWN = remo.stat.filt[remo.stat.filt$Diff_SelCont<0,]

remo.stat.filt.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE",]
remo.stat.filt.signUP = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont>0,]
remo.stat.filt.signDOWN = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont<0,]

nrow(remo.stat.filt) #110
nrow(remo.stat.filt.UP) #86
nrow(remo.stat.filt.DOWN) #24

nrow(remo.stat.filt.sign) #76
nrow(remo.stat.filt.signUP) #63
nrow(remo.stat.filt.signDOWN) #13

#Approach #3: Consistent differences
remo.stat.filt$ConsSC = apply(remo.stat.filt[,c('Cont_1','Cont_2' ,'Cont_3')], 1, FUN=max) < apply(remo.stat.filt[,c('Sel_1','Sel_2' ,'Sel_3')], 1, FUN=min)
remo.stat.filt$ConsCS = apply(remo.stat.filt[,c('Cont_1','Cont_2' ,'Cont_3')], 1, FUN=min) > apply(remo.stat.filt[,c('Sel_1','Sel_2' ,'Sel_3')], 1, FUN=max)
table(remo.stat.filt$ConsCS) #0
table(remo.stat.filt$ConsSC) #3

remo.consistent.TEs = remo.stat.filt[remo.stat.filt$ConsCS == "TRUE" | remo.stat.filt$ConsSC == "TRUE",]
write.table(remo.stat.filt, paste(output_path1,"remolina_TE_stat_covfilter_withConsistent.txt",sep=""),quote = F, row.names = F,sep="\t")

#Check correlation between absolute difference in insertions and mean number of insertions
remo.stat.filt$MeanIns = rowMeans(remo.stat.filt[,c('Mean_Sel','Mean_Cont')])
cor.test(abs(remo.stat.filt$Diff_SelCont), remo.stat.filt$MeanIns)
#t = 7.7064, df = 108, p-value = 6.7e-12
#r = 0.5956444

#Calculate total TE content
total_cont = t(t(colSums(remo.stat.filt[3:8])))
total_cont.tab = as.data.frame(cbind(matrix(unlist(strsplit(row.names(total_cont),"_")),ncol=2,byrow=T)[,1],
                    row.names(total_cont),
                    total_cont))
total_cont.tab[,3] = numfact(total_cont.tab[,3])
names(total_cont.tab) = c("Regime","Pop","Ins")
head(total_cont.tab)

#Classes contents
total_cont.tab$LTR = t(t(colSums(remo.stat.filt[remo.stat.filt$TEsubclass == "LTR",][3:8])))
total_cont.tab$TIR = t(t(colSums(remo.stat.filt[remo.stat.filt$TEsubclass == "TIR",][3:8])))
total_cont.tab$nonLTR = t(t(colSums(remo.stat.filt[remo.stat.filt$TEsubclass == "non-LTR",][3:8])))
total_cont.tab$RNA = t(t(colSums(remo.stat.filt[remo.stat.filt$TEclass == "RNA",][3:8])))
total_cont.tab$DNA = t(t(colSums(remo.stat.filt[remo.stat.filt$TEclass == "DNA",][3:8])))
head(total_cont.tab)
write.table(total_cont.tab,paste(output_path1,"remolina_TE_total_content.txt",sep=""),quote = F, row.names = F,sep="\t")

#Enrichment of significant TEs with fisher test
table(remo.stat.filt$TEsubclass) / nrow(remo.stat.filt)
table(remo.stat.filt.sign$TEsubclass) / nrow(remo.stat.filt.sign)

#Enrichment for all significant TEs
for(i in 1:length(unique(remo.stat.filt$TEsubclass)) ){
  type = names(table(remo.stat.filt.sign$TEsubclass))[i]
  m = matrix(c(table(remo.stat.filt.sign$TEsubclass)[i],
               sum(table(remo.stat.filt.sign$TEsubclass)) - table(remo.stat.filt.sign$TEsubclass)[i],
               table(remo.stat.filt$TEsubclass)[i],
               sum(table(remo.stat.filt$TEsubclass)) - table(remo.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(remo.stat.filt$TEsubclass)) ){
  type = names(table(remo.stat.filt.signUP$TEsubclass))[i]
  
  m = matrix(c(table(remo.stat.filt.signUP$TEsubclass)[i],
               sum(table(remo.stat.filt.signUP$TEsubclass)) - table(remo.stat.filt.signUP$TEsubclass)[i],
               table(remo.stat.filt$TEsubclass)[i],
               sum(table(remo.stat.filt$TEsubclass)) - table(remo.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(remo.stat.filt$TEsubclass)) ){
  type = names(table(remo.stat.filt.signDOWN$TEsubclass))[i]
  
  m = matrix(c(table(remo.stat.filt.signDOWN$TEsubclass)[i],
               sum(table(remo.stat.filt.signDOWN$TEsubclass)) - table(remo.stat.filt.signDOWN$TEsubclass)[i],
               table(remo.stat.filt$TEsubclass)[i],
               sum(table(remo.stat.filt$TEsubclass)) - table(remo.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}

#
table(remo.stat.filt.signDOWN$TEsubclass) / sum(table(remo.stat.filt.signDOWN$TEsubclass))
table(remo.stat.filt.signUP$TEsubclass) / sum(table(remo.stat.filt.signUP$TEsubclass))
table(remo.stat.filt.sign$TEsubclass) / sum(table(remo.stat.filt.sign$TEsubclass))
table(remo.stat.filt$TEsubclass) / sum(table(remo.stat.filt$TEsubclass))

##### SAME FOR DNA VS RNA #######
#Enrichment for all significant
for(i in 1:length(unique(remo.stat.filt$TEclass)) ){
  type = names(table(remo.stat.filt.sign$TEclass))[i]
  m = matrix(c(table(remo.stat.filt.sign$TEclass)[i],
               sum(table(remo.stat.filt.sign$TEclass)) - table(remo.stat.filt.sign$TEclass)[i],
               table(remo.stat.filt$TEclass)[i],
               sum(table(remo.stat.filt$TEclass)) - table(remo.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(remo.stat.filt$TEclass)) ){
  type = names(table(remo.stat.filt.signUP$TEclass))[i]
  
  m = matrix(c(table(remo.stat.filt.signUP$TEclass)[i],
               sum(table(remo.stat.filt.signUP$TEclass)) - table(remo.stat.filt.signUP$TEclass)[i],
               table(remo.stat.filt$TEclass)[i],
               sum(table(remo.stat.filt$TEclass)) - table(remo.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(remo.stat.filt$TEclass)) ){
  type = names(table(remo.stat.filt.signDOWN$TEclass))[i]
  
  m = matrix(c(table(remo.stat.filt.signDOWN$TEclass)[i],
               sum(table(remo.stat.filt.signDOWN$TEclass)) - table(remo.stat.filt.signDOWN$TEclass)[i],
               table(remo.stat.filt$TEclass)[i],
               sum(table(remo.stat.filt$TEclass)) - table(remo.stat.filt$TEclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#
table(remo.stat.filt.signDOWN$TEclass) / sum(table(remo.stat.filt.signDOWN$TEclass))
table(remo.stat.filt.signUP$TEclass) / sum(table(remo.stat.filt.signUP$TEclass))
table(remo.stat.filt.sign$TEclass) / sum(table(remo.stat.filt.sign$TEclass))
table(remo.stat.filt$TEclass) / sum(table(remo.stat.filt$TEclass))
