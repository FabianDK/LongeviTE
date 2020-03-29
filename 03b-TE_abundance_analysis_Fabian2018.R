#TE abundance analysis for approach 1 (all consensus positions considered) and approach 3 (consistent differences between regimes using mean of all consensus positions)

#For Table 1, and Table S2, and Figure S1, S3

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/")
dir.create("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/sign_TE_plots/")
dir.create("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/")

#Settings
output_path = "/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/sign_TE_plots/"
output_path1 = "/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/"

#Settings for filtering
minAvCov = 0.5
propCovCut = 0.8 

#Bonferroni alpha level
alphaBonf = 0.01

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)

#Load modified table from DeviaTE
tab = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/fabian_all_TE_edited.txt") #File on dryad: https://doi.org/10.5061/dryad.s7h44j13r
names(tab) = c("TEfam","Rep","pos","cov","physcov","hq_cov")
head(tab)

#Edit table for analysis
sample.split = strsplit(x = as.character(tab[,2]),"_")
sample.split.mat = matrix(unlist(sample.split),ncol=2,byrow=T)
tab.edit = cbind(tab[,1],as.data.frame(sample.split.mat[,1]),tab[,-1])
head(tab.edit)

tab.edit$hqphys_cov = tab.edit[,6] + tab.edit[,7]
names(tab.edit) = c("TEfam","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab.edit)

#Subset into TE families in Loop and create output stat table
te.fams = as.vector(unique(tab.edit$TEfam))
stat_output = NULL
for (i in 1:length(te.fams)) { #length(te.fams)
  te.subset = te.fams[i]
  tab.subset = tab.edit[tab.edit$TEfam == te.subset,]
  
  #check how much of the TE length span is covered
  mean_Pos = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$pos, mean) #insertions per pos, all populations pooled
  length(mean_Pos) == max(tab.subset$pos)+1
  
  pos.cov = length(mean_Pos[mean_Pos >= minAvCov]) / length(mean_Pos) #proportion of positions with abundance larger minAvCov (0)
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

  #Need Bonferroni
  
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
names(stat_output) = c("TEfam","TEpropCov",'Cont_Ra','Cont_Rb','Sel_2La','Sel_2Lb','Sel_La','Sel_Lb',
                       "Mean_Cont","Mean_Sel","Diff_SelCont","log2RatSelCont","Df","F_value","P_value","Df_nested","F_value_nested","P_value_nested","Df_resid")

nrow(stat_output) == length(te.fams)

#factor columns to numeric columns
for (i in 2:ncol(stat_output)) {
  stat_output[,i] = numfact(stat_output[,i])
}

#Filter table for TEs that have sufficient coverage
stat_output.filter = stat_output[stat_output$TEpropCov >= propCovCut,] #Take only TEs with [value propCovCut] proportion of positions larger minAvCov 
nrow(stat_output.filter) #110

#Plot filter
pdf(file = paste(output_path1,"FigS3_fabian_coverage_filter.pdf", sep =""),width = 8,height = 6)
par(font = 2,font.axis =2,font.lab=2,cex.axis=1.4,cex.lab = 1.4, cex.main = 1.5)
plot(sort(stat_output$TEpropCov), xlab = "TE families sorted by proportion", 
     ylab=paste("Prop. of Positions with NormCov >=",minAvCov,sep=""), 
     main = "Fabian2018", xaxt="n", cex = 1.5) 
axis(1,seq(0,180,20))
abline(h=propCovCut,col="red")
text(130,propCovCut-0.04, paste(nrow(stat_output) - nrow(stat_output.filter), " TE families removed"),col = "red",cex=1.4)
text(130,propCovCut+0.04, paste(nrow(stat_output.filter), "TE families retained"),col = "red",cex=1.4)
dev.off()

#Add Bonferroni column and Filter for Bonferroni cut-off
stat_output.filter$Bonf = stat_output.filter$P_value < alphaBonf/nrow(stat_output.filter)
stat_output.filter = merge(stat_output.filter,ensembl_casey,by="TEfam")
stat_output.filter.sign = stat_output.filter[stat_output.filter$Bonf == "TRUE",]
nrow(stat_output.filter.sign) #85

#Save tables
#write.table(stat_output.filter.sign,paste(output_path1,"fabian_TE_stat_covfilter_sign.txt",sep=""),row.names = F,sep = "\t",quote = F)
write.table(stat_output.filter,paste(output_path1,"fabian_TE_stat_covfilter.txt",sep=""),row.names = F,sep = "\t",quote = F)

#Check results with/without filter 
stat_output.noFilt = stat_output[!stat_output$P_value == "NaN",]
stat_output.noFilt.sign = stat_output.noFilt[stat_output.noFilt$P_value < alphaBonf/nrow(stat_output.noFilt),]

#pdf(file = paste(output_path1,"insertion_diff_Fabian_Filter_vs_NoFilter.pdf", sep =""),width = 7,height = 9)
par(mfrow=c(2,1),font = 2,font.axis =2,font.lab = 2)
plot(sort(stat_output.noFilt$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Fabian2018 (no filters)",ylim = c(-10,10))
abline(h = 0, lty = 2)
diff_count_prop.noFilt.sign = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt.sign) 
diff_count_prop.noFilt.all = table(stat_output.noFilt.sign$Diff_SelCont > 0) / nrow(stat_output.noFilt)
diff_count_prop.noFilt = table(stat_output$Diff_SelCont > 0) / nrow(stat_output)
text(x = 55, y = 9, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.noFilt.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 55, y = 7, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 55, y = 5, paste(nrow(stat_output.noFilt.sign),"/",nrow(stat_output.noFilt)," TEs significant",sep =""),cex=1)

plot(sort(stat_output.filter$Diff_SelCont), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Fabian2018 (min. 0.5 NormCov at 80% of TE family positions)",ylim = c(-10,10), xlim=c(0,110))
abline(h = 0, lty = 2)
diff_count_prop.sign = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter.sign) 
diff_count_prop.all = table(stat_output.filter.sign$Diff_SelCont > 0) / nrow(stat_output.filter)
text(x = 43, y = 9, paste("Sel>Cont: ",100*round(diff_count_prop.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 43, y = 7, paste("Sel>Cont: ",100*round(diff_count_prop.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 43, y = 5, paste(nrow(stat_output.filter.sign),"/",nrow(stat_output.filter)," TEs significant",sep =""),cex=1)
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
  PLOT_TITLE = ggtitle(paste(flybase.name, "in Fabian et al. (2018)"))

  #COVERAGE ACROSS POSITIONS
  covplot = ggplot(data = tab.subset.sign, aes(x=pos, y=hqphys_cov, group = Type, color = Type)) + geom_line(alpha = 0.2) + theme_bw() #geom_density(alpha = .3)
  covplot2 = covplot + PLOT_TITLE + THEME_DEF + PLOT_LAB + scale_colour_manual(values = c("blue", "red")) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Cont",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('blue')) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Type == "Sel",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('red')) +
    guides(colour = guide_legend(override.aes = list(size=3,linetype=1))) + ylim(0, max(tab.subset.sign$hqphys_cov))
  # covplot2

  #BOXPLOT FOR REPLICATES
  boxcovplot = ggplot(data = tab.subset.sign, aes(x=Rep, y=hqphys_cov, color = Rep)) + theme_bw() + geom_boxplot()
  boxcovplot2 = boxcovplot + THEME_DEF_BOX + PLOT_LAB_BOX2 + ylim(0, max(tab.subset.sign$hqphys_cov)) +
    scale_colour_manual(values = c("blue", "blue","red","red","red","red"))
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
  combi.plot.path=paste(output_path,flybase.name,"_fabian.pdf", sep="")

  # Defining save and export
  ggsave(combi.plot, file=combi.plot.path, height=9, width=14)
}

#Plot mean NormCov boxplots to see if differences are driven by outliers
head(stat_output.filter)
names(stat_output.filter)
stat_output.filter.Edit = melt(data = stat_output.filter[,c(23,3:8)],value="flybase_name")
stat_output.filter.Edit$Regime = matrix(unlist(strsplit(x = as.character(stat_output.filter.Edit$variable),split = "_")),byrow=T,ncol=2)[,1]
names(stat_output.filter.Edit)[3] = "Insertions"

boxplot.mean = ggplot(data = stat_output.filter.Edit , aes(x=Regime, y=Insertions, color = Regime)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(col = "black",width = 0.1) + 
  scale_colour_manual(values = c("blue","red")) + 
  labs(x = "", y = "Insertions") + theme(axis.title.y = element_text(vjust=0.5,size=18)) + 
  facet_wrap(~flybase_name,scales = "free") 
ggsave(filename = "all_sign_appr1_fabian.pdf",plot = boxplot.mean,width = 16,height = 16) #on dryad

#Filter original input tab for sufficiently abundant TEs
tab.edit.filtered = tab.edit[tab.edit$TEfam %in% stat_output.filter$TEfam,]
write.table(tab.edit.filtered, paste(output_path1,"fabian_all_TE_edited_filtered.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)

#Statistics on average insertions
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter.txt",header = T)
fabian.stat.filt.UP = fabian.stat.filt[fabian.stat.filt$Diff_SelCont>0,]
fabian.stat.filt.DOWN = fabian.stat.filt[fabian.stat.filt$Diff_SelCont<0,]

fabian.stat.filt.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE",]
fabian.stat.filt.signUP = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont>0,]
fabian.stat.filt.signDOWN = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont<0,]

nrow(fabian.stat.filt) #110
nrow(fabian.stat.filt.UP) #64
nrow(fabian.stat.filt.DOWN) #46

nrow(fabian.stat.filt.sign) #85
nrow(fabian.stat.filt.signUP) #51
nrow(fabian.stat.filt.signDOWN) #34

#Approach #3: Consistent differences
fabian.stat.filt$ConsSC = apply(fabian.stat.filt[c('Cont_Ra','Cont_Rb')], 1, FUN=max) < apply(fabian.stat.filt[,c('Sel_2La','Sel_2Lb' ,'Sel_La','Sel_Lb')], 1, FUN=min)
fabian.stat.filt$ConsCS = apply(fabian.stat.filt[c('Cont_Ra','Cont_Rb')], 1, FUN=min) > apply(fabian.stat.filt[,c('Sel_2La','Sel_2Lb' ,'Sel_La','Sel_Lb')], 1, FUN=max)
table(fabian.stat.filt$ConsCS) #8
table(fabian.stat.filt$ConsSC) #15

fabian.consistent.TEs = fabian.stat.filt[fabian.stat.filt$ConsCS == "TRUE" | fabian.stat.filt$ConsSC == "TRUE",]
write.table(fabian.stat.filt, paste(output_path1,"fabian_TE_stat_covfilter_withConsistent.txt",sep=""),quote = F, row.names = F,sep="\t")

#Check correlation between absolute difference in insertions and mean number of insertions
fabian.stat.filt$MeanIns = rowMeans(fabian.stat.filt[,c("Mean_Cont","Mean_Sel")])
cor.test(abs(fabian.stat.filt$Diff_SelCont), fabian.stat.filt$MeanIns)
#t = 6.1484, df = 108, p-value = 1.343e-08
#r = 0.5091912 

#Check total content
total_cont = t(t(colSums(fabian.stat.filt[3:8])))
total_cont.tab = as.data.frame(cbind(matrix(unlist(strsplit(row.names(total_cont),"_")),ncol=2,byrow=T)[,1],
                                     row.names(total_cont),
                                     total_cont))
total_cont.tab[,3] = numfact(total_cont.tab[,3])
names(total_cont.tab) = c("Regime","Pop","Ins")
head(total_cont.tab)

#Classes contents
total_cont.tab$LTR = t(t(colSums(fabian.stat.filt[fabian.stat.filt$TEsubclass == "LTR",][3:8])))
total_cont.tab$TIR = t(t(colSums(fabian.stat.filt[fabian.stat.filt$TEsubclass == "TIR",][3:8])))
total_cont.tab$nonLTR = t(t(colSums(fabian.stat.filt[fabian.stat.filt$TEsubclass == "non-LTR",][3:8])))
total_cont.tab$RNA = t(t(colSums(fabian.stat.filt[fabian.stat.filt$TEclass == "RNA",][3:8])))
total_cont.tab$DNA = t(t(colSums(fabian.stat.filt[fabian.stat.filt$TEclass == "DNA",][3:8])))
head(total_cont.tab)
write.table(total_cont.tab,paste(output_path1,"fabian_TE_total_content.txt",sep=""),quote = F, row.names = F,sep="\t")

#Enrichment of significant TEs with fisher test - meaningful if almost all TEs significant?
table(fabian.stat.filt$TEsubclass) / nrow(fabian.stat.filt)
table(fabian.stat.filt.sign$TEsubclass) / nrow(fabian.stat.filt.sign)

#Enrichment for all significant
for(i in 1:length(unique(fabian.stat.filt$TEsubclass)) ){
  type = names(table(fabian.stat.filt.sign$TEsubclass))[i]
  m = matrix(c(table(fabian.stat.filt.sign$TEsubclass)[i],
               sum(table(fabian.stat.filt.sign$TEsubclass)) - table(fabian.stat.filt.sign$TEsubclass)[i],
               table(fabian.stat.filt$TEsubclass)[i],
               sum(table(fabian.stat.filt$TEsubclass)) - table(fabian.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(fabian.stat.filt$TEsubclass)) ){
  type = names(table(fabian.stat.filt.signUP$TEsubclass))[i]
  
  m = matrix(c(table(fabian.stat.filt.signUP$TEsubclass)[i],
               sum(table(fabian.stat.filt.signUP$TEsubclass)) - table(fabian.stat.filt.signUP$TEsubclass)[i],
               table(fabian.stat.filt$TEsubclass)[i],
               sum(table(fabian.stat.filt$TEsubclass)) - table(fabian.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(fabian.stat.filt$TEsubclass)) ){
  type = names(table(fabian.stat.filt.signDOWN$TEsubclass))[i]
  
  m = matrix(c(table(fabian.stat.filt.signDOWN$TEsubclass)[i],
               sum(table(fabian.stat.filt.signDOWN$TEsubclass)) - table(fabian.stat.filt.signDOWN$TEsubclass)[i],
               table(fabian.stat.filt$TEsubclass)[i],
               sum(table(fabian.stat.filt$TEsubclass)) - table(fabian.stat.filt$TEsubclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

##### SAME FOR DNA VS RNA #######
#Enrichment for all significant
for(i in 1:length(unique(fabian.stat.filt$TEclass)) ){
  type = names(table(fabian.stat.filt.sign$TEclass))[i]
  m = matrix(c(table(fabian.stat.filt.sign$TEclass)[i],
               sum(table(fabian.stat.filt.sign$TEclass)) - table(fabian.stat.filt.sign$TEclass)[i],
               table(fabian.stat.filt$TEclass)[i],
               sum(table(fabian.stat.filt$TEclass)) - table(fabian.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(fabian.stat.filt$TEclass)) ){
  type = names(table(fabian.stat.filt.signUP$TEclass))[i]
  
  m = matrix(c(table(fabian.stat.filt.signUP$TEclass)[i],
               sum(table(fabian.stat.filt.signUP$TEclass)) - table(fabian.stat.filt.signUP$TEclass)[i],
               table(fabian.stat.filt$TEclass)[i],
               sum(table(fabian.stat.filt$TEclass)) - table(fabian.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(fabian.stat.filt$TEclass)) ){
  type = names(table(fabian.stat.filt.signDOWN$TEclass))[i]
  
  m = matrix(c(table(fabian.stat.filt.signDOWN$TEclass)[i],
               sum(table(fabian.stat.filt.signDOWN$TEclass)) - table(fabian.stat.filt.signDOWN$TEclass)[i],
               table(fabian.stat.filt$TEclass)[i],
               sum(table(fabian.stat.filt$TEclass)) - table(fabian.stat.filt$TEclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#
table(fabian.stat.filt.signDOWN$TEclass) / sum(table(fabian.stat.filt.signDOWN$TEclass))
table(fabian.stat.filt.signUP$TEclass) / sum(table(fabian.stat.filt.signUP$TEclass))
table(fabian.stat.filt.sign$TEclass) / sum(table(fabian.stat.filt.sign$TEclass))
table(fabian.stat.filt$TEclass) / sum(table(fabian.stat.filt$TEclass))
