#TE abundance analysis for approach 1 (all consensus positions considered) and approach 3 (consistent differences between regimes using mean of all consensus positions)

#For Table 1, and Table S2, S3, and Figure S1, S3, S4

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/")
dir.create("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/sign_TE_plots/effcut03/")
dir.create("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/")

#Settings
output_path = "/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/sign_TE_plots/effcut03/"
output_path1 = "/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/"

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
tab = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/hoedjes_all_TE_edited.txt") #File on dryad: https://doi.org/10.5061/dryad.s7h44j13r
names(tab) = c("TEfam","Rep","pos","cov","physcov","hq_cov")
head(tab)
table(sort(tapply(tab$cov,tab$TEfam,max)) == 0) #24 without any reads mapping

#Edit table for analysis
Type = gsub("[1-9]","",tab[,2])
unique(Type) #"CE" "CP" "HE" "HP" "LE" "LP"
Breeding = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Early","Postponed","Early","Postponed","Early","Postponed"),Type)
unique(Breeding)
Diet = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Control","Control","High","High","Low","Low"),Type)
unique(Diet)

tab$Type = Type
tab$Breeding = Breeding
tab$Diet = Diet

tab.edit = cbind(as.character(tab$TEfam),
                 tab$Breeding,
                 tab$Diet,
                 tab$Type,
                 as.character(tab$Rep),
                 tab$pos,
                 tab$cov,
                 tab$physcov,
                 tab$hq_cov)

tab.edit = as.data.frame(tab.edit)
ncol(tab.edit) #9

tab.edit[,6] = numfact(tab.edit[,6])
tab.edit[,7] = numfact(tab.edit[,7])
tab.edit[,8] = numfact(tab.edit[,8])
tab.edit[,9] = numfact(tab.edit[,9])

tab.edit$hqphys_cov = tab.edit[,8] + tab.edit[,9]
ncol(tab.edit) #10
names(tab.edit) = c("TEfam","Breeding","Diet","Type","Rep","pos","cov","physcov","hq_cov","hqphys_cov")
head(tab.edit)

#Subset into TE families in Loop and create output stat table
te.fams = as.vector(unique(tab.edit$TEfam))
stat_output = NULL
for (i in 1:length(te.fams)) {
  te.subset = te.fams[i]
  tab.subset = tab.edit[tab.edit$TEfam == te.subset,]
  
  #check how much of the TE length span is covered
  mean_Pos = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$pos, mean) #insertions per pos, all populations pooled
  length(mean_Pos) == max(tab.subset$pos)+1 #+1 because starts with 0
  
  pos.cov = length(mean_Pos[mean_Pos >= minAvCov]) / length(mean_Pos) #proportion of positions with abundance larger minAvCov (0)
  pos.cov
  
  #averages of types (i.e. 6 regimes)
  mean_Type = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Type, mean)
  diff_PostEarlyCont = mean_Type[2] - mean_Type[1]
  diff_PostEarlyHigh = mean_Type[4] - mean_Type[3]
  diff_PostEarlyLow = mean_Type[6] - mean_Type[5]
  
  log2rat_PostEarlyCont = log2(mean_Type[2] / mean_Type[1])
  log2rat_PostEarlyHigh = log2(mean_Type[4] / mean_Type[3])
  log2rat_PostEarlyLow = log2(mean_Type[6] / mean_Type[5])
  
  #average of replicate populations
  mean_Pop = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Rep, mean)
  
  #averages differences between BREEDING regimes
  mean_Breeding = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Breeding, mean)
  diff_PostEarly =  mean_Breeding[2] - mean_Breeding[1]
  diff_PostEarly
  
  log2rat_PostEarly =  log2(mean_Breeding[2] / mean_Breeding[1])
  log2rat_PostEarly
  
  #averages differences between DIET regimes
  mean_Diet = tapply(tab.subset$hqphys_cov, INDEX = tab.subset$Diet, mean)
  diff_HighCont =  mean_Diet[2] - mean_Diet[1]
  diff_HighLow =  mean_Diet[2] - mean_Diet[3]
  diff_ContLow =  mean_Diet[1] - mean_Diet[3]
  
  log2rat_HighCont =  log2(mean_Diet[2] / mean_Diet[1])
  log2rat_HighLow =  log2(mean_Diet[2] / mean_Diet[3])
  log2rat_ContLow =  log2(mean_Diet[1] / mean_Diet[3])
  
  #FULL - stat model
  model.lm = lm(hqphys_cov ~ Diet + Breeding + Diet * Breeding, data=tab.subset)
  anova(model.lm)
  
  pvalFull_Diet = anova(model.lm)[["Pr(>F)"]][1]
  pvalFull_Breeding = anova(model.lm)[["Pr(>F)"]][2]
  pvalFull_inter = anova(model.lm)[["Pr(>F)"]][3]
  
  fvalue_Diet = round(anova(model.lm)[["F value"]][1],2)
  fvalue_Breeding = round(anova(model.lm)[["F value"]][2],2)
  fvalue_inter = round(anova(model.lm)[["F value"]][3],2)
  
  df_Diet = anova(model.lm)[["Df"]][1]
  df_Breeding = anova(model.lm)[["Df"]][2]
  df_inter = anova(model.lm)[["Df"]][3]
  df_resid = anova(model.lm)[["Df"]][4]
  
  #Effect of Breeding within Diet - stat models
  tab.subset.Cont = tab.subset[tab.subset$Type %in% c("CE","CP"),]
  tab.subset.Cont$Breeding = droplevels(tab.subset.Cont$Breeding)
  tab.subset.Cont$Rep = droplevels(tab.subset.Cont$Rep)
  
  tab.subset.High = tab.subset[tab.subset$Type %in% c("HE","HP"),]
  tab.subset.High$Breeding = droplevels(tab.subset.High$Breeding)
  tab.subset.High$Rep = droplevels(tab.subset.High$Rep)
  
  tab.subset.Low = tab.subset[tab.subset$Type %in% c("LE","LP"),]
  tab.subset.Low$Breeding = droplevels(tab.subset.Low$Breeding)
  tab.subset.Low$Rep = droplevels(tab.subset.Low$Rep)
  
  model.lm.cont = lm(hqphys_cov ~ Breeding + Breeding/Rep, data=tab.subset.Cont)
  model.lm.high = lm(hqphys_cov ~ Breeding + Breeding/Rep, data=tab.subset.High)
  model.lm.low = lm(hqphys_cov ~ Breeding + Breeding/Rep, data=tab.subset.Low)
  
  Df_Cont = anova(model.lm.cont)[["Df"]][1]
  Df_High = anova(model.lm.high)[["Df"]][1]
  Df_Low = anova(model.lm.low)[["Df"]][1]
  Df_Cont.nested = anova(model.lm.cont)[["Df"]][2]
  Df_High.nested = anova(model.lm.high)[["Df"]][2]
  Df_Low.nested = anova(model.lm.low)[["Df"]][2]
  Df_Cont.resid = anova(model.lm.cont)[["Df"]][3]
  Df_High.resid = anova(model.lm.high)[["Df"]][3]
  Df_Low.resid = anova(model.lm.low)[["Df"]][3]
  
  fvalue_Cont = anova(model.lm.cont)[["F value"]][1]
  fvalue_High = anova(model.lm.high)[["F value"]][1]
  fvalue_Low = anova(model.lm.low)[["F value"]][1]
  fvalue_Cont.nested = anova(model.lm.cont)[["F value"]][2]
  fvalue_High.nested = anova(model.lm.high)[["F value"]][2]
  fvalue_Low.nested = anova(model.lm.low)[["F value"]][2]
  
  pval_Cont = anova(model.lm.cont)[["Pr(>F)"]][1]
  pval_High = anova(model.lm.high)[["Pr(>F)"]][1]
  pval_Low = anova(model.lm.low)[["Pr(>F)"]][1]
  pval_Cont.nested = anova(model.lm.cont)[["Pr(>F)"]][2]
  pval_High.nested = anova(model.lm.high)[["Pr(>F)"]][2]
  pval_Low.nested = anova(model.lm.low)[["Pr(>F)"]][2]
  
  int.tab = cbind(te.subset, 
                  pos.cov,
                  round(matrix(mean_Type,ncol=length(mean_Type)),6), #CE       CP       HE       HP       LE       LP 
                  round(matrix(mean_Pop,ncol=length(mean_Pop)),24), 
                  round(matrix(mean_Breeding,ncol=2),3), #Early Postponed 
                  round(matrix(mean_Diet,ncol=3),3), #Control     High      Low 
                  round(diff_PostEarly,3),
                  round(log2rat_PostEarly,3), 
                  round(diff_HighCont,3),
                  round(diff_HighLow,3),
                  round(diff_ContLow,3),
                  round(log2rat_HighCont,3),
                  round(log2rat_HighLow,3),
                  round(log2rat_ContLow,3),
                  
                  round(diff_PostEarlyCont,3),
                  round(diff_PostEarlyHigh,3),
                  round(diff_PostEarlyLow,3),
                  round(log2rat_PostEarlyCont,3),
                  round(log2rat_PostEarlyHigh,3),
                  round(log2rat_PostEarlyLow,3),
                  
                  df_Diet,
                  fvalue_Diet,
                  pvalFull_Diet,
                  df_Breeding,
                  fvalue_Breeding,
                  pvalFull_Breeding,
                  df_inter,
                  fvalue_inter,
                  pvalFull_inter,
                  df_resid,
                  
                  Df_Low,
                  fvalue_Low,
                  pval_Low,
                  Df_Low.nested,
                  fvalue_Low.nested,
                  pval_Low.nested,
                  Df_Low.resid,
                  
                  Df_Cont,
                  fvalue_Cont,
                  pval_Cont,
                  Df_Cont.nested,
                  fvalue_Cont.nested,
                  pval_Cont.nested,
                  Df_Cont.resid,
                  
                  Df_High,
                  fvalue_High,
                  pval_High,
                  Df_High.nested,
                  fvalue_High.nested,
                  pval_High.nested,
                  Df_High.resid)
  
  stat_output = rbind(stat_output,int.tab)
}
stat_output = as.data.frame(stat_output)
names(stat_output) = c("TEfam",
                       "TEpropCov",
                       'CE', 'CP', 'HE', 'HP', 'LE', 'LP',
                       names(mean_Pop),
                       "Mean_Early","Mean_Postponed",
                       "Mean_Control","Mean_High","Mean_Low",
                       "Diff_PostEarly", "log2Rat_PostEarly",
                       "Diff_HighCont", "Diff_HighLow", "Diff_ContLow",
                       "log2Rat_HighCont", "log2Rat_HighLow", "log2Rat_ContLow",
                       "Diff_PostEarlyCont","Diff_PostEarlyHigh","Diff_PostEarlyLow",
                       "log2Rat_PostEarlyCont","log2Rat_PostEarlyHigh","log2Rat_PostEarlyLow",
                       
                       "Df_Diet","F_value_Diet","PvalFull_Diet",
                       "Df_Breeding","F_value_Breeding","PvalFull_Breeding",
                       "Df_Inter","F_value_Inter","PvalFull_Inter",
                       "Df_resid",
                       
                       "Df_Low","F_value_Low","Pval_Low",
                       "Df_Low_nested","F_value_Low_nested","Pval_Low_nested",
                       "Df_Low_resid",
                       
                       "Df_Cont", "F_value_Cont","Pval_Cont",
                       "Df_Cont_nested","F_value_Cont_nested", "Pval_Cont_nested",
                       "Df_Cont_resid",
                       
                       "Df_High","F_value_High","Pval_High",
                       "Df_High_nested","F_value_High_nested","Pval_High_nested",
                       "Df_High_resid")

nrow(stat_output) == length(te.fams)

#factor columns to numeric columns
for (i in 2:ncol(stat_output)) {
  stat_output[,i] = numfact(stat_output[,i])
}

head(stat_output)

#Filter table for TEs that have sufficient coverage
stat_output.filter = stat_output[stat_output$TEpropCov >= propCovCut,]
nrow(stat_output.filter) #115

#Plot filter
pdf(file = paste(output_path1,"FigS3_hoedjes_coverage_filter_minAvCov",minAvCov,"_propCovCut",propCovCut,".pdf", sep =""),width = 8,height = 6)
par(font = 2,font.axis =2,font.lab=2,cex.axis=1.4,cex.lab = 1.4, cex.main = 1.5)
plot(sort(stat_output$TEpropCov), xlab = "TE families sorted by proportion", 
     ylab=paste("Prop. of Positions with NormCov >=",minAvCov,sep=""), 
     main = "Hoedjes2019", xaxt="n", cex = 1.5) 
axis(1,seq(0,180,20))
abline(h=propCovCut,col="red")
text(130,propCovCut-0.04, paste(nrow(stat_output) - nrow(stat_output.filter), " TE families removed"),col = "red",cex=1.4)
text(130,propCovCut+0.04, paste(nrow(stat_output.filter), "TE families retained"),col = "red",cex=1.4)
dev.off()

#Add Bonferroni column and Filter for Bonferroni cut-off
BonfCut = alphaBonf / nrow(stat_output.filter)
stat_output.filter$Bonf_Full = stat_output.filter$PvalFull_Breed < BonfCut
stat_output.filter$Bonf_Low = stat_output.filter$Pval_Low < BonfCut
stat_output.filter$Bonf_Cont = stat_output.filter$Pval_Cont < BonfCut
stat_output.filter$Bonf_High = stat_output.filter$Pval_High < BonfCut
stat_output.filter = merge(stat_output.filter,ensembl_casey,by="TEfam")

stat_output.filter.sign = stat_output.filter[stat_output.filter$Bonf_Full == "TRUE" & abs(stat_output.filter$Diff_PostEarly) > effcut,]
nrow(stat_output.filter.sign) #52

#Full model
table(stat_output.filter[stat_output.filter$PvalFull_Breed < BonfCut & abs(stat_output.filter$Diff_PostEarly) > effcut,]$Mean_Postponed > 
        stat_output.filter[stat_output.filter$PvalFull_Breed < BonfCut & abs(stat_output.filter$Diff_PostEarly) > effcut,]$Mean_Early) #41 more in selected, 11 less
nrow(stat_output.filter[stat_output.filter$PvalFull_Diet < BonfCut &
                          abs(stat_output.filter$Diff_HighCont) > effcut & abs(stat_output.filter$Diff_HighLow) > effcut & abs(stat_output.filter$Diff_ContLow) > effcut,]) #35 for diet

table(stat_output.filter[stat_output.filter$PvalFull_Diet < BonfCut & abs(stat_output.filter$Diff_HighCont) > effcut,]$Diff_HighCont > 0) #44 more in high, 16 less compared to control
table(stat_output.filter[stat_output.filter$PvalFull_Diet < BonfCut & abs(stat_output.filter$Diff_HighLow) > effcut,]$Diff_HighLow > 0) #52 more in high, 8 less compared to low
table(stat_output.filter[stat_output.filter$PvalFull_Diet < BonfCut & abs(stat_output.filter$Diff_ContLow) > effcut,]$Diff_ContLow > 0) #36 more in control, 25 less compared to low
nrow(stat_output.filter[stat_output.filter$PvalFull_Inter< BonfCut,]) #113 significant for interaction

#Within-Diet models for early vs late reproduction (separated by Diet regimes)
table(stat_output.filter[stat_output.filter$Pval_Low < BonfCut & abs(stat_output.filter$Diff_PostEarlyLow) > effcut,]$Diff_PostEarlyLow > 0) #In Low: Post (80) more than early (6)
table(stat_output.filter[stat_output.filter$Pval_Cont < BonfCut & abs(stat_output.filter$Diff_PostEarlyCont) > effcut,]$Diff_PostEarlyCont > 0) #In control: Post (57) more than early (20)
table(stat_output.filter[stat_output.filter$Pval_High < BonfCut & abs(stat_output.filter$Diff_PostEarlyHigh) > effcut,]$Diff_PostEarlyHigh > 0) #In High: Post (3) LESS than early (64)

round(table(stat_output.filter[stat_output.filter$Pval_Low < BonfCut & abs(stat_output.filter$Diff_PostEarlyLow) > effcut,]$Diff_PostEarlyLow > 0) / nrow(stat_output.filter),2) #
round(table(stat_output.filter[stat_output.filter$Pval_Cont < BonfCut & abs(stat_output.filter$Diff_PostEarlyCont) > effcut,]$Diff_PostEarlyCont > 0) / nrow(stat_output.filter),2) #
round(table(stat_output.filter[stat_output.filter$Pval_High < BonfCut & abs(stat_output.filter$Diff_PostEarlyHigh) > effcut,]$Diff_PostEarlyHigh > 0) / nrow(stat_output.filter),2) #

#Save tables
#write.table(stat_output.filter.sign,paste(output_path1,"hoedjesTE_stat_covfilter_sign_Full.txt",sep=""),row.names = F,sep = "\t",quote = F)
write.table(stat_output.filter,paste(output_path1,"hoedjes_TE_stat_covfilter.txt",sep=""),row.names = F,sep = "\t",quote = F)

#Check results with/without filter (only Full model)
stat_output.noFilt = stat_output[!stat_output$PvalFull_Breeding == "NaN",]
stat_output.noFilt.sign = stat_output.noFilt[stat_output.noFilt$PvalFull_Breeding < alphaBonf/nrow(stat_output.noFilt) & abs(stat_output.noFilt$Diff_PostEarly) > effcut,]

#pdf(file = paste(output_path1,"insertion_diff_Hoedjes_Filter_vs_NoFilter_effcut03.pdf", sep =""),width = 7,height = 9)
par(mfrow=c(2,1),font = 2,font.axis =2,font.lab = 2)
plot(sort(stat_output.noFilt$Diff_PostEarly), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Hoedjes2019 (no filters)",ylim = c(-10,10))
abline(h = 0, lty = 2)
diff_count_prop.noFilt.sign = table(stat_output.noFilt.sign$Diff_PostEarly > 0) / nrow(stat_output.noFilt.sign) 
diff_count_prop.noFilt.all = table(stat_output.noFilt.sign$Diff_PostEarly > 0) / nrow(stat_output.noFilt)
diff_count_prop.noFilt = table(stat_output$Diff_PostEarly > 0) / nrow(stat_output)
text(x = 55, y = 9, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.noFilt.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 55, y = 7, paste("Sel>Cont: ",100*round(diff_count_prop.noFilt.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 55, y = 5, paste(nrow(stat_output.noFilt.sign),"/",nrow(stat_output.noFilt)," TEs significant",sep =""),cex=1)

plot(sort(stat_output.filter$Diff_PostEarly), xlab = "TE families sorted by effect size", ylab= expression(bold(paste(delta, "Insertions (S - C)"))), main = "Hoedjes2019 (min. 0.5 NormCov at 80% of TE family positions)",ylim = c(-10,10))
abline(h = 0, lty = 2)
diff_count_prop.sign = table(stat_output.filter.sign$Diff_PostEarly > 0) / nrow(stat_output.filter.sign) 
diff_count_prop.all = table(stat_output.filter.sign$Diff_PostEarly > 0) / nrow(stat_output.filter)
text(x = 45, y = 9, paste("Sel>Cont: ",100*round(diff_count_prop.all,2)[2],"% & Cont>Sel: ",100*round(diff_count_prop.all,2)[1],"% of all tested TE families",sep =""),cex=1)
text(x = 45, y = 7, paste("Sel>Cont: ",100*round(diff_count_prop.sign,2)[2],"% of significant TE families",sep =""),cex=1)
text(x = 45, y = 5, paste(nrow(stat_output.filter.sign),"/",nrow(stat_output.filter)," TEs significant",sep =""),cex=1)
dev.off()

#only plot significant ones from full model
dev.off()
for (i in 1:nrow(stat_output.filter.sign)) {
  tab.subset.sign = tab.edit[tab.edit$TEfam %in% stat_output.filter.sign$TEfam[i],]
  tab.subset.sign.TEfam = as.vector(unique(tab.subset.sign$TEfam))
  tab.subset.sign$Breeding = mgsub(c("Early","Postponed"),c("Cont","Sel"),tab.subset.sign$Breeding)
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
  PLOT_TITLE = ggtitle(paste(flybase.name , "in Hoedjes et al. 2019"))
  
  #COVERAGE ACROSS POSITIONS
  covplot = ggplot(data = tab.subset.sign, aes(x=pos, y=hqphys_cov, group = Breeding, color = Breeding)) + geom_line(alpha=0.2) + theme_bw() #geom_density(alpha = .3)
  covplot2 = covplot + PLOT_TITLE + THEME_DEF + PLOT_LAB + scale_colour_manual(values = c("blue", "red")) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Breeding == "Cont",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('blue')) +
    geom_line(data = tab.subset.sign[tab.subset.sign$Breeding == "Sel",], stat = 'summary', fun.y = mean, linetype = 1, colour = c('red')) +
    guides(colour = guide_legend(override.aes = list(size=3,linetype=1))) + ylim(0, max(tab.subset.sign$hqphys_cov))
  #covplot2
  
  #BOXPLOT FOR REPLICATES
  boxcovplot = ggplot(data = tab.subset.sign, aes(x=Rep, y=hqphys_cov, color = Rep)) + theme_bw() + geom_boxplot()
  boxcovplot2 = boxcovplot + THEME_DEF_BOX + PLOT_LAB_BOX2 + ylim(0, max(tab.subset.sign$hqphys_cov)) + 
    scale_colour_manual(values = rep(c("blue", "blue","blue","blue","red","red","red","red"),4)) + theme(axis.text.x = element_text(size=10, colour = "black"))
  #boxcovplot2
  
  #BOXPLOT FOR REGIMES
  boxcovplotreg = ggplot(data = tab.subset.sign, aes(x=Breeding, y=hqphys_cov, color = Breeding)) + theme_bw() + geom_boxplot() #outlier.shape=NA
  boxcovplotreg2 = boxcovplotreg + THEME_DEF_BOX + PLOT_LAB_BOX + ylim(0, max(tab.subset.sign$hqphys_cov)) + 
    scale_colour_manual(values = c("blue","red")) #+ geom_jitter(width = 0.2,size=0.001)
  #boxcovplotreg2
  
  # Combine plots
  combi.plot = covplot2 / (boxcovplotreg2 | boxcovplot2)
  # combi.plot
  
  # Defining plot names
  combi.plot.path=paste(output_path,flybase.name,"_hoedjes.jpg", sep="")
  
  # Defining save and export
  ggsave(combi.plot, file=combi.plot.path, height=9, width=14)
}

#Plot mean NormCov boxplots to see if differences are driven by outliers
head(stat_output.filter.sign)
nrow(stat_output.filter.sign)
stat_output.filter.Edit = melt(data = stat_output.filter.sign[,c(89,9:32)],value="flybase_name")
stat_output.filter.Edit$Regime = matrix(unlist(strsplit(x = as.character(stat_output.filter.Edit$variable),split = "_")),byrow=T,ncol=2)[,1]
reg1 = mgsub(c(1:4),c(rep("",4)),stat_output.filter.Edit$variable)
stat_output.filter.Edit$Regime = mgsub(unique(reg1),c("Cont","Sel","Cont","Sel","Cont","Sel"),reg1)
stat_output.filter.Edit$RegimeDiet = mgsub(c("CE","CP"),c("ME","MP"),reg1)
stat_output.filter.Edit$RegimeDiet = as.factor(stat_output.filter.Edit$RegimeDiet)
stat_output.filter.Edit$RegimeDiet = factor(stat_output.filter.Edit$RegimeDiet, c("LE","LP","ME","MP","HE","HP"))
names(stat_output.filter.Edit)[3] = "Insertions"
head(stat_output.filter.Edit)

boxplot.mean = ggplot(data = stat_output.filter.Edit , aes(x=RegimeDiet, y=Insertions, color = Regime)) + theme_bw() + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(col = "black",width = 0.1) + 
  scale_colour_manual(values = c("blue","red")) + 
  labs(x = "", y = "Insertions") + theme(axis.title.y = element_text(vjust=0.5,size=18)) + 
  facet_wrap(~flybase_name,scales = "free") 
ggsave(filename = "all_sign_appr1_hoedjes_effcut03.pdf",plot = boxplot.mean,width = 25,height = 20) #on dryad

# stat_output.filter.Edit$variable = mgsub(c("CE","CP"),c("ME","MP"),stat_output.filter.Edit$variable)
# stat_output.filter.Edit$variable = factor(stat_output.filter.Edit$variable, c("LE1", "LE2" ,"LE3" ,"LE4", "LP1" ,"LP2", "LP3", "LP4",
#                                                                               "ME1" ,"ME2" ,"ME3" ,"ME4", "MP1", "MP2", "MP3", "MP4" ,
#                                                                               "HE1", "HE2" ,"HE3" ,"HE4", "HP1", "HP2" ,"HP3" ,"HP4"))
# boxplot.pop = ggplot(data = stat_output.filter.Edit , aes(x=variable, y=Insertions)) + theme_bw() + 
#   geom_point() + 
#   labs(x = "", y = "Insertions") + theme(axis.title.y = element_text(vjust=0.5,size=18),
#                                          axis.text.x = element_text(size=5, colour = "black")) + 
#   facet_wrap(~flybase_name,scales = "free") 
# ggsave(filename = "all_sign_appr1_hoedjes_pop.pdf",plot = boxplot.pop,width = 49,height = 25)

#Filter original input tab for sufficiently abundant TEs
tab.edit.filtered = tab.edit[tab.edit$TEfam %in% stat_output.filter$TEfam,]
write.table(tab.edit.filtered, paste(output_path1,"hoedjes_all_TE_edited_filtered.txt",sep=""),quote=F,sep="\t",col.names = T,row.names = F)

#Statistics on average insertions
hoedjes.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter.txt",header = T)
hoedjes.stat.filt.sign = hoedjes.stat.filt[hoedjes.stat.filt$Bonf_Full == "TRUE" & abs(hoedjes.stat.filt$Diff_PostEarly) > effcut,]

hoedjes.stat.filt.UP = hoedjes.stat.filt[hoedjes.stat.filt$Diff_PostEarly>0,]
hoedjes.stat.filt.DOWN = hoedjes.stat.filt[hoedjes.stat.filt$Diff_PostEarly<0,]

hoedjes.stat.filt.signUP = hoedjes.stat.filt.sign[hoedjes.stat.filt.sign$Diff_PostEarly>0,]
hoedjes.stat.filt.signDOWN = hoedjes.stat.filt.sign[hoedjes.stat.filt.sign$Diff_PostEarly<0,]

nrow(hoedjes.stat.filt) #115
nrow(hoedjes.stat.filt.UP) #78
nrow(hoedjes.stat.filt.DOWN) #37

nrow(hoedjes.stat.filt.sign) #52
nrow(hoedjes.stat.filt.signUP) #41
nrow(hoedjes.stat.filt.signDOWN) #11
  
#Consistent differences: Hoedjes - based on all populations
hoedjes.stat.filt$ConsSC = apply(hoedjes.stat.filt[,c('LE1','LE2','LE3','LE4','CE1','CE2','CE3','CE4','HE1','HE2','HE3','HE4')], 1, FUN=max) < apply(hoedjes.stat.filt[,c('LP1','LP2','LP3','LP4','CP1','CP2','CP3','CP4','HP1','HP2','HP3','HP4')], 1, FUN=min)
hoedjes.stat.filt$ConsCS = apply(hoedjes.stat.filt[,c('LE1','LE2','LE3','LE4','CE1','CE2','CE3','CE4','HE1','HE2','HE3','HE4')], 1, FUN=min) > apply(hoedjes.stat.filt[,c('LP1','LP2','LP3','LP4','CP1','CP2','CP3','CP4','HP1','HP2','HP3','HP4')], 1, FUN=max)
table(hoedjes.stat.filt$ConsSC) #0
table(hoedjes.stat.filt$ConsCS) #0

#Consistent differences: Hoedjes - based on averages of treatments
hoedjes.stat.filt$ConsSCav = apply(hoedjes.stat.filt[,c('LE','CE','HE')], 1, FUN=max) < apply(hoedjes.stat.filt[,c('LP','CP','HP')], 1, FUN=min)
hoedjes.stat.filt$ConsCSav = apply(hoedjes.stat.filt[,c('LE','CE','HE')], 1, FUN=min) > apply(hoedjes.stat.filt[,c('LP','CP','HP')], 1, FUN=max)
table(hoedjes.stat.filt$ConsSCav) #2
table(hoedjes.stat.filt$ConsCSav) #2

#Consistent differences: Hoedjes - within diet
hoedjes.stat.filt$ConsSC_L = apply(hoedjes.stat.filt[,c('LE1','LE2','LE3','LE4')], 1, FUN=max) < apply(hoedjes.stat.filt[,c('LP1','LP2','LP3','LP4')], 1, FUN=min)
hoedjes.stat.filt$ConsCS_L = apply(hoedjes.stat.filt[,c('LE1','LE2','LE3','LE4')], 1, FUN=min) > apply(hoedjes.stat.filt[,c('LP1','LP2','LP3','LP4')], 1, FUN=max)
table(hoedjes.stat.filt$ConsCS_L) #2
table(hoedjes.stat.filt$ConsSC_L) #43

hoedjes.stat.filt$ConsSC_C = apply(hoedjes.stat.filt[,c('CE1','CE2','CE3','CE4')], 1, FUN=max) < apply(hoedjes.stat.filt[,c('CP1','CP2','CP3','CP4')], 1, FUN=min)
hoedjes.stat.filt$ConsCS_C = apply(hoedjes.stat.filt[,c('CE1','CE2','CE3','CE4')], 1, FUN=min) > apply(hoedjes.stat.filt[,c('CP1','CP2','CP3','CP4')], 1, FUN=max)
table(hoedjes.stat.filt$ConsCS_C) #0
table(hoedjes.stat.filt$ConsSC_C) #4

hoedjes.stat.filt$ConsSC_H = apply(hoedjes.stat.filt[,c('HE1','HE2','HE3','HE4')], 1, FUN=max) < apply(hoedjes.stat.filt[,c('HP1','HP2','HP3','HP4')], 1, FUN=min)
hoedjes.stat.filt$ConsCS_H = apply(hoedjes.stat.filt[,c('HE1','HE2','HE3','HE4')], 1, FUN=min) > apply(hoedjes.stat.filt[,c('HP1','HP2','HP3','HP4')], 1, FUN=max)
table(hoedjes.stat.filt$ConsCS_H) #33
table(hoedjes.stat.filt$ConsSC_H) #3

write.table(hoedjes.stat.filt, paste(output_path1,"hoedjes_TE_stat_covfilter_withConsistent.txt",sep=""),quote = F, row.names = F,sep="\t")

#Check correlation between absolute difference in insertions and mean number of insertions
hoedjes.stat.filt$MeanIns = rowMeans(hoedjes.stat.filt[,c("Mean_Early","Mean_Postponed")])
cor.test(abs(hoedjes.stat.filt$Diff_PostEarly), hoedjes.stat.filt$MeanIns)
#t = 12.807, df = 113, p-value < 2.2e-16
#r = 0.7694693

#Check total content
total_cont = t(t(colSums(hoedjes.stat.filt[9:32])))
Pop = row.names(total_cont)
Type = gsub("[1-9]","",Pop)
Breeding = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Early","Postponed","Early","Postponed","Early","Postponed"),Type)
Diet = mgsub(c("CE","CP","HE","HP","LE","LP"),c("Control","Control","High","High","Low","Low"),Type)
total_cont.tab = as.data.frame(cbind(Type,Pop,Breeding,Diet,total_cont))
total_cont.tab[,5] = numfact(total_cont.tab[,5])
names(total_cont.tab)[5] = c("Ins")

#Classes contents
total_cont.tab$LTR = t(t(colSums(hoedjes.stat.filt[hoedjes.stat.filt$TEsubclass == "LTR",][9:32])))
total_cont.tab$TIR = t(t(colSums(hoedjes.stat.filt[hoedjes.stat.filt$TEsubclass == "TIR",][9:32])))
total_cont.tab$nonLTR = t(t(colSums(hoedjes.stat.filt[hoedjes.stat.filt$TEsubclass == "non-LTR",][9:32])))
total_cont.tab$RNA = t(t(colSums(hoedjes.stat.filt[hoedjes.stat.filt$TEclass == "RNA",][9:32])))
total_cont.tab$DNA = t(t(colSums(hoedjes.stat.filt[hoedjes.stat.filt$TEclass == "DNA",][9:32])))
head(total_cont.tab)
write.table(total_cont.tab,paste(output_path1,"hoedjes_TE_total_content.txt",sep=""),quote = F, row.names = F,sep="\t")

#Enrichment of significant TEs with fisher test
table(hoedjes.stat.filt$TEsubclass) / nrow(hoedjes.stat.filt)
table(hoedjes.stat.filt.sign$TEsubclass) / nrow(hoedjes.stat.filt.sign)

#Get rid of single FB before!
hoedjes.stat.filt.noFB = hoedjes.stat.filt[!hoedjes.stat.filt$TEsubclass == "FB",]
hoedjes.stat.filt.noFB$TEsubclass = droplevels(hoedjes.stat.filt.noFB$TEsubclass)
hoedjes.stat.filt.sign = hoedjes.stat.filt.noFB[hoedjes.stat.filt.noFB$Bonf_Full == TRUE & abs(hoedjes.stat.filt.noFB$Diff_PostEarly) > effcut,]
hoedjes.stat.filt.signUP = hoedjes.stat.filt.sign[hoedjes.stat.filt.sign$Diff_PostEarly>0,]
hoedjes.stat.filt.signDOWN = hoedjes.stat.filt.sign[hoedjes.stat.filt.sign$Diff_PostEarly<0,]

#Enrichment for all significant
for(i in 1:length(unique(hoedjes.stat.filt.noFB$TEsubclass)) ){
  type = names(table(hoedjes.stat.filt.sign$TEsubclass))[i]
  m = matrix(c(table(hoedjes.stat.filt.sign$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.sign$TEsubclass)) - table(hoedjes.stat.filt.sign$TEsubclass)[i],
               table(hoedjes.stat.filt.noFB$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.noFB$TEsubclass)) - table(hoedjes.stat.filt.noFB$TEsubclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#ns
table(hoedjes.stat.filt.noFB$TEsubclass) / sum(table(hoedjes.stat.filt.noFB$TEsubclass))
table(hoedjes.stat.filt.sign$TEsubclass) / sum(table(hoedjes.stat.filt.sign$TEsubclass))

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(hoedjes.stat.filt.noFB$TEsubclass)) ){
  type = names(table(hoedjes.stat.filt.signUP$TEsubclass))[i]
  
  m = matrix(c(table(hoedjes.stat.filt.signUP$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.signUP$TEsubclass)) - table(hoedjes.stat.filt.signUP$TEsubclass)[i],
               table(hoedjes.stat.filt.noFB$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.noFB$TEsubclass)) - table(hoedjes.stat.filt.noFB$TEsubclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#only TIRs    twosided: 0.0296161628144587 | greater: 0.01788803865288 | less: 0.9939328550639"
table(hoedjes.stat.filt.signUP$TEsubclass) / sum(table(hoedjes.stat.filt.signUP$TEsubclass))
table(hoedjes.stat.filt.noFB$TEsubclass) / sum(table(hoedjes.stat.filt.noFB$TEsubclass))

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(hoedjes.stat.filt.noFB$TEsubclass)) ){
  type = names(table(hoedjes.stat.filt.signDOWN$TEsubclass))[i]
  
  m = matrix(c(table(hoedjes.stat.filt.signDOWN$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.signDOWN$TEsubclass)) - table(hoedjes.stat.filt.signDOWN$TEsubclass)[i],
               table(hoedjes.stat.filt.noFB$TEsubclass)[i],
               sum(table(hoedjes.stat.filt.noFB$TEsubclass)) - table(hoedjes.stat.filt.noFB$TEsubclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#ns
table(hoedjes.stat.filt.signDOWN$TEsubclass) / sum(table(hoedjes.stat.filt.signDOWN$TEsubclass))
table(hoedjes.stat.filt.noFB$TEsubclass) / sum(table(hoedjes.stat.filt.noFB$TEsubclass))

##### SAME FOR DNA VS RNA #######
#Enrichment for all significant
for(i in 1:length(unique(hoedjes.stat.filt$TEclass)) ){
  type = names(table(hoedjes.stat.filt.sign$TEclass))[i]
  m = matrix(c(table(hoedjes.stat.filt.sign$TEclass)[i],
               sum(table(hoedjes.stat.filt.sign$TEclass)) - table(hoedjes.stat.filt.sign$TEclass)[i],
               table(hoedjes.stat.filt$TEclass)[i],
               sum(table(hoedjes.stat.filt$TEclass)) - table(hoedjes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment

#Enrichment for sign up (i.e. S>C)
for(i in 1:length(unique(hoedjes.stat.filt$TEclass)) ){
  type = names(table(hoedjes.stat.filt.signUP$TEclass))[i]
  
  m = matrix(c(table(hoedjes.stat.filt.signUP$TEclass)[i],
               sum(table(hoedjes.stat.filt.signUP$TEclass)) - table(hoedjes.stat.filt.signUP$TEclass)[i],
               table(hoedjes.stat.filt$TEclass)[i],
               sum(table(hoedjes.stat.filt$TEclass)) - table(hoedjes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2)
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#"DNA    twosided: 0.0323060760740314 | greater: 0.0226928599625851 | less: 0.991962706111103"
table(hoedjes.stat.filt.signUP$TEclass) / sum(table(hoedjes.stat.filt.signUP$TEclass))

#Enrichment for sign down (i.e. C>S)
for(i in 1:length(unique(hoedjes.stat.filt$TEclass)) ){
  type = names(table(hoedjes.stat.filt.signDOWN$TEclass))[i]
  
  m = matrix(c(table(hoedjes.stat.filt.signDOWN$TEclass)[i],
               sum(table(hoedjes.stat.filt.signDOWN$TEclass)) - table(hoedjes.stat.filt.signDOWN$TEclass)[i],
               table(hoedjes.stat.filt$TEclass)[i],
               sum(table(hoedjes.stat.filt$TEclass)) - table(hoedjes.stat.filt$TEclass)[i]),
             ncol=2,nrow=2,dimnames = list(c(),c("Sign","Not sign")))
  
  pval.2 = fisher.test(m,alternative="two.sided")$p.value
  pval.greater = fisher.test(m,alternative="greater")$p.value
  pval.less = fisher.test(m,alternative="less")$p.value
  print(paste(type,"   twosided:",pval.2,"| greater:",pval.greater,"| less:",pval.less))
}
#no enrichment
table(hoedjes.stat.filt.signDOWN$TEclass) / sum(table(hoedjes.stat.filt.signDOWN$TEclass))
#all RNA
