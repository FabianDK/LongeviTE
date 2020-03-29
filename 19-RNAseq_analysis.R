#Differential expression analysis of filtered RNA-seq read counts from Carnes et al. 2015 (PLoS One)
#With comparison to genomic copy number estimates
#For Figure 4; Table S12 to S17, and Figure S11 to S15

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(DESeq2)

setwd("/Users/danfab/Carnes/RNAseq/")
dir.create("/Users/danfab/Carnes/RNAseq/inter")

ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
TE.detect = ensembl_casey[ensembl_casey$species == "Dmel",]
TE.detect$TEfam = droplevels(TE.detect$TEfam)
TE.detect$flybase_name = gsub("Dmel/","",TE.detect[,c(4)],fixed = T)

#TE regulation genes
piRNA_genes = read.table("/Users/danfab/extra_files/piRNA_genes.txt",header=T)
epigen_genes = read.table("/Users/danfab/extra_files/epigenetics_genes.txt",header=T)
transpos_genes = read.table("/Users/danfab/extra_files//transposition_genes_edit.txt")
transpos_genes.private = transpos_genes[!transpos_genes$V2 %in% piRNA_genes$GeneId & !transpos_genes$V2 %in% epigen_genes$GeneId,]
TE.regul.genes = unique(c(as.character(piRNA_genes$FBGN),as.character(epigen_genes$FBGN),as.character(transpos_genes$V1)))

#Abundance data from DeviaTE
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == TRUE,]

#Filtered read counts
tab.mean = read.table("/Users/danfab/Carnes/RNAseq/rawCounts_filtered_Chrom_min400_allDmel_mean.txt",header=T,check.names = F)
ncol(tab.mean) #13394

#Exclude 5 non-Dmel TEs
all.TEs = colnames(tab.mean)[grep(pattern = "FBgn",x = colnames(tab.mean),invert = T)][-c(1:4)]
non.Dmel.TEs = all.TEs[!all.TEs %in% TE.detect$TEfam]
non.Dmel.TEs #"DSRN"      "TARTYAK"   "DNTOMRETA" "DDBARI1"   "DTEII" 
tab.mean = tab.mean[,!colnames(tab.mean) %in% non.Dmel.TEs]
ncol(tab.mean) #13389
#write.table(tab.mean,"rawCounts_filtered_Chrom_min400_onlyDmel_mean.txt",quote=F,sep="\t",row.names=F)

tab.mean.t = as.data.frame(t(tab.mean[,5:ncol(tab.mean)]))
names(tab.mean.t) = paste(tab.mean$Regime,
                          tab.mean$Pop,
                          tab.mean$Age,
                          tab.mean$Sex,sep="_")
ncol(tab.mean.t) #40
nrow(tab.mean.t) #13385
tab.mean.t.round = round(tab.mean.t)

#Check if males have generaly higher count
summed.counts = as.data.frame(cbind(tab.mean[,1:4],colSums(tab.mean.t)))
names(summed.counts)[5] = "CountSum"
head(summed.counts)
countplot = ggplot(data = summed.counts, aes(x=Sex, y=CountSum/1000000, group = Sex, color = Sex)) + geom_boxplot(color="black",outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Sex")), 
       y = expression(bold("Sum of Read Counts in Millions"))) + geom_jitter(width = 0.2)
countplot
t.test(summed.counts[summed.counts$Sex == "M",]$CountSum,
       summed.counts[summed.counts$Sex == "F",]$CountSum) #0.00239
#Females even more, males do not have generally higher count!


#######################################################
#DEseq2 analysis - Main effects for TEs and genes
#######################################################

variables = tab.mean[,1:4]
variables$Age = relevel(variables$Age,ref = "young")
variables$Regime = relevel(variables$Regime,ref = "Cont")
variables$Sex

dds = DESeqDataSetFromMatrix(countData = tab.mean.t.round, 
                             colData = variables, 
                             design = ~ Sex + Age + Regime)

dds1 = DESeq(dds)

res.TE.summary = deseq.results(DEseqobj = dds1, 
                               TEvect = TE.detect$TEfam,
                               cutoff = 0.05,
                               type = "summary")
res.TE.summary
res.TE.list = deseq.results(DEseqobj = dds1, 
                               TEvect = TE.detect$TEfam,
                               cutoff = 0.05,
                               type = "table")
res.TE.list

res.Genes.summary = deseq.results(DEseqobj = dds1, 
                                  TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                  cutoff = 0.05,
                                  type = "summary")
res.Genes.summary

res.Genes.list = deseq.results(DEseqobj = dds1, 
                                  TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                  cutoff = 0.05,
                                  type = "table")


write.table(as.data.frame(res.TE.list),file="TE_expr_all_MAIN_model.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list),file="Genes_all_MAIN_model.txt",quote=F,sep="\t")

#Enrichment of Sex: Genes vs TEs
mat.sex1 = matrix(c(110,122-112,2,
                    7372,13245-12544,5172),ncol=3,byrow = T)
colnames(mat.sex1) = c("M>F","n.s.","F>M")
row.names(mat.sex1) = c("TEs","Genes")
mat.sex1
chisq.test(mat.sex1) #p-value = 3.231e-16
#X-squared = 71.337, df = 2, p-value = 3.231e-16
mat.sex.prop = rbind(mat.sex1[1,] / sum(mat.sex1[1,]),
                      mat.sex1[2,] / sum(mat.sex1[2,]))
row.names(mat.sex.prop) = c("TE families","Genes")

pdf("FigS14_Sex_TEs_vs_Genes_for_Legend1.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.sex.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
legend("topright",legend = row.names(mat.sex.prop), fill= c("brown","grey"),cex=1.5)
#text(x = 2.6, y = 0.98,labels = expression(paste(chi[2],' = 71.34, df = 2, P = 3.2e-16')),cex = 1.6)
dev.off()

pdf("FigS14_Sex_TEs_vs_Genes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.sex.prop,beside = T, ylim = c(0,1), ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
#legend("topright",legend = row.names(mat.sex.prop), fill= c("brown","grey"),cex=1.5)
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.sex.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0),
        add=T)
text(x = 4, y = 0.96,labels = expression(paste(chi^{2},' = 71.34, df = 2, P = 3.2e-16')),cex = 1.8)
dev.off()

##Plot log2FC of sex effect
res.TE.list.sex = na.omit(res.TE.list$Sex_M_vs_F[order(res.TE.list$Sex_M_vs_F$log2FoldChange,decreasing = F),])
res.TE.list.sex$col = mgsub(c("TRUE","FALSE"),c("red","darkgrey"),res.TE.list.sex$padj<0.05)
res.TE.list.sex$TEfam = row.names(res.TE.list.sex)
res.TE.list.sex = merge(res.TE.list.sex,TE.detect[,c(1,4)],by="TEfam")
res.TE.list.sex = res.TE.list.sex[order(res.TE.list.sex$log2FoldChange,decreasing = F),]

pdf("FigS12_Sex_fold_change_allDmel_revised.pdf", width=20, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.8,mar=c(8,4.5,2,2))
plot(1:nrow(res.TE.list.sex),res.TE.list.sex$log2FoldChange,
     pch=21,bg="grey",
     ylim=c(-2,4),xlim=c(3,119),
     axes=F,
     cex=0.01,
     xlab="",ylab="log2 FC (M/F)",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     frame=T)
abline(v = 1:nrow(res.TE.list.sex),col="grey90")
abline(h = 0, lty=2,lwd=2)
arrows(x0 = 1:nrow(res.TE.list.sex),
       y0 = res.TE.list.sex$log2FoldChange,
       x1 = 1:nrow(res.TE.list.sex),
       y1 = c(res.TE.list.sex$log2FoldChange-res.TE.list.sex$lfcSE,
              res.TE.list.sex$log2FoldChange+res.TE.list.sex$lfcSE),
       length = 0.05,angle = 90)
points(1:nrow(res.TE.list.sex),res.TE.list.sex$log2FoldChange,pch=21,bg=res.TE.list.sex$col,cex=2)
axis(2,seq(-2,4,1),seq(-2,4,1),lwd=2,line=0)
axis(1,seq(1,nrow(res.TE.list.sex),1),labels = F,las=2,cex.axis=1,lwd=2)
text(seq(1,nrow(res.TE.list.sex),1), par("usr")[3]-0.15, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = res.TE.list.sex$flybase_name, cex=1)
dev.off()

##Plot log2FC of age effect
res.TE.list.age = na.omit(res.TE.list$Age_old_vs_young[order(res.TE.list$Age_old_vs_young$log2FoldChange,decreasing = F),])
res.TE.list.age$col = mgsub(c("TRUE","FALSE"),c("red","darkgrey"),res.TE.list.age$padj<0.05)
res.TE.list.age$TEfam = row.names(res.TE.list.age)
res.TE.list.age = merge(res.TE.list.age,TE.detect[,c(1,4)],by="TEfam")
res.TE.list.age = res.TE.list.age[order(res.TE.list.age$log2FoldChange,decreasing = F),]

pdf("Age_fold_change_allDmel_revised.pdf", width=20, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.8,mar=c(8,4.5,2,2))
plot(1:nrow(res.TE.list.age),res.TE.list.age$log2FoldChange,
     pch=21,bg="grey",
     ylim=c(-2,4),xlim=c(3,119),
     axes=F,
     cex=0.01,
     xlab="",ylab="log2 FC (Old/Young)",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     frame=T)
abline(v = 1:nrow(res.TE.list.age),col="grey90")
abline(h = 0, lty=2,lwd=2)
arrows(x0 = 1:nrow(res.TE.list.age),
       y0 = res.TE.list.age$log2FoldChange,
       x1 = 1:nrow(res.TE.list.age),
       y1 = c(res.TE.list.age$log2FoldChange-res.TE.list.age$lfcSE,
              res.TE.list.age$log2FoldChange+res.TE.list.age$lfcSE),
       length = 0.05,angle = 90)
points(1:nrow(res.TE.list.age),res.TE.list.age$log2FoldChange,pch=21,bg=res.TE.list.age$col,cex=2)
axis(2,seq(-2,4,1),seq(-2,4,1),lwd=2,line=0)
axis(1,seq(1,nrow(res.TE.list.age),1),labels = F,las=2,cex.axis=1,lwd=2)
text(seq(1,nrow(res.TE.list.age),1), par("usr")[3]-0.15, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = res.TE.list.age$flybase_name, cex=1)
dev.off()

##Plot log2FC of regime effect
res.TE.list.regime = na.omit(res.TE.list$Regime_Sel_vs_Cont[order(res.TE.list$Regime_Sel_vs_Cont$log2FoldChange,decreasing = F),])
res.TE.list.regime$col = mgsub(c("TRUE","FALSE"),c("red","darkgrey"),res.TE.list.regime$padj<0.05)
res.TE.list.regime$TEfam = row.names(res.TE.list.regime)
res.TE.list.regime = merge(res.TE.list.regime,TE.detect[,c(1,4)],by="TEfam")
res.TE.list.regime = res.TE.list.regime[order(res.TE.list.regime$log2FoldChange,decreasing = F),]

pdf("Regime_fold_change_allDmel_revised.pdf", width=20, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.8,mar=c(8,4.5,2,2))
plot(1:nrow(res.TE.list.regime),res.TE.list.regime$log2FoldChange,
     pch=21,bg="grey",
     ylim=c(-2,4),xlim=c(3,119),
     axes=F,
     cex=0.01,
     xlab="",ylab="log2 FC (S/C)",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     frame=T)
abline(v = 1:nrow(res.TE.list.regime),col="grey90")
abline(h = 0, lty=2,lwd=2)
arrows(x0 = 1:nrow(res.TE.list.regime),
       y0 = res.TE.list.regime$log2FoldChange,
       x1 = 1:nrow(res.TE.list.regime),
       y1 = c(res.TE.list.regime$log2FoldChange-res.TE.list.regime$lfcSE,
              res.TE.list.regime$log2FoldChange+res.TE.list.regime$lfcSE),
       length = 0.05,angle = 90)
points(1:nrow(res.TE.list.regime),res.TE.list.regime$log2FoldChange,pch=21,bg=res.TE.list.regime$col,cex=2)
axis(2,seq(-2,4,1),seq(-2,4,1),lwd=2,line=0)
axis(1,seq(1,nrow(res.TE.list.regime),1),labels = F,las=2,cex.axis=1,lwd=2)
text(seq(1,nrow(res.TE.list.regime),1), par("usr")[3]-0.15, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = res.TE.list.regime$flybase_name, cex=1)
dev.off()


#complete model
norm.counts = counts(dds1, normalized=TRUE)
norm.counts.TE = norm.counts[row.names(norm.counts) %in% TE.detect$TEfam,]
head(norm.counts.TE)
nrow(norm.counts.TE) #123
ncol(norm.counts.TE) #40
norm.counts.TE.fin = cbind(variables,t(norm.counts.TE))
norm.counts.TE.fin$Age = factor(norm.counts.TE.fin$Age, 
                                levels = unique(norm.counts.TE.fin$Age[order(norm.counts.TE.fin$Age,
                                                                             decreasing = T)]))
head(norm.counts.TE.fin)
#write.table(norm.counts.TE.fin,"normalized_counts_full_model_allDmel_REVISED_onlyMainFactors.txt",quote=F,sep="\t")


###############################################################
#ANALYZE SEXES SEPARATELY WITH DESEQ2 - INTERACTIONS
###############################################################

tab.mean.F = tab.mean[tab.mean$Sex == "F",]
tab.mean.M = tab.mean[tab.mean$Sex == "M",]
tab.mean.F.t = as.data.frame(t(tab.mean.F[,5:ncol(tab.mean)]))
tab.mean.M.t = as.data.frame(t(tab.mean.M[,5:ncol(tab.mean)]))

names(tab.mean.F.t) = paste(tab.mean.F$Regime,
                            tab.mean.F$Pop,
                            tab.mean.F$Age,
                            tab.mean.F$Sex,sep="_")
names(tab.mean.M.t) = paste(tab.mean.M$Regime,
                            tab.mean.M$Pop,
                            tab.mean.M$Age,
                            tab.mean.M$Sex,sep="_")

ncol(tab.mean.F.t)
ncol(tab.mean.M.t)
tab.mean.F.t.round = round(tab.mean.F.t)
tab.mean.M.t.round = round(tab.mean.M.t)

variables.F = tab.mean.F[,1:3]
variables.M = tab.mean.M[,1:3]
variables.F$Age = relevel(variables.F$Age,ref = "young")
variables.F$Regime = relevel(variables.F$Regime,ref = "Cont")

variables.M$Age = relevel(variables.M$Age,ref = "young")
variables.M$Regime = relevel(variables.M$Regime,ref = "Cont")

dds.M = DESeqDataSetFromMatrix(countData = tab.mean.M.t.round, 
                               colData = variables.M, 
                               design = ~ Age + Regime + Regime*Age)
dds.F = DESeqDataSetFromMatrix(countData = tab.mean.F.t.round, 
                               colData = variables.F, 
                               design = ~ Age + Regime + Regime*Age)

dds1.M = DESeq(dds.M)
dds1.F = DESeq(dds.F)
res.TE.summary.M = deseq.results(DEseqobj = dds1.M, 
                                 TEvect = as.character(TE.detect$TEfam),
                                 cutoff = 0.05, type = "summary")
res.Genes.summary.M = deseq.results(DEseqobj = dds1.M, 
                                    TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "summary")
res.TE.summary.M
res.Genes.summary.M
res.TE.list.M = deseq.results(DEseqobj = dds1.M,   TEvect = TE.detect$TEfam,
                              cutoff = 0.05, type = "table")
res.Genes.list.M = deseq.results(DEseqobj = dds1.M, 
                                    TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "table")

res.TE.summary.F = deseq.results(DEseqobj = dds1.F,  TEvect = TE.detect$TEfam, cutoff = 0.05,type = "summary")
res.Genes.summary.F= deseq.results(DEseqobj = dds1.F,   TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                   cutoff = 0.05, type = "summary")
res.TE.summary.F
res.Genes.summary.F
res.TE.list.F = deseq.results(DEseqobj = dds1.F,   TEvect = TE.detect$TEfam,
                                   cutoff = 0.05, type = "table")
res.Genes.list.F = deseq.results(DEseqobj = dds1.F, 
                                 TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                 cutoff = 0.05, type = "table")

write.table(as.data.frame(res.TE.list.F),file="TE_expr_INTERACTION_model_Females.txt",quote=F,sep="\t")
write.table(as.data.frame(res.TE.list.M),file="TE_expr_INTERACTION_model_Males.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list.F),file="Genes_expr_INTERACTION_model_Females.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list.M),file="Genes_expr_INTERACTION_model_Males.txt",quote=F,sep="\t")

#TE regulation genes - interaction effects
TE.regul.genes
res.TEregulGenesInter.summary.M = deseq.results(DEseqobj = dds1.M,  TEvect = TE.regul.genes, cutoff = 0.05,  type = "summary")
res.TEregulGenesInter.summary.F = deseq.results(DEseqobj = dds1.F, TEvect =  TE.regul.genes, cutoff = 0.05, type = "summary")
res.TEregulGenesInter.summary.M #8
res.TEregulGenesInter.summary.F #11
res.TEregulGenesInter.table.M = deseq.results(DEseqobj = dds1.M,  TEvect = TE.regul.genes, cutoff = 0.05,  type = "table")
res.TEregulGenesInter.table.F = deseq.results(DEseqobj = dds1.F, TEvect =  TE.regul.genes, cutoff = 0.05, type = "table")
as.data.frame(row.names(na.omit(res.TEregulGenesInter.table.M$Ageold.RegimeSel[res.TEregulGenesInter.table.M$Ageold.RegimeSel$padj<0.05,])))
as.data.frame(row.names(na.omit(res.TEregulGenesInter.table.F$Ageold.RegimeSel[res.TEregulGenesInter.table.F$Ageold.RegimeSel$padj<0.05,])))

#Shared genes of interaction
res.Genes.list.M = deseq.results(DEseqobj = dds1.M, 
                                    TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "table")
res.Genes.list.F = deseq.results(DEseqobj = dds1.F, 
                                    TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "table")
shared.Genes.inter = merge( na.omit(res.Genes.list.M$Ageold.RegimeSel[res.Genes.list.M$Ageold.RegimeSel$padj < 0.05,]),
                          na.omit(res.Genes.list.F$Ageold.RegimeSel[res.Genes.list.F$Ageold.RegimeSel$padj < 0.05,]), by ="row.names")
nrow(shared.Genes.inter) #262

###PLOT AGE x REGIME INTERACTIONS USING NORMALIZED COUNTS IN FEMALES (males not significant for interaction)
te.inter.sign.F = na.omit(res.TE.list.F$Ageold.RegimeSel[res.TE.list.F$Ageold.RegimeSel$padj < 0.05,])
norm.counts.F = counts(dds1.F, normalized=TRUE)
norm.counts.F.int = norm.counts.F[row.names(norm.counts.F) %in% row.names(te.inter.sign.F),]
head(norm.counts.F.int)
norm.counts.F.int.TE = cbind(variables.F,t(norm.counts.F.int))
norm.counts.F.int.TE$Age = factor(norm.counts.F.int.TE$Age, 
                                  levels = unique(norm.counts.F.int.TE$Age[order(norm.counts.F.int.TE$Age,
                                                                                 decreasing = T)]))
ncol(norm.counts.F.int.TE ) - 3
norm.counts.F.int.TE.melt = melt(norm.counts.F.int.TE,id.vars = names(norm.counts.F.int.TE)[1:3])

#Figure S13
for (i in 4:ncol(norm.counts.F.int.TE)) {
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
  
  PLOT_LAB = labs(x = expression(bold("")), 
                  y = expression(bold("Normalized Counts")))
  
  Y_LIMITS = c(0,max(norm.counts.F.int.TE[,i],max(norm.counts.F.int.TE[,i]))*1.5)
  norm.counts.F.int.TE$Age = relevel(norm.counts.F.int.TE$Age,ref = "young")
  LEGEND_interactplot.f = ggplot(norm.counts.F.int.TE, aes(x = as.factor(Age), color = Regime, fill=Regime, group = Regime, y = norm.counts.F.int.TE[,i])) + theme_bw() +
    PLOT_LAB + THEME_DEF + ggtitle(names(norm.counts.F.int.TE[i])) +
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.y = mean, geom = "line",lwd=1.2) + 
    scale_colour_manual(name = "Regime",values = c("blue", "red"),guide = guide_legend(override.aes = list(shape = c(21,22)))) + 
    scale_fill_manual(values = c("blue", "red")) + 
    geom_jitter(position=position_jitter(0.1),alpha=1,size=2.8,
                shape=as.numeric(mgsub(c("Cont","Sel"),c(21,22),norm.counts.F.int.TE$Regime)),color="black") + 
    labs(fill = "Regime")
  plot.name = paste("inter/",names(norm.counts.F.int.TE[i]),".pdf",sep="")
  ggsave(LEGEND_interactplot.f, file=plot.name, height=5, width=7)
}

norm.counts.F = counts(dds1.F, normalized=TRUE)
norm.counts.F.TE = norm.counts.F[row.names(norm.counts.F) %in% TE.detect$TEfam,]
head(norm.counts.F.TE)
nrow(norm.counts.F.TE) #123
ncol(norm.counts.F.TE) #20
norm.counts.F.TE.fin = cbind(variables.F,t(norm.counts.F.TE))
norm.counts.F.TE.fin$Age = factor(norm.counts.F.TE.fin$Age, 
                                  levels = unique(norm.counts.F.TE.fin$Age[order(norm.counts.F.TE.fin$Age,
                                                                                 decreasing = T)]))
head(norm.counts.F.TE.fin)
write.table(norm.counts.F.TE.fin,"normalized_counts_female_model_allDmel_withInteraction.txt",quote=F,sep="\t")

#######################################################
#ANALYZE SEXES SEPARATELY WITH DESEQ2 - MAIN EFFECTS
#######################################################

dds.M = DESeqDataSetFromMatrix(countData = tab.mean.M.t.round, 
                               colData = variables.M, 
                               design = ~ Age + Regime)
dds.F = DESeqDataSetFromMatrix(countData = tab.mean.F.t.round, 
                               colData = variables.F, 
                               design = ~ Age + Regime)

dds1.M = DESeq(dds.M)
dds1.F = DESeq(dds.F)

res.TE.summary.M = deseq.results(DEseqobj = dds1.M,  TEvect = as.character(TE.detect$TEfam),
                                 cutoff = 0.05,  type = "summary")
res.Genes.summary.M = deseq.results(DEseqobj = dds1.M, TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "summary")
res.TE.summary.M
res.Genes.summary.M

res.TE.summary.F = deseq.results(DEseqobj = dds1.F,  TEvect = TE.detect$TEfam, cutoff = 0.05,type = "summary")
res.Genes.summary.F= deseq.results(DEseqobj = dds1.F,   TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                   cutoff = 0.05, type = "summary")
res.TE.summary.F
res.Genes.summary.F

res.TE.list.M = deseq.results(DEseqobj = dds1.M,  TEvect = as.character(TE.detect$TEfam), cutoff = 0.05,  type = "table")
res.Genes.list.M = deseq.results(DEseqobj = dds1.M, TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                    cutoff = 0.05, type = "table")
res.TE.list.F = deseq.results(DEseqobj = dds1.F,  TEvect = TE.detect$TEfam, cutoff = 0.05,type = "table")
res.Genes.list.F= deseq.results(DEseqobj = dds1.F,   TEvect = row.names(tab.mean.t.round)[!row.names(tab.mean.t.round) %in% TE.detect$TEfam],
                                   cutoff = 0.05, type = "table")

res.TE.list.M.sign = lapply(res.TE.list.M, function(x) na.omit(x[x$padj<0.05,]))
res.TE.list.F.sign = lapply(res.TE.list.F, function(x) na.omit(x[x$padj<0.05,]))

write.table(as.data.frame(res.TE.list.F),file="TE_expr_MAIN_model_Females.txt",quote=F,sep="\t")
write.table(as.data.frame(res.TE.list.M),file="TE_expr_MAIN_model_Males.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list.M),file="Genes_expr_MAIN_model_Females.txt",quote=F,sep="\t")
write.table(as.data.frame(res.Genes.list.F),file="Genes_expr_MAIN_model_Males.txt",quote=F,sep="\t")

#TE regulation genes - main effects
TE.regul.genes
res.TEregulGenes.summary.M = deseq.results(DEseqobj = dds1.M,  TEvect = TE.regul.genes, cutoff = 0.05,  type = "summary")
res.TEregulGenes.summary.F = deseq.results(DEseqobj = dds1.F, TEvect =  TE.regul.genes, cutoff = 0.05, type = "summary")
res.TEregulGenes.summary.M
res.TEregulGenes.summary.F

res.TEregulGenes.table.M = deseq.results(DEseqobj = dds1.M,  TEvect = TE.regul.genes, cutoff = 0.05,  type = "table")
res.TEregulGenes.table.F = deseq.results(DEseqobj = dds1.F, TEvect =  TE.regul.genes, cutoff = 0.05, type = "table")
res.TEregulGenes.table.M = lapply(res.TEregulGenes.table.M, function(x) cbind(row.names(x),x))
res.TEregulGenes.table.F = lapply(res.TEregulGenes.table.F, function(x) cbind(row.names(x),x))

#Sharing of TEs and Genes between males and females: Regime
shared.TE.reg = merge( na.omit(res.TE.list.M$Regime_Sel_vs_Cont[res.TE.list.M$Regime_Sel_vs_Cont$padj < 0.05,]),
       na.omit(res.TE.list.F$Regime_Sel_vs_Cont[res.TE.list.F$Regime_Sel_vs_Cont$padj < 0.05,]), by ="row.names")
nrow(shared.TE.reg) #19
table(shared.TE.reg$log2FoldChange.x >0 & shared.TE.reg$log2FoldChange.y >0) #3
table(shared.TE.reg$log2FoldChange.x <0 & shared.TE.reg$log2FoldChange.y <0) #16
plot(shared.TE.reg$log2FoldChange.x,shared.TE.reg$log2FoldChange.y)
abline(h=0,v=0)

shared.Genes.reg = merge( na.omit(res.Genes.list.M$Regime_Sel_vs_Cont[res.Genes.list.M$Regime_Sel_vs_Cont$padj < 0.05,]),
                       na.omit(res.Genes.list.F$Regime_Sel_vs_Cont[res.Genes.list.F$Regime_Sel_vs_Cont$padj < 0.05,]), by ="row.names")
nrow(shared.Genes.reg) #1038
table(shared.Genes.reg$log2FoldChange.x >0 & shared.Genes.reg$log2FoldChange.y >0) #469
table(shared.Genes.reg$log2FoldChange.x <0 & shared.Genes.reg$log2FoldChange.y <0) #533
469+533
plot(shared.Genes.reg$log2FoldChange.x,shared.Genes.reg$log2FoldChange.y)
abline(h=0,v=0)

#Table with shared TEs
names(shared.TE.reg)[1] = "TEfam"
shared.TE.reg.annot = merge(shared.TE.reg,TE.detect[,c(1,4,17)],by="TEfam")   
shared.TE.reg.annot[shared.TE.reg.annot$log2FoldChange.x>0,] #TART-A, TART-B- TAHRE the only shared between males and feamles that go up
shared.TE.reg.annot[shared.TE.reg.annot$log2FoldChange.x<0,]
shared.TE.reg.annot[,c("flybase_name", "TEsubclass","log2FoldChange.x","log2FoldChange.y","padj.x","padj.y")]

#Sharing of TEs and Genes between males and females: Age
shared.TE.age = merge( na.omit(res.TE.list.M$Age_old_vs_young[res.TE.list.M$Age_old_vs_young$padj < 0.05,]),
                       na.omit(res.TE.list.F$Age_old_vs_young[res.TE.list.F$Age_old_vs_young$padj < 0.05,]), by ="row.names")
nrow(shared.TE.age) #13
table(shared.TE.age$log2FoldChange.x >0 & shared.TE.age$log2FoldChange.y >0) #9
table(shared.TE.age$log2FoldChange.x <0 & shared.TE.age$log2FoldChange.y <0) #1
plot(shared.TE.age$log2FoldChange.x,shared.TE.age$log2FoldChange.y)
abline(h=0,v=0)

shared.Genes.age = merge( na.omit(res.Genes.list.M$Age_old_vs_young[res.Genes.list.M$Age_old_vs_young$padj < 0.05,]),
                          na.omit(res.Genes.list.F$Age_old_vs_young[res.Genes.list.F$Age_old_vs_young$padj < 0.05,]), by ="row.names")
nrow(shared.Genes.age) #1479
table(shared.Genes.age$log2FoldChange.x >0 & shared.Genes.age$log2FoldChange.y >0) #677
table(shared.Genes.age$log2FoldChange.x <0 & shared.Genes.age$log2FoldChange.y <0) #679
677+679
plot(shared.Genes.age$log2FoldChange.x,shared.Genes.age$log2FoldChange.y)
abline(h=0,v=0)


names(shared.TE.age)[1] = "TEfam"
shared.TE.age.annot = merge(shared.TE.age,TE.detect[,c(1,4,17)],by="TEfam")   
nrow(shared.TE.age.annot[shared.TE.age.annot$log2FoldChange.x>0,]) #12
nrow(shared.TE.age.annot[shared.TE.age.annot$log2FoldChange.x<0,]) #1
shared.TE.age.annot[,c("flybase_name", "TEsubclass","log2FoldChange.x","log2FoldChange.y")]

#Check presence of telomeric TEs
res.TE.list.M.reg.sign = na.omit(res.TE.list.M$Regime_Sel_vs_Cont[res.TE.list.M$Regime_Sel_vs_Cont$padj < 0.05,])
res.TE.list.M.reg.sign$TEfam = row.names(res.TE.list.M.reg.sign)
res.TE.list.M.reg.sign = merge(res.TE.list.M.reg.sign,TE.detect[,c(1,4)],by="TEfam")   
res.TE.list.M.reg.sign = res.TE.list.M.reg.sign[order(res.TE.list.M.reg.sign$log2FoldChange),]
res.TE.list.F.reg.sign = na.omit(res.TE.list.F$Regime_Sel_vs_Cont[res.TE.list.F$Regime_Sel_vs_Cont$padj < 0.05,])
res.TE.list.F.reg.sign$TEfam = row.names(res.TE.list.F.reg.sign)
res.TE.list.F.reg.sign = merge(res.TE.list.F.reg.sign,TE.detect[,c(1,4)],by="TEfam")    
res.TE.list.F.reg.sign = res.TE.list.F.reg.sign[order(res.TE.list.F.reg.sign$log2FoldChange),]

res.TE.list.reg.sign.merge = merge(res.TE.list.M.reg.sign,res.TE.list.F.reg.sign,by="TEfam")
res.TE.list.reg.sign.merge[res.TE.list.reg.sign.merge$log2FoldChange.x>0 & res.TE.list.reg.sign.merge$log2FoldChange.y>0,]
#TART-A, TART-B, TAHRE upregulated in S in both M and F
res.TE.list.reg.sign.merge[res.TE.list.reg.sign.merge$log2FoldChange.x<0 & res.TE.list.reg.sign.merge$log2FoldChange.y<0,]

#Are TEs de-repressed with age: check in males and females separately?
res.TE.summary.M #
res.Genes.summary.M #

res.TE.summary.F #
res.Genes.summary.F #

mat.M.age = matrix(c(107,122-107-1,1,
                     3515,13248-3515-3556,3556),ncol=3,byrow=T)
colnames(mat.M.age) = c("Old>Young","n.s.","Young>Old")
row.names(mat.M.age) = c("TE families","Genes")
mat.M.age
mat.M.age[1,] / sum(mat.M.age[1,])
mat.M.age[2,] / sum(mat.M.age[2,])
chisq.test(mat.M.age) #p-value < 2.2e-16
round(mat.M.age/(rowSums(mat.M.age)),2)

mat.ageM.prop = rbind(mat.M.age[1,] / sum(mat.M.age[1,]),
                      mat.M.age[2,] / sum(mat.M.age[2,]))
row.names(mat.ageM.prop) = c("TE families","Genes")

pdf("FigS14_Age_males_TEs_vs_Genes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.ageM.prop,beside = T, ylim = c(0,1), ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.ageM.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0),
        add=T)
text(x = 4, y = 0.96,labels = expression(paste(chi^{2},' = 230.01, df = 2, P < 2.2e-16')),cex = 1.8)
dev.off()


mat.F.age = matrix(c(10,123-10-5,5,
                     930,11620-930-1091,1091),ncol=3,byrow=T)
colnames(mat.F.age) = c("Old>Young","n.s.","Young>Old")
row.names(mat.F.age) = c("TE families","Genes")
mat.F.age
mat.F.age[1,] / sum(mat.F.age[1,])
mat.F.age[2,] / sum(mat.F.age[2,])
chisq.test(mat.F.age) #p-value = 0.129
round(mat.F.age/(rowSums(mat.F.age)),2)
round(mat.M.age/(rowSums(mat.M.age)),2)

mat.ageF.prop = rbind(mat.F.age[1,] / sum(mat.F.age[1,]),
                      mat.F.age[2,] / sum(mat.F.age[2,]))
row.names(mat.ageF.prop) = c("TE families","Genes")

pdf("FigS14_Age_females_TEs_vs_Genes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.ageF.prop,beside = T, ylim = c(0,1), ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.ageF.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0),
        add=T)
text(x = 3.6, y = 0.96,labels = expression(paste(chi^{2},' = 4.1, df = 2, P = 0.129')),cex = 1.8)
dev.off()


#
mat.M.reg = matrix(c(5,122-5-36,36,
                     1168,13248-1168-1690,1690),ncol=3,byrow=T)
colnames(mat.M.reg) = c("S>C","n.s.","C>S")
row.names(mat.M.reg) = c("TE families","Genes")
mat.M.reg
mat.M.reg[1,] / sum(mat.M.reg[1,])
mat.M.reg[2,] / sum(mat.M.reg[2,])
chisq.test(mat.M.reg) #p-value < 1.389e-07

mat.regM.prop = rbind(mat.M.reg[1,] / sum(mat.M.reg[1,]),
                      mat.M.reg[2,] / sum(mat.M.reg[2,]))
row.names(mat.regM.prop) = c("TE families","Genes")

pdf("FigS14_Regime_males_TEs_vs_Genes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.regM.prop,beside = T, ylim = c(0,1), ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.regM.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0),
        add=T)
text(x = 4, y = 0.96,labels = expression(paste(chi^{2},' = 31.58, df = 2, P = 1.4e-07')),cex = 1.8)
#legend("topright",legend = row.names(mat.ageM.prop), fill= c("brown","grey"),cex=1.5)
#text(x = 6, y = 0.95,labels = expression(chi[2]))
dev.off()


mat.F.reg = matrix(c(4,123-4-23,23,
                     1201,12394-1201-1215,1215),ncol=3,byrow=T)
colnames(mat.F.reg) = c("S>C","n.s.","C>S")
row.names(mat.F.reg) = c("TEs","Genes")
mat.F.reg
mat.F.reg[1,] / sum(mat.F.reg[1,])
mat.F.reg[2,] / sum(mat.F.reg[2,])
chisq.test(mat.F.reg) #p-value = 0.0005313
round(mat.F.reg/(rowSums(mat.F.reg)),2)
round(mat.M.reg/(rowSums(mat.M.reg)),2)

mat.regF.prop = rbind(mat.F.reg[1,] / sum(mat.F.reg[1,]),
                      mat.F.reg[2,] / sum(mat.F.reg[2,]))
row.names(mat.regF.prop) = c("TE families","Genes")

pdf("FigS14_Regime_females_TEs_vs_Genes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=2,mar=c(8,4.5,2,2))
barplot(mat.regF.prop,beside = T, ylim = c(0,1), ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.regF.prop,beside = T, ylim = c(0,1), 
        ylab = "Proportion", main = "", col=c("brown", "grey"),space=c(0.6,0,0.6,0,0.6,0),
        add=T)
#legend("topright",legend = row.names(mat.ageM.prop), fill= c("brown","grey"),cex=1.5)
text(x = 4, y = 0.96,labels = expression(paste(chi^{2},' = 15.08, df = 2, P = 0.0005')),cex = 1.8)
dev.off()

###############################################
#Plot Age and Regime log2 of males vs females, marking significant TEs
#######################################################

row.names(res.TE.list.F$Age_old_vs_young) == row.names(res.TE.list.M$Age_old_vs_young)
nrow(res.TE.list.F$Age_old_vs_young) #123
nrow(res.TE.list.M$Age_old_vs_young) #123

res.TE.list.age.mf = merge(res.TE.list.F$Age_old_vs_young,
                           res.TE.list.M$Age_old_vs_young,by="row.names")
res.TE.list.reg.mf = merge(res.TE.list.F$Regime_Sel_vs_Cont,
                           res.TE.list.M$Regime_Sel_vs_Cont,by="row.names")

#Pearson's correlation
cor.test(res.TE.list.age.mf$log2FoldChange.x,
         res.TE.list.age.mf$log2FoldChange.y,method="pearson") #r = 0.4980267, P = 4.596e-09
cor.test(res.TE.list.reg.mf$log2FoldChange.x,
         res.TE.list.reg.mf$log2FoldChange.y,method="pearson") #r = 0.728792 , p =  2.2e-16

pdf("Fig4_Regime_Age_M_F_change_allDmel1_revised_MainEffectOnly1.pdf", width=16, height=8)
par(mfrow = c(1,2),font=1,font.lab=1,font.axis=1,cex.lab=1.9,cex.axis=1.9,mar=c(8,5,2,1.5))
plot(res.TE.list.reg.mf$log2FoldChange.x,
     res.TE.list.reg.mf$log2FoldChange.y,
     ylab=expression(paste('Log'[2], ' FC (S/C males)')),xlab=expression(paste('Log'[2], ' FC (S/C females)')),
     xlim=c(-3,2),ylim=c(-3,2),
     cex = 1.8,
     pch = 21, bg="darkgrey", col = "darkgrey",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     axes=F,frame=T)
axis(1,seq(-3,2,1),seq(-3,2,1),lwd=2)
axis(2,seq(-3,2,1),seq(-3,2,1),lwd=2)
abline(h=0,lty = 2,lwd = 2)
abline(v=0,lty = 2,lwd = 2)
points(res.TE.list.reg.mf[res.TE.list.reg.mf$padj.x<0.05,]$log2FoldChange.x,
       res.TE.list.reg.mf[res.TE.list.reg.mf$padj.x<0.05,]$log2FoldChange.y, 
       pch=21,bg="red",cex=2.2)
points(res.TE.list.reg.mf[res.TE.list.reg.mf$padj.y<0.05,]$log2FoldChange.x,
       res.TE.list.reg.mf[res.TE.list.reg.mf$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="blue",cex=2.2)
points(res.TE.list.reg.mf[res.TE.list.reg.mf$padj.x<0.05 & res.TE.list.reg.mf$padj.y<0.05,]$log2FoldChange.x,
       res.TE.list.reg.mf[res.TE.list.reg.mf$padj.x<0.05 & res.TE.list.reg.mf$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="orange",cex=2.2)
text(1.5,1.6,labels = "r = 0.73***",col = "black",cex=1.9)
text(-1.8,2,labels = "Padj<0.05 in males only",col = "blue",cex=1.5)
text(-1.8,1.7,labels = "Padj<0.05 in females only",col = "red",cex=1.5)
text(-1.8,1.4,labels = "Padj<0.05 in both sexes",col = "orange",cex=1.5)

#SWITCH SIGN TO HAVE YOUNG/OLD NOT OLD/YOUNG (thisi sb etter for the comparison of age vs regime changes)
plot(-res.TE.list.age.mf$log2FoldChange.x,
     -res.TE.list.age.mf$log2FoldChange.y,
     ylab=expression(paste('Log'[2], ' FC (Young/Old males)')),xlab=expression(paste('Log'[2], ' FC (Young/Old females)')),
     xlim=c(-2.3,2),ylim=c(-2.3,2),
     cex = 1.8,
     pch = 21, bg="darkgrey", col = "darkgrey",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     axes=F,frame=T)
axis(1,seq(-3,2,1),seq(-3,2,1),lwd=2)
axis(2,seq(-3,2,1),seq(-3,2,1),lwd=2)
abline(h=0,lty = 2,lwd = 2)
abline(v=0,lty = 2,lwd = 2)
points(-res.TE.list.age.mf[res.TE.list.age.mf$padj.x<0.05,]$log2FoldChange.x,
       -res.TE.list.age.mf[res.TE.list.age.mf$padj.x<0.05,]$log2FoldChange.y, 
       pch=21,bg="red",cex=2.2)
points(-res.TE.list.age.mf[res.TE.list.age.mf$padj.y<0.05,]$log2FoldChange.x,
       -res.TE.list.age.mf[res.TE.list.age.mf$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="blue",cex=2.2)
points(-res.TE.list.age.mf[res.TE.list.age.mf$padj.x<0.05 & res.TE.list.age.mf$padj.y<0.05,]$log2FoldChange.x,
       -res.TE.list.age.mf[res.TE.list.age.mf$padj.x<0.05 & res.TE.list.age.mf$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="orange",cex=2.2)
text(1.5,1.6,labels = "r = 0.5***",col = "black",cex=1.9)
text(-1.2,2,labels = "Padj<0.05 in males only",col = "blue",cex=1.5)
text(-1.2,1.7,labels = "Padj<0.05 in females only",col = "red",cex=1.5)
text(-1.2,1.4,labels = "Padj<0.05 in both sexes",col = "orange",cex=1.5)
dev.off()

#######################################################
#Barcharts for number of DE TEs
#######################################################

de_te = read.table("/Users/danfab/extra_files/number_prop_DE_TEs_Main_factors.txt",header=T) #Created this table in excel
de_te
de_te$Level = factor(de_te$Level, c("young","old","control","selected","f","m","f1","m1"))
de_te$Level

# Stacked barplot with multiple groups
barplot.DE= ggplot(data=de_te, aes(x=paste(Factor,Sex,sep="_"), y=Prop*c(1,-1), fill=Level)) +
  geom_bar(stat="identity", color = "black") + 
  scale_fill_manual(values=c("goldenrod2", "goldenrod4","blue","red","black","lightblue","green4","white")) + 
  geom_hline(yintercept = 0) +
  THEME_DEF + theme(axis.text.y = element_text(size=28, colour = "black"),
                    axis.text.x = element_text(size=28, colour = "black"),
                    axis.title.y = element_text(vjust=0.5,size=30)) +
  ylim(-1,1) + ylab("Proportion D.E. TEs") + xlab("") +
  scale_x_discrete(limits=c("Sex_NA", "Regime_m", "Regime_f","Age_m","Age_f","RxA_RxA"), 
                   labels=c("Sex_NA" = "Sex", 
                            "Regime_m" = "Regime 
M", "Regime_f" = "Regime
F",
                            "Age_m" = "Age 
M", "Age_f" = "Age 
F","RxA_RxA" = "RxA 
interaction")) + annotate("text", x = 6, y=1, label = "n = 123",size = 9)
barplot.DE    
ggsave(barplot.DE, file="Fig4_DE_TE_Numbers_MainEffectOnly.pdf", height=10, width=12)


########################################
#Age change vs Regime change - TEs
########################################

cor.test(-res.TE.list.F$Age_old_vs_young$log2FoldChange,
         res.TE.list.F$Regime_Sel_vs_Cont$log2FoldChange) #0.2082257 , p = 0.02082
cor.test(-res.TE.list.M$Age_old_vs_youn$log2FoldChange,
         res.TE.list.M$Regime_Sel_vs_Cont$log2FoldChange) #-0.01 , p = 0.8749
#TEs that are upregulated with age (i.e. more in older), are downregulated in selected, in females

#Age change vs Regime change - Genes
cor.test(-res.Genes.list.F$Age_old_vs_young$log2FoldChange,
         res.Genes.list.F$Regime_Sel_vs_Cont$log2FoldChange) #-0.1019146  , p = 2.2e-16
cor.test(-res.Genes.list.M$Age_old_vs_young$log2FoldChange,
         res.Genes.list.M$Regime_Sel_vs_Cont$log2FoldChange) #0.07859942 , p = < 2.2e-16

#boostrap correlation for genes
gene.fbgn = intersect(row.names(res.Genes.list.F$Age_old_vs_young),row.names(res.Genes.list.F$Regime_Sel_vs_Cont))
allcor.f = NULL
allcor.m = NULL
for (i in 1:1000){
  gene.fbgn.random = sample(gene.fbgn,123) #same length as TEs
  int.age.f = res.Genes.list.F$Age_old_vs_young[row.names(res.Genes.list.F$Age_old_vs_young) %in% gene.fbgn.random,]
  int.reg.f = res.Genes.list.F$Regime_Sel_vs_Cont[row.names(res.Genes.list.F$Regime_Sel_vs_Cont) %in% gene.fbgn.random,]
  cor.f = cor.test(-int.age.f$log2FoldChange,int.reg.f$log2FoldChange)$estimate
  allcor.f = c(allcor.f, cor.f)
  
  int.age.m = res.Genes.list.M$Age_old_vs_young[row.names(res.Genes.list.M$Age_old_vs_young) %in% gene.fbgn.random,]
  int.reg.m = res.Genes.list.M$Regime_Sel_vs_Cont[row.names(res.Genes.list.M$Regime_Sel_vs_Cont) %in% gene.fbgn.random,]
  cor.m = cor.test(-int.age.m$log2FoldChange,int.reg.m$log2FoldChange)$estimate
  allcor.m = c(allcor.m, cor.m)
}
hist(allcor.f,breaks = 10)
abline(v = -0.1019146, col="red")
abline(v = mean(allcor.f), col="blue")
hist(allcor.m,breaks = 10)
abline(v = 0.07859942, col="red")
abline(v = mean(allcor.m), col="blue")
library(Rmisc)
CI(allcor.f,ci=0.95)
CI(allcor.m,ci=0.95)

res.TE.list.age.reg.F = merge(res.TE.list.F$Age_old_vs_young,
                              res.TE.list.F$Regime_Sel_vs_Cont,by="row.names")
res.TE.list.age.reg.M = merge(res.TE.list.M$Age_old_vs_young,
                              res.TE.list.M$Regime_Sel_vs_Cont,by="row.names")

########################################
#PLOT - AGE VS REGIME CHANGE
########################################

pdf("Fig4_Regime_vs_Age_change_allDmel_revised_MainEffectOnly1.pdf", width=16, height=8)
par(mfrow = c(1,2),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.9,mar=c(8,5,2,1.5))
plot(-res.TE.list.F$Age_old_vs_young$log2FoldChange,
     res.TE.list.F$Regime_Sel_vs_Cont$log2FoldChange,
     ylab=expression(paste('Log'[2], ' FC (S/C)')),xlab=expression(paste('Log'[2], ' FC (Young/Old)')),
     xlim=c(-3,2),ylim=c(-3,2),
     cex = 1.8,
     pch = 21, bg="darkgrey", col = "darkgrey",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     axes=F,frame=T)
axis(1,seq(-3,2,1),seq(-3,2,1),lwd=2)
axis(2,seq(-3,2,1),seq(-3,2,1),lwd=2)
abline(h=0,lty = 2,lwd = 2)
abline(v=0,lty = 2,lwd = 2)
points(-res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.x<0.05,]$log2FoldChange.x, #sign age
       res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.x<0.05,]$log2FoldChange.y, 
       pch=21,bg="red",cex=2.2)
points(-res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.y<0.05,]$log2FoldChange.x, #sign regime
       res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="blue",cex=2.2)
points(-res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.x<0.05 & res.TE.list.age.reg.F$padj.y<0.05,]$log2FoldChange.x,
       res.TE.list.age.reg.F[res.TE.list.age.reg.F$padj.x<0.05 & res.TE.list.age.reg.F$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="orange",cex=2.2) #sign both

text(1.5,1.6,labels = "r = 0.21*",col = "black",cex=1.9)
text(-1.8,2,labels = "Padj<0.05 for Regime",col = "blue",cex=1.5)
text(-1.8,1.7,labels = "Padj<0.05 for Age",col = "red",cex=1.5)
text(-1.8,1.4,labels = "Padj<0.05 for both",col = "orange",cex=1.5)

plot(-res.TE.list.M$Age_old_vs_young$log2FoldChange,
     res.TE.list.M$Regime_Sel_vs_Cont$log2FoldChange,
     ylab="",xlab=expression(paste('Log'[2], ' FC (Young/Old)')),
     xlim=c(-3,2),ylim=c(-3,2),
     cex = 1.8,
     pch = 21, bg="darkgrey", col = "darkgrey",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE),
     axes=F,frame=T)
axis(1,seq(-3,2,1),seq(-3,2,1),lwd=2)
axis(2,seq(-3,2,1),seq(-3,2,1),lwd=2)
abline(h=0,lty = 2,lwd = 2)
abline(v=0,lty = 2,lwd = 2)
points(-res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.x<0.05,]$log2FoldChange.x, #sign age
       res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.x<0.05,]$log2FoldChange.y, 
       pch=21,bg="red",cex=2.2)
points(-res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.y<0.05,]$log2FoldChange.x, #sign regime
       res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="blue",cex=2.2)
points(-res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.x<0.05 & res.TE.list.age.reg.M$padj.y<0.05,]$log2FoldChange.x,
       res.TE.list.age.reg.M[res.TE.list.age.reg.M$padj.x<0.05 & res.TE.list.age.reg.M$padj.y<0.05,]$log2FoldChange.y, 
       pch=21,bg="orange",cex=2.2) #sign both
text(1.4,1.6,labels = "r = -0.01 (ns)",col = "black",cex=1.9)
text(-1.8,2,labels = "Padj<0.05 for Regime",col = "blue",cex=1.5)
text(-1.8,1.7,labels = "Padj<0.05 for Age",col = "red",cex=1.5)
text(-1.8,1.4,labels = "Padj<0.05 for both",col = "orange",cex=1.5)
dev.off()

#Check expression of TE regulation genes
res.TEregul.table.F = deseq.results(DEseqobj = dds1.F, 
                                    TEvect = TE.regul.genes,
                                    cutoff = 0.05,
                                    type = "table")
res.TEregul.table.M = deseq.results(DEseqobj = dds1.M, 
                                    TEvect = TE.regul.genes,
                                    cutoff = 0.05,
                                    type = "table")

################################################
#Enrichment of TE classes
########################################

row.names(ensembl_casey) = ensembl_casey$TEfam
RNA = 100
DNA = 23 
LTR = 61
nonLTR = 39
TIR = 22

#Enrichment for Age factor:
#class males
table(merge(ensembl_casey,res.TE.list.M.sign$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEclass)
table(merge(ensembl_casey,res.TE.list.M$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEclass)
M = as.table(rbind(c(85, 23), c(RNA-85, DNA-23)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("RNA","DNA"))
fisher.test(M) #ns
M
#subclass males
table(merge(ensembl_casey,res.TE.list.M.sign$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEsubclass)
table(merge(ensembl_casey,res.TE.list.M$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEsubclass)
sum(table(merge(ensembl_casey,res.TE.list.M.sign$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEsubclass)) - 39
M = as.table(rbind(c(53, 55), c(LTR-53, 123-LTR-55)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("LTR","Others"))
M
fisher.test(M) #ns
M = as.table(rbind(c(32, 76), c(nonLTR-32, 123-nonLTR-76)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("nonLTR","Others"))
M
fisher.test(M) #p-value = 0.2369
M = as.table(rbind(c(22, 86), c(TIR-22, 123-TIR-86)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("TIR","Others"))
M
fisher.test(M) #p-value = 0.07045

#class females
table(merge(ensembl_casey,res.TE.list.F.sign$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEclass)
table(merge(ensembl_casey,res.TE.list.F$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEclass)
M = as.table(rbind(c(15, 0), c(RNA-15, DNA-0)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("RNA","DNA"))
M
fisher.test(M) #0.07129, amrginally n.s., retrotransposons change more in females!
#subclass females
table(merge(ensembl_casey,res.TE.list.F.sign$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEsubclass)
table(merge(ensembl_casey,res.TE.list.F$Age_old_vs_young,by="row.names")[,c(2,18,19,20)]$TEsubclass)
M = as.table(rbind(c(10, 5), c(LTR-10, 123-LTR-5)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("LTR","Others"))
sum(M)
M
fisher.test(M) #0.1792
M = as.table(rbind(c(5, 10), c(nonLTR-5, 123-nonLTR-10)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("nonLTR","Others"))
M
sum(M)
fisher.test(M) #p-value = 1
M = as.table(rbind(c(0, 15), c(TIR-0, 123-TIR-15)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("TIR","Others"))
M
fisher.test(M) #p-value = 0.07045, TIR marginally ns


#Enrichment for Regime factor:
ensembl_casey[ensembl_casey$TEfam %in% intersect(row.names(res.TE.list.M.sign$Regime_Sel_vs_Cont),row.names(res.TE.list.F.sign$Regime_Sel_vs_Cont)),]
#all RNA
ensembl_casey[ensembl_casey$TEfam %in% row.names(res.TE.list.M.sign$Regime_Sel_vs_Cont),] 

#males
table(merge(ensembl_casey,res.TE.list.M.sign$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEclass)
table(merge(ensembl_casey,res.TE.list.M$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEclass)
M = as.table(rbind(c(34, 7), c(RNA-34, DNA-7)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("RNA","DNA"))
M
fisher.test(M) #0.8108 no enrichment of RNA or DNA
table(merge(ensembl_casey,res.TE.list.M.sign$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEsubclass)
table(merge(ensembl_casey,res.TE.list.M$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEsubclass)
M = as.table(rbind(c(16, 25), c(nonLTR-16, 123-nonLTR-25)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("nonLTR","others"))
M
fisher.test(M) #0.2256 no enrichment of non-LTR
M = as.table(rbind(c(18, 23), c(LTR-18, 123-LTR-23)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("LTR","others"))
M
fisher.test(M) #0.4454 no enrichment of LTR
M = as.table(rbind(c(6, 35), c(TIR-6, 123-TIR-35)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("LTR","others"))
M
fisher.test(M) #0.6213 no enrichment of TIR

#females
table(merge(ensembl_casey,res.TE.list.F.sign$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEclass)
table(merge(ensembl_casey,res.TE.list.F$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEclass)
M = as.table(rbind(c(26, 1), c(RNA-26, DNA-1)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("RNA","DNA"))
M
fisher.test(M) #0.02475 enriched for retrotransposons!

table(merge(ensembl_casey,res.TE.list.F.sign$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEsubclass)
table(merge(ensembl_casey,res.TE.list.F$Regime_Sel_vs_Cont,by="row.names")[,c(2,18,19,20)]$TEsubclass)
M = as.table(rbind(c(12, 15), c(nonLTR-12, 123-nonLTR-15)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("nonLTR","others"))
M
fisher.test(M) # 0.1589 no enrichment of non-LTR
M = as.table(rbind(c(14, 13), c(LTR-14, 123-LTR-13)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("LTR","others"))
M
fisher.test(M) #0.8303 no enrichment of LTR
M = as.table(rbind(c(1, 26), c(TIR-1, 123-TIR-26)))
dimnames(M) <- list(stat = c("Diff", "NotDiff"),class = c("TIR","others"))
M
fisher.test(M) #0.04305 TIR are significantly underrepresented

########################################
#### EFFECT SIZES: magnitude of expression change across regime and age ####
########################################
#females age
t.test(abs(res.TE.list.F.sign$Age_old_vs_young[res.TE.list.F.sign$Age_old_vs_young$log2FoldChange<0,]$log2FoldChange),
       res.TE.list.F.sign$Age_old_vs_young[res.TE.list.F.sign$Age_old_vs_young$log2FoldChange>0,]$log2FoldChange)
#t = -2.6947, df = 12.982, p-value = 0.0184 #still significant
fem.inc.dec = 
  as.data.frame( rbind(
    cbind("decrease",-res.TE.list.F.sign$Age_old_vs_young[res.TE.list.F.sign$Age_old_vs_young$log2FoldChange<0,]$log2FoldChange),
    cbind("increase",res.TE.list.F.sign$Age_old_vs_young[res.TE.list.F.sign$Age_old_vs_young$log2FoldChange>0,]$log2FoldChange)  ) )   
fem.inc.dec[,2] = numfact(fem.inc.dec[,2])
names(fem.inc.dec) = c("AgeChange","log2FC")
p1a = ggplot(data = fem.inc.dec, aes(x=AgeChange, y=log2FC, group = AgeChange, color = AgeChange)) + geom_boxplot(color="black",outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Age (female)")), y = expression(bold(""))) +
  ylim(c(0,3)) + geom_jitter(width = 0.15,size=2.8,shape=21,col="black",fill="white") + 
  annotate("text", x = 2.3, y=3, label = "P = 0.018",size = 5)
p1a

#males age
t.test(abs(res.TE.list.M.sign$Age_old_vs_young[res.TE.list.M.sign$Age_old_vs_young$log2FoldChange<0,]$log2FoldChange),
       res.TE.list.M.sign$Age_old_vs_young[res.TE.list.M.sign$Age_old_vs_young$log2FoldChange>0,]$log2FoldChange)
#not enough observations
mal.inc.dec = 
  as.data.frame( rbind(
    cbind("decrease",-res.TE.list.M.sign$Age_old_vs_young[res.TE.list.M.sign$Age_old_vs_young$log2FoldChange<0,]$log2FoldChange),
    cbind("increase",res.TE.list.M.sign$Age_old_vs_young[res.TE.list.M.sign$Age_old_vs_young$log2FoldChange>0,]$log2FoldChange)  ) )   
mal.inc.dec[,2] = numfact(mal.inc.dec[,2])
names(mal.inc.dec) = c("AgeChange","log2FC")
p1b = ggplot(data = mal.inc.dec, aes(x=AgeChange, y=log2FC, group = AgeChange, color = AgeChange)) + geom_boxplot(color="black",outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Age (male)")), y = expression(bold(""))) +
  ylim(c(0,3)) + geom_jitter(width = 0.15,size=2.8,shape=21,col="black",fill="white") + ylab(expression(bold("log2FC") ))
p1b

#females regime
t.test(abs(res.TE.list.F.sign$Regime_Sel_vs_Cont[res.TE.list.F.sign$Regime_Sel_vs_Cont$log2FoldChange<0,]$log2FoldChange),
       res.TE.list.F.sign$Regime_Sel_vs_Cont[res.TE.list.F.sign$Regime_Sel_vs_Cont$log2FoldChange>0,]$log2FoldChange) #
#t = 0.18209, df = 8.2091, p-value = 0.8599

#males regime
t.test(abs(res.TE.list.M.sign$Regime_Sel_vs_Cont[res.TE.list.M.sign$Regime_Sel_vs_Cont$log2FoldChange<0,]$log2FoldChange),
       res.TE.list.M.sign$Regime_Sel_vs_Cont[res.TE.list.M.sign$Regime_Sel_vs_Cont$log2FoldChange>0,]$log2FoldChange) # Effect size difference p =  0.0004661
#Males t = -0.064267, df = 14.397, p-value = 0.9496

fem.inc.dec.reg = 
  as.data.frame( rbind(
    cbind("C>S",-res.TE.list.F.sign$Regime_Sel_vs_Cont[res.TE.list.F.sign$Regime_Sel_vs_Cont$log2FoldChange<0,]$log2FoldChange),
    cbind("S>C",res.TE.list.F.sign$Regime_Sel_vs_Cont[res.TE.list.F.sign$Regime_Sel_vs_Cont$log2FoldChange>0,]$log2FoldChange)  ) )   
fem.inc.dec.reg[,2] = numfact(fem.inc.dec.reg[,2])
names(fem.inc.dec.reg) = c("RegimeChange","log2FC")
p2 = ggplot(data = fem.inc.dec.reg, aes(x=RegimeChange, y=log2FC, group = RegimeChange, color = RegimeChange)) + geom_boxplot(color="black",outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Regime (female)")), y = expression(bold(""))) +
  ylim(c(0,3)) + geom_jitter(width = 0.15,size=2.8,shape=21,col="black",fill="white") + 
  annotate("text", x = 2.3, y=3, label = "P = 0.86",size = 5)
p2

mal.inc.dec.reg = 
  as.data.frame( rbind(
    cbind("C>S",-res.TE.list.M.sign$Regime_Sel_vs_Cont[res.TE.list.M.sign$Regime_Sel_vs_Cont$log2FoldChange<0,]$log2FoldChange),
    cbind("S>C",res.TE.list.M.sign$Regime_Sel_vs_Cont[res.TE.list.M.sign$Regime_Sel_vs_Cont$log2FoldChange>0,]$log2FoldChange)  ) )   
mal.inc.dec.reg[,2] = numfact(mal.inc.dec.reg[,2])
names(mal.inc.dec.reg) = c("RegimeChange","log2FC")
p3 = ggplot(data = mal.inc.dec.reg, aes(x=RegimeChange, y=log2FC, group = RegimeChange, color = RegimeChange)) + geom_boxplot(color="black",outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Regime (male)")), y = expression(bold(""))) +
  ylim(c(0,3)) + geom_jitter(width = 0.15,size=2.8,shape=21,col="black",fill="white") + ylab(expression(bold("log2FC"))) +
  annotate("text", x = 2.3, y=3, label = "P = 0.95",size = 5)

p3_p2_p1 = p3 + p2 + p1b + p1a
p3_p2_p1 
ggsave(p3_p2_p1, file="FigS12_degree_age_regime_change_allDmel_diffOrder_revised_MainEffectOnly.pdf", height=9, width=9)


#COMPARISON: Expression vs Genomic Abundance in Carnes2015
#############################################################
#Check if significant TEs in expression are same as in insertion analysis
#############################################################

ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T) 
convert.tab = ensembl_casey[c(1,4)]
names(convert.tab) = c("TEfam", "Name")

female.results = read.table("/Users/danfab/Carnes/RNAseq/TE_expr_MAIN_model_Females.txt",header = T,check.names = F)
female.results$TEfam = row.names(female.results)

male.results = read.table("/Users/danfab/Carnes/RNAseq/TE_expr_MAIN_model_Males.txt",header=T,check.names = F)
male.results$TEfam = row.names(male.results)

carnes.stat.filt =read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = merge(convert.tab,carnes.stat.filt,by="TEfam")
carnes.stat.filt = carnes.stat.filt[order(carnes.stat.filt$Diff_SelCont),]

#subset to significant TEs
carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE",]
nrow(carnes.stat.filt.sign)

female.results.reg = na.omit(female.results[female.results$Regime_Sel_vs_Cont.padj < 0.05 & !is.na(female.results$Regime_Sel_vs_Cont.padj),])
male.results.reg = male.results[male.results$Regime_Sel_vs_Cont.padj < 0.05 & !is.na(male.results$Regime_Sel_vs_Cont.padj),]
female.results.reg.SC = female.results.reg[female.results.reg$Regime_Sel_vs_Cont.log2FoldChange>0,]$TEfam
female.results.reg.CS = female.results.reg[female.results.reg$Regime_Sel_vs_Cont.log2FoldChange<0,]$TEfam
male.results.reg.SC = male.results.reg[male.results.reg$Regime_Sel_vs_Cont.log2FoldChange>0,]$TEfam #5
male.results.reg.CS = male.results.reg[male.results.reg$Regime_Sel_vs_Cont.log2FoldChange<0,]$TEfam
nrow(male.results.reg) #41
nrow(female.results.reg) #27

female.results.age = na.omit(female.results[female.results$Age_old_vs_young.padj < 0.05 & !is.na(female.results$Age_old_vs_young.padj),])
female.results.age.OY = female.results.age[female.results.age$Age_old_vs_young.log2FoldChange>0,]$TEfam
female.results.age.YO = female.results.age[female.results.age$Age_old_vs_young.log2FoldChange<0,]$TEfam
male.results.age = male.results[male.results$Age_old_vs_young.padj < 0.05 & !is.na(male.results$Age_old_vs_young.padj),]
nrow(male.results.age) #108
nrow(female.results.age) #15

#REMOVE MALES, AS INSERTIONS WERE FROM FEMALES
table(carnes.stat.filt.sign[carnes.stat.filt.sign$TEfam %in% female.results.reg$TEfam,]$Diff_SelCont > 0) #15 S>C and 8 CS
carnes.stat.filt.sign[carnes.stat.filt.sign$TEfam %in% female.results.reg.SC,]
carnes.stat.filt.sign[carnes.stat.filt.sign$TEfam %in% female.results.reg.CS,]
table(carnes.stat.filt.sign[carnes.stat.filt.sign$TEfam %in% female.results.reg.SC,]$Diff_SelCont > 0)
table(carnes.stat.filt.sign[carnes.stat.filt.sign$TEfam %in% female.results.reg.CS,]$Diff_SelCont < 0) 

#create table with combined log2fc and abundance
female.results.sign = female.results[female.results$Regime_Sel_vs_Cont.padj < 0.05,]
female.results.sign1 = female.results.sign[,c('TEfam','Regime_Sel_vs_Cont.log2FoldChange',"Regime_Sel_vs_Cont.padj")]
carnes.ins.expr.sign = merge(carnes.stat.filt.sign[,c("TEfam","Name","Mean_Cont","Mean_Sel",'Diff_SelCont', 'log2RatSelCont')],
                             female.results.sign1,by="TEfam",all=T)
write.table(carnes.ins.expr.sign,"TabS17_Abundance_Females_MainEffect.txt", quote=F,sep="\t",row.names = F)

carnes.ins.expr.sign.filter = na.omit(carnes.ins.expr.sign[carnes.ins.expr.sign$Regime_Sel_vs_Cont.padj<0.05,])
nrow(carnes.ins.expr.sign.filter) #23
carnes.ins.expr.sign.filter.sameDirS = carnes.ins.expr.sign.filter[carnes.ins.expr.sign.filter$log2RatSelCont>0 & carnes.ins.expr.sign.filter$Regime_Sel_vs_Cont.log2FoldChange>0,]
carnes.ins.expr.sign.filter.sameDirC = carnes.ins.expr.sign.filter[carnes.ins.expr.sign.filter$log2RatSelCont<0 & carnes.ins.expr.sign.filter$Regime_Sel_vs_Cont.log2FoldChange<0,]

#
######################################################################
#CORRELATION BETWEEN NORMCOUNTS AND INSERTIONS - ALL TEs TOGETHER
######################################################################

norm.counts.F.TE.fin = read.table("/Users/danfab/Carnes/RNAseq/normalized_counts_female_model_allDmel_withInteraction.txt",header=T)
TE.detect = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
nrow(TE.detect)
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_to_flybase  = as.data.frame(cbind(as.character(ensembl_casey$TEfam),gsub("Dmel/","",ensembl_casey$flybase_name)))
names(ensembl_to_flybase ) = c("TEfam","TEname")

#edit genomic insertions tab
TE.detect.melt = melt(TE.detect[,c(1,3:12)], id.vars=c('TEfam'),var='Pop')
TE.detect.melt$Regime = unlist(lapply(strsplit(as.character(TE.detect.melt$Pop),"_"),function(x) x[1]))
names(TE.detect.melt)[3] = "Ins"
TE.detect.melt$InsRank = rank(TE.detect.melt$Ins)
TE.detect.melt[order(TE.detect.melt$InsRank,decreasing = T),]
TE.detect.melt$ID = paste(TE.detect.melt$TEfam,"_",TE.detect.melt$Pop,sep="")
head(TE.detect.melt)

norm.counts.melt = melt(norm.counts.TE.fin,id.vars = names(norm.counts.TE.fin[1:4]))
norm.counts.melt$Pop = paste(norm.counts.melt$Regime,"_",norm.counts.melt$Pop,sep="")
names(norm.counts.melt)[5:6] = c("TEfam","NormCounts")
norm.counts.melt$ID = paste(norm.counts.melt$TEfam,"_",norm.counts.melt$Pop,sep="")
head(norm.counts.melt)

norm.counts.ins = merge(norm.counts.melt,TE.detect.melt,by="ID")
norm.counts.ins = norm.counts.ins[order(norm.counts.ins$Ins, decreasing = T),]
head(norm.counts.ins)
names(norm.counts.ins)[c(2:3,6)] = c("Regime","Pop","TEfam")
norm.counts.ins = norm.counts.ins[order(norm.counts.ins$TEfam),]

#Subset into males/females, old/young, regime
TE.detect.new = TE.detect[,c(1,13,14,15,16,19)]
#TE.detect.new$AbsRank = rank(abs(TE.detect.new$Diff_SelCont))
nrow(norm.counts.ins) #4280
norm.counts.ins = merge(norm.counts.ins,TE.detect.new,by = "TEfam")
nrow(norm.counts.ins) #4280

norm.counts.ins.f = norm.counts.ins[norm.counts.ins$Sex == "F",]
norm.counts.ins.f.y = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "young",]
norm.counts.ins.f.o = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "old",]
norm.counts.ins.f.y.c = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "young" & norm.counts.ins$Regime == "Cont",]
norm.counts.ins.f.y.s = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "young" & norm.counts.ins$Regime == "Sel",]
norm.counts.ins.f.o.c = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "old" & norm.counts.ins$Regime == "Cont",]
norm.counts.ins.f.o.s = norm.counts.ins[norm.counts.ins$Sex == "F" & norm.counts.ins$Age == "old" & norm.counts.ins$Regime == "Sel",]

#correlation (Table S16)
cor.test(norm.counts.ins.f.y.c$NormCounts, norm.counts.ins.f.y.c$Ins,method="spearman" ) #0.7126622 
cor.test(norm.counts.ins.f.y.s$NormCounts, norm.counts.ins.f.y.s$Ins,method="spearman" ) #0.7377342 
cor.test(norm.counts.ins.f.o.c$NormCounts, norm.counts.ins.f.o.c$Ins,method="spearman" ) #0.6840876
cor.test(norm.counts.ins.f.o.s$NormCounts, norm.counts.ins.f.o.s$Ins,method="spearman" ) #0.7129242 
0.01/8
#all correlated, more insertions = more expression!

#average
norm.counts.mean = cbind( as.data.frame(tapply(norm.counts.ins.f.y.c$NormCounts,norm.counts.ins.f.y.c$TEfam,mean)),
                          as.data.frame(tapply(norm.counts.ins.f.y.s$NormCounts,norm.counts.ins.f.y.s$TEfam,mean)),
                          as.data.frame(tapply(norm.counts.ins.f.o.c$NormCounts,norm.counts.ins.f.o.c$TEfam,mean)),
                          as.data.frame(tapply(norm.counts.ins.f.o.s$NormCounts,norm.counts.ins.f.o.s$TEfam,mean)))

names(norm.counts.mean) = c("F_Y_C","F_Y_S", "F_O_C","F_O_S")
norm.counts.mean$mean_all = rowMeans(norm.counts.mean)
norm.counts.mean$TEfam = row.names(norm.counts.mean)
head(norm.counts.mean)
norm.counts.mean = merge(norm.counts.mean, ensembl_to_flybase, by="TEfam") #add common names
norm.counts.mean.diffins = merge(norm.counts.mean,TE.detect.new,by="TEfam")
norm.counts.mean.diffins$MeanIns = rowMeans(norm.counts.mean.diffins[,c("Mean_Cont","Mean_Sel")])
carnes.ins.expr.sign.filter.sameDirS$col = "red"
carnes.ins.expr.sign.filter.sameDirC$col = "blue"
ndir.col = "grey80"
norm.counts.mean.diffins$col = ndir.col
norm.counts.mean.diffins[norm.counts.mean.diffins$TEfam %in% carnes.ins.expr.sign.filter.sameDirS$TEfam,]$col = "red"
norm.counts.mean.diffins[norm.counts.mean.diffins$TEfam %in% carnes.ins.expr.sign.filter.sameDirC$TEfam,]$col = "blue"
#norm.counts.mean.diffins$col = mgsub(c("TRUE","FALSE"), c("red","blue"),norm.counts.mean.diffins$Diff_SelCont>0)
norm.counts.mean.diffins = norm.counts.mean.diffins[order(norm.counts.mean.diffins$col,decreasing=T),]

#Average correlation across all conditions
ins.expr.cor = cor.test(norm.counts.mean.diffins$mean_all,
                        norm.counts.mean.diffins$MeanIns, method="spearman")
ins.expr.cor #0.7209206 , P< 2.2e-16

pdf("FigS15_Ins_vs_Expr_all_TEs_log10_mean_allDmel_REVISED_onlyFemales.pdf", width=6.4, height=6)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2,cex.axis=1.8,mar=c(4.1,4.5,2,2))
plot(x = log10(norm.counts.mean.diffins$MeanIns),
     y = log10(norm.counts.mean.diffins$mean_all),
     pch = 21, bg = "darkgrey", cex = 1.4,
     xlab = expression(paste("Genomic Insertions")), ylab = "Normalized Expression Read Counts",
     axes=F, xlim=c(-0.05,2.2),ylim=c(0,6), frame=T,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(1,seq(0,2.5,0.5),labels = round(10^seq(0,2.5,0.5),1))
axis(2,seq(0,6,2),labels = round(10^seq(0,6,2),2))
text(0.22,max(log10(norm.counts.mean.diffins$mean_all))+0.4, labels = expression(paste(italic(rho), ' = 0.72***')),cex=1.5)
dev.off()


######################################################################
#CORRELATION BETWEEN FOLD CHANGE AND GENOMIC INSERTIONS
######################################################################

TE.expr.f = read.table("/Users/danfab/Carnes/RNAseq/TE_expr_MAIN_model_Females.txt",header=T)
TE.detect = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
TE.expr.f$TEfam = row.names(TE.expr.f)
TE.detect.sign = TE.detect[TE.detect$Bonf %in% c("TRUE","FALSE"),] #take all
TE.expr.Diff.f = merge(TE.detect.sign[,c(1,15,16,19)],TE.expr.f,by="TEfam")
TE.expr.Diff.f$grps = mgsub(c("TRUE","FALSE"),c("S>C","C>S"),TE.expr.Diff.f$Diff_SelCont>0)

cor.test(TE.expr.Diff.f$log2RatSelCont, TE.expr.Diff.f$Regime_Sel_vs_Cont.log2FoldChange,method="spearman") #r = 0.139; p-value = 0.1487

#Colors
TE.expr.f.sign = TE.expr.f[TE.expr.f$Regime_Sel_vs_Cont.padj < 0.05,]
TE.detect.sign = TE.detect[TE.detect$Bonf == TRUE,]
both.sign = cbind(intersect(row.names(TE.expr.f.sign),TE.detect.sign$TEfam),"forestgreen")
expr.sign = cbind(as.character(TE.expr.f.sign[!TE.expr.f.sign$TEfam %in% TE.detect.sign$TEfam,]$TEfam),"magenta")
ins.sign = cbind(as.character(TE.detect.sign[!TE.detect.sign$TEfam %in% TE.expr.f.sign$TEfam,]$TEfam),"orange2")
not.sign = cbind(intersect(TE.expr.f[TE.expr.f$Regime_Sel_vs_Cont.padj >= 0.05,]$TEfam, TE.detect[TE.detect$Bonf == FALSE,]$TEfam),"darkgrey")
sign.cols = rbind(both.sign,expr.sign,ins.sign,not.sign)
colnames(sign.cols) = c("TEfam","col")
sign.cols = sign.cols[sign.cols[,1] %in% intersect(sign.cols[,1],TE.expr.Diff.f$TEfam),]
TE.expr.Diff.f = merge(TE.expr.Diff.f,sign.cols, by = "TEfam")
TE.expr.Diff.f$col = as.character(TE.expr.Diff.f$col)

#output_path = "Carnes/RNAseq/map2ref/"
pdf("Fig4_log2_insertions_vs_log2_RegimeExpr_corr_allDmel_MainEffectsOnly.pdf", width=6, height=5)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=1.6,cex.axis=1.3,mar=c(4.1,5,0.8,2))
plot(TE.expr.Diff.f$log2RatSelCont, TE.expr.Diff.f$Regime_Sel_vs_Cont.log2FoldChange, 
     ylab = "", xlab = "",
     cex=1.8, pch=21, bg = "darkgrey", col ="darkgrey",
     ylim=c(-3,1.2),axes=F,frame=T,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-3,1,1),seq(-3,1,1),lwd=1,line=0,mgp=c(3, .6, 0))
axis(1,seq(-1,1,0.5),seq(-1,1,0.5),lwd=1,line=0,mgp=c(3, .6, 0))

abline(h=0,lty=2)
abline(v=0,lty=2)
points(TE.expr.Diff.f[!TE.expr.Diff.f$col == "darkgrey",]$log2RatSelCont, 
       TE.expr.Diff.f[!TE.expr.Diff.f$col == "darkgrey",]$Regime_Sel_vs_Cont.log2FoldChange,cex=1.8, pch=21, bg = TE.expr.Diff.f[!TE.expr.Diff.f$col == "darkgrey",]$col, col ="black")

title(ylab=expression(paste('Expression: Log'[2], ' FC (S/C)')), line=1.7, cex.lab=1.6)
title(xlab=expression(paste('Genomic Insertions: Log'[2], ' FC (S/C)')), line=2.2, cex.lab=1.6)

text(0.77,1.2, labels = "Females",cex=1.6)
text(0.77,0.88, labels = expression(paste(italic(rho), ' = 0.14 (ns)')),cex=1.3)

text(-0.55,1.2, labels = "Only expression sign.",col = "magenta",cex=1.2)
text(-0.545,0.95,labels = "Only abundance sign.",col = "orange2",cex=1.2)
text(-0.775,0.7,labels = "Both sign.",col = "forestgreen",cex=1.2)
dev.off()

