#Quality checks, check difference in replicates, and filtering of RNA-seq read counts from Carnes et al. 2015 (PLoS One)

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(DESeq2)

setwd("Carnes/RNAseq/")

#Load tables
counts = read.table("/Users/danfab/Carnes/RNAseq/readCounts_carnes.txt",header=T)
ribosomal_genes = read.table("/Users/danfab/extra_files/ribosomal_gene_IDs.txt")
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_casey.dmel = ensembl_casey[ensembl_casey$species == "Dmel",]
nrow(ensembl_casey.dmel) #126
TE.detect = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
length(TE.detect$TEfam) #112

head(counts)
names(counts)

#Edit read counts table
namelist = strsplit(names(counts)[8:ncol(counts)] ,"_")
names(counts)[8:ncol(counts)] = unlist(lapply(namelist, function(x) paste(x[1:3],collapse = "_")))
head(counts)

#Edit chromosome column so that chrom is only showing once
chrom.uniq = lapply(X = strsplit(as.character(counts$Chr),";"),FUN = unique)
counts$Chr = unlist(chrom.uniq)

#Same for strands
strand.uniq = lapply(X = strsplit(as.character(counts$Strand),";"),FUN = unique)
strand.uniq[[3933]] = paste(strand.uniq[[3933]],collapse = "/")
counts$Strand = unlist(strand.uniq)

counts.edit = counts[,c(1,2,5,6,7,8:ncol(counts))]
head(counts.edit)

############################################################
counts.edit.detect = counts.edit[rowSums(counts.edit[,6:ncol(counts.edit)]) >=400,]
nrow(counts.edit.detect) #15094
counts.edit.detect = counts.edit.detect[order(rowSums(counts.edit.detect[,6:ncol(counts.edit.detect)]),decreasing = T),]

#check ribosomal gene counts
rib.genes.detect = counts.edit.detect[counts.edit.detect$Geneid %in% ribosomal_genes$V1,]
nrow(rib.genes.detect) #12
colSums(rib.genes.detect[,6:ncol(rib.genes.detect)])
mean(colSums(rib.genes.detect[,6:ncol(rib.genes.detect)])) #34471.86
prop.rib.genes = colSums(rib.genes.detect[,6:ncol(rib.genes.detect)]) / colSums(counts.edit.detect[,6:ncol(counts.edit.detect)])
min(prop.rib.genes) #0.00052
max(prop.rib.genes) #0.00700
#only small number of ribosomal genes expressed!

#check Y chrom genes
as.data.frame(colSums(counts.edit[counts.edit$Chr == "Y",][,6:ncol(counts.edit[counts.edit$Chr == "Y",])]))
as.data.frame(colSums(counts.edit.detect[counts.edit.detect$Chr == "Y",][,6:ncol(counts.edit.detect[counts.edit.detect$Chr == "Y",])]))
#females only few wrongly mapped reads, males have many as expected

#Filter for TEs and major chromosomal arms
counts.edit.detect.filter = counts.edit.detect[counts.edit.detect$Chr %in% c("2L","2R","3L","3R","X","4",as.character(ensembl_casey$TEfam)),]

#Check strands
table(counts.edit.detect.filter$Strand) / nrow(counts.edit.detect.filter) #about 50/50 strands

counts.edit.detect.filter.m = melt(counts.edit.detect.filter,id.var=names(counts.edit.detect.filter)[1:5])
names(counts.edit.detect.filter.m) = c("Geneid","Chr","Strand","Length","gene_symbol","Sample","ReadCounts")
counts.edit.detect.filter.m$log10Counts = log10(counts.edit.detect.filter.m$ReadCounts)
counts.edit.detect.filter.m[counts.edit.detect.filter.m$log10Counts == -Inf,]
counts.edit.detect.filter.m$log10Counts = gsub(pattern = "-Inf",replacement = NA,counts.edit.detect.filter.m$log10Counts)
counts.edit.detect.filter.m.narm = counts.edit.detect.filter.m[!is.na(counts.edit.detect.filter.m$log10Counts),]
counts.edit.detect.filter.m.narm$log10Counts = as.numeric(counts.edit.detect.filter.m.narm$log10Counts)
counts.edit.detect.filter.m.narm.TE = counts.edit.detect.filter.m.narm[counts.edit.detect.filter.m.narm$Geneid %in% ensembl_casey$TEfam,]
head(counts.edit.detect.filter.m.narm.TE)

#transpose so that genes as columns
counts.edit.detect.filter.t = t(counts.edit.detect.filter[,-c(2:5)])
colnames(counts.edit.detect.filter.t) = counts.edit.detect.filter.t[1,]
counts.edit.detect.filter.t = as.data.frame(counts.edit.detect.filter.t[-1,])

labmeta = matrix(unlist(strsplit(mgsub(c("B","O"),c("Cont_B","Sel_O"),row.names(counts.edit.detect.filter.t)),"_")),ncol=4,byrow=T)
labmeta[,3] = mgsub(c("W1","W5"),c("young","old"),labmeta[,3])
labmeta[,4] = mgsub(c("1","2"),c("",""),labmeta[,4])

counts.edit.detect.filter.t = as.data.frame(cbind(labmeta,
                                                  row.names(counts.edit.detect.filter.t),
                                                  counts.edit.detect.filter.t))
names(counts.edit.detect.filter.t)[1:5] = c("Regime","Pop","Age","Sex","ID")
for (i in 6:ncol(counts.edit.detect.filter.t)){
  counts.edit.detect.filter.t[,i] = numfact(counts.edit.detect.filter.t[,i])
}
bighead(counts.edit.detect.filter.t,1,10,1,10)
nrow(counts.edit.detect.filter.t) #80
table(colnames(counts.edit.detect.filter.t) == "AY561850")

#Some more checks
strand.dist.plot = ggplot(data = counts.edit.detect.filter.m, 
                          aes(x=Sample, y=ReadCounts, group = Sample, color = Strand,fill=Strand)) + geom_boxplot() + theme_bw()
strand.dist.plot

#TE plots with counts
TE.dist.plot = ggplot(data = counts.edit.detect.filter.m.narm.TE,
                      aes(x=reorder(Chr, log10Counts, FUN = median),
                          y=log10Counts, group = Chr, color = Chr)) + geom_boxplot(color="black") + theme_bw() +
  theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)) +
  labs(x = expression(bold("TE")), y = expression(bold("log10 ReadCounts")))
TE.dist.plot

TE.dist.plot1 = ggplot(data = counts.edit.detect.filter.m.narm.TE,
                      aes(x=reorder(Chr, ReadCounts, FUN = median),
                          y=ReadCounts, group = Chr, color = Chr)) + geom_boxplot(color="black") + theme_bw() +
  theme(legend.position = "none") + theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)) +
  labs(x = expression(bold("TE")), y = expression(bold("ReadCounts")))
TE.dist.plot1

#PCA
pca.rawcounts = PCA(counts.edit.detect.filter.t[,6:ncol(counts.edit.detect.filter.t)], graph = F)
summary(pca.rawcounts) #%var explained, 63.700   13.225    3.517

dimtab = as.data.frame(cbind(counts.edit.detect.filter.t[,1:5],pca.rawcounts$ind$coord))
for (i in 6:ncol(dimtab)){
  dimtab[,i] = numfact(dimtab[,i])
}
class(dimtab$Dim.1)

colors.for.plot = mgsub(c("Sel_young","Sel_old","Cont_young","Cont_old"), c("orange","red","lightblue","blue"),paste(dimtab$Regime,dimtab$Age,sep="_"))
shape.for.plot = as.numeric(mgsub(c("F","M"), c("21","22"),dimtab$Sex))

# pdf("PCA_RNAseq_rawCounts.pdf", width=7, height=7)
par(font=2,font.axis=2,font.lab=2)
plot(x = dimtab$Dim.1, y = -dimtab$Dim.2, bg = colors.for.plot, cex = 1.5, pch = shape.for.plot,
     xlab = "Dimension 1 (63.7%)", ylab = "Dimension 2 (13.2%)",
     xlim=c(-180,180),
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
text(x = dimtab$Dim.1+20, y = -dimtab$Dim.2,dimtab$ID,cex=0.4,col="grey")
# legend("topleft",legend = unique(paste(dimtab$Regime,dimtab$Age,dimtab$Sex,sep="_")), 
#        pt.bg = c("lightblue","lightblue","blue","blue","orange","orange","red","red"), 
#        pch = rep(c(21,22),8),
#        cex = 1, pt.cex = 1.5,
#        x.intersp = 0.7)
dev.off()
#Biggest difference between males and females

#############
#Check if sign difference between replicates using the whole TE distribution
tab.repCheck = counts.edit.detect.filter.t
tab.repCheck[,5] = mgsub(c("F1","F2","M1","M2"),c("F","F","M","M"),tab.repCheck$ID)
bighead(tab.repCheck,1,10,1,10)

repCheck.summary = NULL
for (i in 1:length(unique(tab.repCheck$ID))){
  subset.repCheck = tab.repCheck[tab.repCheck$ID == unique(tab.repCheck$ID)[i],]
  
  t.model = t.test(as.numeric(subset.repCheck[1,6:ncol(subset.repCheck)]),
                   as.numeric(subset.repCheck[2,6:ncol(subset.repCheck)]))
  
  cor.model = cor.test(as.numeric(subset.repCheck[1,6:ncol(subset.repCheck)]),
                   as.numeric(subset.repCheck[2,6:ncol(subset.repCheck)]),method="pearson")
  
  wilcox = wilcox.test(as.numeric(subset.repCheck[1,6:ncol(subset.repCheck)]),
              as.numeric(subset.repCheck[2,6:ncol(subset.repCheck)]),paired=F, alternative="two.sided")
  
  stat.text = paste("t = ",round(t.model$statistic,1), ", df = ",round(t.model$parameter,1), ", P = ",round(t.model$p.value,3),sep="")
  
  axis.limit = max(ceiling(max(log10(as.numeric(subset.repCheck[1,6:ncol(subset.repCheck)])))),
                   ceiling(max(log10(as.numeric(subset.repCheck[2,6:ncol(subset.repCheck)])))))
  #pdf(paste("replicates_var/",unique(tab.repCheck$ID)[i],".pdf",sep=""), width=6, height=6)
  par(font=2,font.lab=2,font.axis=2)
  plot(log10(as.numeric(subset.repCheck[1,6:ncol(subset.repCheck)])),
       log10(as.numeric(subset.repCheck[2,6:ncol(subset.repCheck)])),
       main = unique(tab.repCheck$ID)[i], 
       xlab = "log10 Counts (Rep 1)",
       ylab = "log10 Counts (Rep 2)",pch=21,bg="grey",cex=1.2,
       xlim=c(0,axis.limit),
       ylim=c(0,axis.limit))
  abline(0,1,col="red",lty=2,lwd=2)
  legend("topleft",legend = stat.text,x.intersp = 0)
  dev.off()
  
  int.tab = cbind(as.character(unique(tab.repCheck$ID)[i]),
                  round(t.model$p.value,4),
                  round(wilcox$p.value,4),
                  round(cor.model$estimate,4),
                  round(cor.model$p.value,5))
  repCheck.summary = rbind(repCheck.summary,int.tab)
  
}
repCheck.summary = as.data.frame(repCheck.summary)
repCheck.summary
names(repCheck.summary) = c("Sample","Pval","WilcoxP","PearsCor","CorPval")
repCheck.summary$Pval = numfact(repCheck.summary$Pval)
repCheck.summary$WilcoxP = numfact(repCheck.summary$WilcoxP)
repCheck.summary$PearsCor = numfact(repCheck.summary$PearsCor)
repCheck.summary$CorPval = numfact(repCheck.summary$CorPval)
repCheck.summary$FDR = p.adjust(repCheck.summary$Pval,"fdr")
repCheck.summary$FDRcor = p.adjust(repCheck.summary$CorPval,"fdr")
repCheck.summary$FDRwilcox = p.adjust(repCheck.summary$WilcoxP,"fdr")
#all replicates highly correlated

# pdf(paste("replicates_var/FDR_replicates.pdf",sep=""), width=8, height=5)
par(font=2,font.lab=2,font.axis=2,cex.lab=1.2,cex.axis=1.2)
plot(1:length(sort(repCheck.summary$FDR)),
     sort(repCheck.summary$FDR), cex=1.2,
     ylab="FDR between 2 replicates", xlab = "Samples (sorted by FDR)",
     pch=21,bg="grey",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = T))
abline(h=0.01,col="red",lty=2)
abline(h=0.05,col="orange",lty=2)
dev.off()
#
table(repCheck.summary$Pval<0.05) #8
table(repCheck.summary$FDR<0.05) #7 sign
table(repCheck.summary$Pval<0.01/40) #5 signat Bonferroni
table(repCheck.summary$WilcoxP<0.01/40) #19
table(repCheck.summary$CorPval<0.01/40) #40 all sign correlated
minmax(repCheck.summary$PearsCor) #0.9512 0.9999
minmax(repCheck.summary$CorPval)
table(repCheck.summary$FDRwilcox<0.05) #26
repCheck.summary[repCheck.summary$FDR<0.05,] #7
repCheck.summary[repCheck.summary$FDR<0.01,] #5
repCheck.summary$FDR001 = repCheck.summary$FDR<0.01
repCheck.summary$Bonf001 = repCheck.summary$Pval<0.01/40
minmax(repCheck.summary$PearsCor) #0.9512 to 0.9999
#most replicates not significantly different, all of them highly correlated

#To simplify analysis, continue with analyzing mean between replicates
tab.mean = aggregate(counts.edit.detect.filter.t[,6:ncol(counts.edit.detect.filter.t)],
                         by=list(counts.edit.detect.filter.t$Regime,
                                 counts.edit.detect.filter.t$Pop,
                                 counts.edit.detect.filter.t$Age,
                                 counts.edit.detect.filter.t$Sex),mean)
names(tab.mean)[1:4] = c("Regime","Pop","Age","Sex")
allTEs = colnames( tab.mean)[grep(pattern = "FBgn",x = colnames(tab.mean),invert=T)]
write.table(tab.mean,"rawCounts_filtered_Chrom_min400_allDmel_mean.txt",quote=F,row.names=F,sep="\t")

#PCA using means between replicates
#PCA for TE families
tab.mean.TE = tab.mean[,names(tab.mean) %in% TE.detect$TEfam]
head(tab.mean.TE)
tab.mean.TE = cbind(tab.mean[,1:4],tab.mean.TE)
pca.rawcounts.mean = PCA(tab.mean.TE[,5:ncol(tab.mean.TE)], graph = F)
summary(pca.rawcounts.mean)

dimtab = as.data.frame(cbind(tab.mean.TE[,1:4],pca.rawcounts.mean$ind$coord))
for (i in 6:ncol(dimtab)){
  dimtab[,i] = numfact(dimtab[,i])
}
class(dimtab$Dim.1)

colors.for.plot = mgsub(c("Sel_young","Sel_old","Cont_young","Cont_old"), c("orange","red","lightblue","blue"),paste(dimtab$Regime,dimtab$Age,sep="_"))
shape.for.plot = as.numeric(mgsub(c("F","M"), c("21","22"),dimtab$Sex))

#pdf("PCA_RNAseq_rawCounts_mean_TEs.pdf", width=7, height=7)
par(font=2,font.axis=2,font.lab=2)
plot(x = dimtab$Dim.1, y = -dimtab$Dim.2, bg = colors.for.plot, cex = 1.5, pch = shape.for.plot,
     xlab = "Dimension 1 (54.8%)", ylab = "Dimension 2 (8%)",
     xlim=c(-20,20),
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
text(x = dimtab$Dim.1+20, y = -dimtab$Dim.2,dimtab$ID,cex=0.4,col="grey")
legend("topleft",legend = c("Cont_young_F","Cont_young_M","Cont_old_F","Cont_old_M",
                            "Sel_young_F","Sel_young_M","Sel_old_F","Sel_old_M"), 
       pt.bg = c("lightblue","lightblue","blue","blue","orange","orange","red","red"), 
       pch = rep(c(21,22),8),
       cex = 1, pt.cex = 1.5,
       x.intersp = 0.7)
dev.off()

#PCA for genes
tab.mean.genes = tab.mean[,!names(tab.mean) %in% TE.detect$TEfam]
pca.rawcounts.mean = PCA(tab.mean.genes[,5:ncol(tab.mean.genes)], graph = F)
summary(pca.rawcounts.mean)

dimtab = as.data.frame(cbind(tab.mean.genes[,1:4],pca.rawcounts.mean$ind$coord))
for (i in 6:ncol(dimtab)){
  dimtab[,i] = numfact(dimtab[,i])
}
class(dimtab$Dim.1)

colors.for.plot = mgsub(c("Sel_young","Sel_old","Cont_young","Cont_old"), c("orange","red","lightblue","blue"),paste(dimtab$Regime,dimtab$Age,sep="_"))
shape.for.plot = as.numeric(mgsub(c("F","M"), c("21","22"),dimtab$Sex))

#pdf("PCA_RNAseq_rawCounts_mean_Genes.pdf", width=7, height=7)
par(font=2,font.axis=2,font.lab=2)
plot(x = dimtab$Dim.1, y = -dimtab$Dim.2, bg = colors.for.plot, cex = 1.5, pch = shape.for.plot,
     xlab = "Dimension 1 (67.7%)", ylab = "Dimension 2 (12.4%)",
     xlim=c(-125,125),
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
text(x = dimtab$Dim.1+20, y = -dimtab$Dim.2,dimtab$ID,cex=0.4,col="grey")
legend("topleft",legend = c("Cont_young_F","Cont_young_M","Cont_old_F","Cont_old_M",
                            "Sel_young_F","Sel_young_M","Sel_old_F","Sel_old_M"), 
       pt.bg = c("lightblue","lightblue","blue","blue","orange","orange","red","red"), 
       pch = rep(c(21,22),8),
       cex = 1, pt.cex = 1.5,
       x.intersp = 0.7)
dev.off()

