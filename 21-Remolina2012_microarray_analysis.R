#Analysis of microarray data from Remolina et al. 2012

#For Figure 5

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

dir.create("/Users/danfab/Remolina/microarray")
setwd("/Users/danfab/Remolina/microarray")

source("https://bioconductor.org/biocLite.R")

library(GPL15486)
library(GEOmetadb)
library(lme4)
library(lmerTest)

#Load TE regulation genes table
piRNA_genes = read.table("/Users/danfab/extra_files/piRNA_genes.txt",header=T)
epigen_genes = read.table("/Users/danfab/extra_files/epigenetics_genes.txt",header=T)
transpos_genes = read.table("/Users/danfab/extra_files/transposition_genes_edit.txt")
transpos_genes.private = transpos_genes[!transpos_genes$V2 %in% piRNA_genes$GeneId & !transpos_genes$V2 %in% epigen_genes$GeneId,]
TE.regul.genes = union(transpos_genes.private$V1,union(piRNA_genes$FBGN,epigen_genes$FBGN))
  
#load expressio ndata
gset <- getGEO("GSE38106", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL15486", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


#Download expression set
gset = getGEO("GSE38106", GSEMatrix =TRUE, AnnotGPL=TRUE,destdir="/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/remolina_all/microarray/")

if (length(gset) > 1) idx = grep("GPL15486", attr(gset, "names")) else idx = 1
gset = gset[[idx]]

ncol(exprs(gset))
# 120
nrow(exprs(gset))
# 16637

boxplot(exprs(gset))
#if means are 

regime = gset$characteristics_ch1.1
tissue = gset$characteristics_ch1.2
pop = gset$characteristics_ch1.3
age = gset$characteristics_ch1.4
rep = gset$characteristics_ch1.5
samplenames = sampleNames(gset)

table(age)
table(regime)
table(tissue)

colnames(fData(gset)) #gene names not here
fData(gset)$ORF
fData(gset)

#Expression data
table(exprs(gset) > 0)
row.names(expr.dat) #these are already the genes
table(duplicated(row.names(expr.dat))) #none duplciated



#########################
#multiple probes per gene (different transcripts or replicates) - average them
#########################

exprs(gset)[1:2,1:2]
expr.probe.annot = merge(fData(gset)[,2:3],exprs(gset),by="row.names")
row.names(expr.probe.annot) = expr.probe.annot$Row.names
expr.probe.annot = expr.probe.annot[,-c(1,2)]
head(expr.probe.annot)
expr.probe.mean = aggregate(expr.probe.annot[,2:ncol(expr.probe.annot)],
          by=list(expr.probe.annot$CLONE_ID),mean) 
nrow(expr.probe.mean) #13995
row.names(expr.probe.mean) = expr.probe.mean$Group.1
expr.probe.mean1 = expr.probe.mean[,-1]
head(expr.probe.mean1)

t(exprs(gset))[1:2,1:2]
variables = cbind(gsub("tissue: ","",as.character(tissue)),
      gsub("age (day): ","",as.character(age),fixed = T),
      mgsub(c("group: "," populations"," (for longer reproductive lifespans)"),c("","",""),as.character(regime),fixed=T),
      gsub("population: ","",as.character(pop),fixed = T),
      gsub("biological replicate: ","",as.character(rep),fixed = T))

expr.tab = as.data.frame(cbind(variables,t(expr.probe.mean1)))
names(expr.tab)[1:5] = c("Tissue","Age","Regime","Pop","Rep")
names(expr.tab)[1:10]
expr.tab$Age = numfact(expr.tab$Age)
class(expr.tab$Age)
class(expr.tab$Rep)

#transform to numeric
for (i in 6:ncol(expr.tab)){
  expr.tab[,i] = numfact(expr.tab[,i])
}

bighead(expr.tab,1,10,1,10)
expr.tab$ID
table(do.call(paste, c(expr.tab[c(1:4)], sep="_"))) #2 replicates of each
expr.tab = cbind(ID=do.call(paste, c(expr.tab[c(1:4)], sep="_")),
                expr.tab)

#test if replicates differ in the 60 samples
rep.names = as.character(unique(expr.tab$ID))
nrow(expr.tab[expr.tab$ID == rep.names[1],]) #2

rep.stat.test = NULL
for (i in 1:length(rep.names)){
  testtab = expr.tab[expr.tab$ID == rep.names[i],]
  pval = t.test(as.numeric(testtab[1,7:ncol(testtab)]),
                as.numeric(testtab[2,7:ncol(testtab)]))$p.value
  rep.int = cbind(rep.names[i],pval)
  rep.stat.test = rbind(rep.stat.test,rep.int) 
}
rep.stat.test
rep.stat.test1 = as.data.frame(rep.stat.test)
rep.stat.test1$pval = numfact(rep.stat.test1$pval)
rep.stat.test1$FDR = p.adjust(rep.stat.test1$pval,method="fdr")
rep.stat.test1 = rep.stat.test1[order(rep.stat.test1$pval),]
table(rep.stat.test1$pval<0.01) #true for 6, 54 equal
table(rep.stat.test1$FDR<0.01) #true for 4, 56 equal

######################################################################
#STATS
######################################################################

#Simplify and just take mean of replicates
expr.tab.mean = aggregate(expr.tab[,7:ncol(expr.tab)],by=list(expr.tab$Tissue,expr.tab$Age,expr.tab$Regime,expr.tab$Pop), mean)
nrow(expr.tab.mean) #60
nrow(expr.tab) #120
expr.tab.mean[1:10,1:10]
names(expr.tab.mean)[1:4] = c("Tissue","Age","Regime","Pop")

############
#with tissues separate 
expr.tab.ab = expr.tab[expr.tab$Tissue == "Abdomen",]
expr.tab.hd = expr.tab[expr.tab$Tissue == "Head",]

expr.tab.ab$pop_age = paste(expr.tab.ab$Pop,expr.tab.ab$Age,sep="_")
expr.tab.hd$pop_age = paste(expr.tab.hd$Pop,expr.tab.hd$Age,sep="_")

expr.pvals = NULL
end.loop = ncol(expr.tab.ab)-1 #-1 becuase last column is pop_age
for (i in 7:end.loop){ 
  gene.name = names(expr.tab.ab)[i]
  print(i)
  
  model.hd = lmer(expr.tab.hd[,i] ~ Age + Regime + Age*Regime + (1|pop_age), data=expr.tab.hd)
  anova( model.hd)
  F.hd = anova(model.hd)[["F value"]][1:3]
  numDF.hd = anova(model.hd)[["NumDF"]][1:3]
  denDF.hd = anova(model.hd)[["DenDF"]][1:3]
  pval.hd = anova(model.hd)[["Pr(>F)"]][1:3]
  
  model.ab = lmer(expr.tab.ab[,i] ~ Age + Regime + Age*Regime + (1|pop_age), data=expr.tab.ab)
  anova( model.ab)
  F.ab= anova(model.ab)[["F value"]][1:3]
  numDF.ab = anova(model.ab)[["NumDF"]][1:3]
  denDF.ab = anova(model.ab)[["DenDF"]][1:3]
  pval.ab = anova(model.ab)[["Pr(>F)"]][1:3]
  
  fin.pval = cbind(gene.name, 
                   F.hd[1],numDF.hd[1],denDF.hd[1],pval.hd[1],
                   F.hd[2],numDF.hd[2],denDF.hd[2],pval.hd[2],
                   F.hd[3],numDF.hd[3],denDF.hd[3],pval.hd[3],
                   F.ab[1],numDF.ab[1],denDF.ab[1],pval.ab[1],
                   F.ab[2],numDF.ab[2],denDF.ab[2],pval.ab[2],
                   F.ab[3],numDF.ab[3],denDF.ab[3],pval.ab[3])
  expr.pvals = rbind(fin.pval,
                     expr.pvals)
}
expr.pvals = as.data.frame(expr.pvals)
nrow(expr.pvals) #13995
head(expr.pvals)
head(expr.pvals)

names(expr.pvals) = c("Gene",
                      "F_Age_H","numDF_Age_H","denomDF_Age_H","Pval_Age_H",
                      "F_Regime_H","numDF_Regime_H","denomDF_Regime_H","Pval_Regime_H",
                      "F_AgeReg_H","numDF_AgeReg_H","denomDF_AgeReg_H","Pval_AgeReg_H",
                      "F_Age_A","numDF_Age_A","denomDF_Age_A","Pval_Age_A",
                      "F_Regime_A","numDF_Regime_A","denomDF_Regime_A","Pval_Regime_A",
                      "F_AgeReg_A","numDF_AgeReg_A","denomDF_AgeReg_A","Pval_AgeReg_A")
ncol(expr.pvals)

for (i in c(5,9,13,17,21,25)){
  expr.pvals[,i] = numfact(expr.pvals[,i])
  expr.pvals = cbind(expr.pvals,p.adjust(expr.pvals[,i],method="fdr"))
}
names(expr.pvals) = c("FBGN",
                      "F_Age_H","numDF_Age_H","denomDF_Age_H","Pval_Age_H",
                      "F_Regime_H","numDF_Regime_H","denomDF_Regime_H","Pval_Regime_H",
                      "F_AgeReg_H","numDF_AgeReg_H","denomDF_AgeReg_H","Pval_AgeReg_H",
                      "F_Age_A","numDF_Age_A","denomDF_Age_A","Pval_Age_A",
                      "F_Regime_A","numDF_Regime_A","denomDF_Regime_A","Pval_Regime_A",
                      "F_AgeReg_A","numDF_AgeReg_A","denomDF_AgeReg_A","Pval_AgeReg_A",
                      "Age_H_fdr","Regime_H_fdr","AgeReg_H_fdr",
                      "Age_A_fdr","Regime_A_fdr","AgeReg_A_fdr")
head(expr.pvals)

for (i in 2:ncol(expr.pvals)){
  expr.pvals[,i] = numfact(expr.pvals[,i])
}

#save output (on dryad)
write.table(expr.pvals,"Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep.txt",quote=F,sep="\t",col.names = T,row.names = F) 

#Check how many genes significant

#FDR<0.05
table(expr.pvals$Age_H_fdr < 0.05) #463
table(expr.pvals$Age_A_fdr < 0.05) #2449
table(expr.pvals$Regime_H_fdr < 0.05) #4
table(expr.pvals$Regime_A_fdr < 0.05) #0
table(expr.pvals$AgeReg_H_fdr < 0.05) #0
table(expr.pvals$AgeReg_A_fdr < 0.05) #0

#P<0.01
table(expr.pvals$Pval_Regime_H< 0.01) #491
table(expr.pvals$Pval_Regime_A < 0.01) #8
table(expr.pvals$Pval_AgeReg_H< 0.01) #15
table(expr.pvals$Pval_AgeReg_A < 0.01) #11

intersect(TE.regul.genes,expr.pvals[expr.pvals$Pval_AgeReg_H< 0.01,]$FBGN) #none
intersect(TE.regul.genes,expr.pvals[expr.pvals$Pval_AgeReg_A< 0.01,]$FBGN) #none

#save and then update FBgns in FlyBase (use upload/convert tool)
write.table(expr.pvals[expr.pvals$Pval_Regime_H<0.01,], "Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_REG01_H.txt",quote=F,sep="\t",col.names = T,row.names = F)
write.table(expr.pvals[expr.pvals$Age_H_fdr < 0.05,], "Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_AGE05fdr_H.txt",quote=F,sep="\t",col.names = T,row.names = F)
write.table(expr.pvals[expr.pvals$Age_A_fdr < 0.05,], "Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_AGE05fdr_A.txt",quote=F,sep="\t",col.names = T,row.names = F)

#Get candidates and update IDs using FlyBase
expr.pvals[expr.pvals$Pval_Regime_A<0.01,]$FBGN #all up to date
intersect(expr.pvals[expr.pvals$Pval_Regime_A<0.01,]$FBGN,TE.regul.genes) #0

######################################################################
#Check if TE regulation genes are significant
######################################################################

#Load table with updated gene IDs
reg01H = read.table("/Users/danfab/Remolina/microarray/Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_REG01_H_convert.txt",header=F,fill=F)
age05fdrH = read.table("/Users/danfab/Remolina/microarray/Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_AGE05fdr_H_convert.txt")
age05fdrA = read.table("/Users/danfab/Remolina/microarray/Age_Reg_AgeReg_PopAgeRandom_onReps_TissueSep_AGE05fdr_A_convert.txt",fill=T)
age05fdrA[171,]
age05fdrA[age05fdrA$V3 == "ID",]

#Regime
intersect(reg01H$V3,piRNA_genes$FBGN) #no overlap Regime to piRNA genes
intersect(reg01H$V3,epigen_genes$FBGN) #"FBgn0016756" "FBgn0035829" - two genes
intersect(reg01H$V3,transpos_genes.private$V1) #none
expr.pvals[expr.pvals$FBGN  %in% c("FBgn0016756","FBgn0035829"), ]
reg01H.epigen = expr.tab.hd[,names(expr.tab.hd) %in% c("Regime","FBgn0016756","FBgn0035829")]
tapply(reg01H.epigen$FBgn0016756,reg01H.epigen$Regime,mean) #up in selected
tapply(reg01H.epigen$FBgn0035829,reg01H.epigen$Regime,mean) #down in selected

#Age Head
luniq(intersect(age05fdrH$V3,TE.regul.genes)) #5
intersect(age05fdrH$V3,piRNA_genes$FBGN) #"FBgn0086251" "FBgn0034255" "FBgn0004400" "FBgn0016034"
intersect(age05fdrH$V3,transpos_genes.private$V1) #none tranposition
intersect(age05fdrH$V3,epigen_genes$FBGN) #"FBgn0035608"

#Head sign:
expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,piRNA_genes$FBGN))]
aggregate(expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,piRNA_genes$FBGN))],by=list(expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,piRNA_genes$FBGN))]$Age), mean)
#all 4 piRNA genes in head are increasing with Age
expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,epigen_genes$FBGN))]
aggregate(expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,epigen_genes$FBGN))],by=list(expr.tab.hd[,names(expr.tab.hd) %in% c("Age",intersect(age05fdrH$V3,piRNA_genes$FBGN))]$Age), mean)
#the single epigenetic gene in head is increasing with Age

luniq(intersect(as.character(age05fdrA$V3),TE.regul.genes)) #24
intersect(as.character(age05fdrA$V3),piRNA_genes$FBGN) 
#"FBgn0263974" "FBgn0027499" "FBgn0086908" "FBgn0016034" "FBgn0051202" "FBgn0002962"
intersect(age05fdrA$V3,transpos_genes.private$V1) #"FBgn0037707" "FBgn0037844" "FBgn0041775"
intersect(age05fdrA$V3,epigen_genes$FBGN) 
#[1] "FBgn0033233" "FBgn0003612" "FBgn0038551" "FBgn0004655" "FBgn0013263" "FBgn0027951" "FBgn0004401" "FBgn0263144" "FBgn0020309"
#[10] "FBgn0015396" "FBgn0003598" "FBgn0086908" "FBgn0015805" "FBgn0027567" "FBgn0260397" "FBgn0030301"
luniq(intersect(age05fdrA$V3,epigen_genes$FBGN) ) #16

#check sign of coefficient to know if increasing or decreasing ABDOMEN
genes.age.ab.defence = unique(union(names(expr.tab.ab[,names(expr.tab.ab) %in% age05fdrA[age05fdrA$V3 %in% intersect(age05fdrA$V3,epigen_genes$FBGN),]$V1,]), 
             names(expr.tab.ab[,names(expr.tab.ab) %in% age05fdrA[age05fdrA$V3 %in% intersect(age05fdrA$V3,piRNA_genes$FBGN),]$V1,])))
age05fdrA$V3 == age05fdrA$V1
genes.age.ab.defence =  age05fdrA[ as.character(age05fdrA$V3) %in% intersect(as.character(age05fdrA$V3),TE.regul.genes),]
genes.age.ab.defence[duplicated(genes.age.ab.defence$V3),] #FBgn0260397
genes.age.ab.defence[genes.age.ab.defence$V3 == "FBgn0260397",]
genes.age.ab.defence[!as.character(genes.age.ab.defence$V1) == as.character(genes.age.ab.defence$V3),]
nrow(genes.age.ab.defence)

coeff.summary=NULL
for (i in 1:nrow(genes.age.ab.defence )){
  model.ab = lmer(expr.tab.ab[,names(expr.tab.ab) == as.character(genes.age.ab.defence$V1[i])] ~ Age + Regime + Age*Regime + (1|pop_age), data=expr.tab.ab)
  coeffs = cbind(as.character(genes.age.ab.defence$V1[i]),
                 as.character(genes.age.ab.defence$V3[i]),
                 as.character(genes.age.ab.defence$V4[i]),
                 matrix(summary(model.ab)$coefficients[c(1,2)],ncol=2)) #pos sign
  coeff.summary = rbind(coeff.summary, coeffs)
}
coeff.summary = as.data.frame(coeff.summary)
coeff.summary
#all up except from FBgn0051202
luniq(coeff.summary$V2) #24 after update
#23 up, 1 down with age


