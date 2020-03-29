#Analysis of genome-wide Watterson's theta calculated using PoPoolation
#For Table S7

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")
library(plotrix)
library(ggplot2)

#Load tables: available on dryad at https://doi.org/10.5061/dryad.s7h44j13r
remo = read.table("Remolina_all_theta_100kb.new",header=F,na.strings = "na")
carnes = read.table("Carnes_all_theta_100kb.new",header=F,na.strings = "na")
fabian = read.table("Fabian_all_theta_100kb.new",header=F,na.strings = "na")
hoedjes = read.table("Hoedjes_all_theta_100kb.new",header=F,na.strings = "na")

#Edit tables
remo$Regime = gsub("_[1-9]","",remo$V6)
carnes$Regime = gsub("_[A-Z][1-9]","",carnes$V6)
fabian$Regime = gsub("_.*","",fabian$V6)

breed = mgsub(c('[A-Z]E','[A-Z]P'),c('Cont_','Sel_'),hoedjes$V6)
regime = matrix(unlist(strsplit(x = as.character(breed),split = "_")),ncol=2,byrow=T)[,1]
diet = mgsub(c('C[A-Z][1-9]','L[A-Z][1-9]','H[A-Z][1-9]'),c('Cont','Low','High'),hoedjes$V6)

hoedjes$regime = regime
hoedjes$diet = diet

names(remo) = c("Chrom","WindowPos","SNPs","CovProp","theta","Pop","Regime")
names(carnes) = c("Chrom","WindowPos","SNPs","CovProp","theta","Pop","Regime")
names(fabian) = c("Chrom","WindowPos","SNPs","CovProp","theta","Pop","Regime")
names(hoedjes) = c("Chrom","WindowPos","SNPs","CovProp","theta","Pop","Regime","Diet")

#Select chromosomes
remo = remo[remo$Chrom %in% c("X","2L","2R","3L","3R","4"),]
carnes = carnes[carnes$Chrom %in% c("X","2L","2R","3L","3R","4"),]
fabian = fabian[fabian$Chrom %in% c("X","2L","2R","3L","3R","4"),]
hoedjes = hoedjes[hoedjes$Chrom %in% c("X","2L","2R","3L","3R","4"),]

head(remo)
head(carnes)
head(fabian)
head(hoedjes)

#REMOLINA
round(tapply(X = remo$theta,remo$Regime,mean,na.rm = T),5)
tapply(X = remo$theta,remo$Regime,std.error,na.rm = T)

model = lm(theta ~ Regime/Pop + Chrom + Regime,data=remo)
summary(model)
anova(model)

p.remo = ggplot(data = remo, aes(x=Pop, y=theta, group = Pop, color = Regime)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Population")), y = expression(bold("Nucleotide diversity theta"))) + 
  scale_color_manual(values=c("blue","red")) + 
  geom_jitter(width = 0.1, pch = 21,size=0.5) + 
  theme(axis.text.x = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(vjust=0.5,size=26),
        axis.title.x = element_text(vjust=0.5,size=26))
p.remo

#CARNES
round(tapply(X = carnes$theta,carnes$Regime,mean,na.rm = T),5)
tapply(X = carnes$theta,carnes$Regime,std.error,na.rm = T)

model = lm(theta ~ Regime/Pop + Chrom + Regime,data=carnes)
summary(model)
anova(model)

p.carnes = ggplot(data = carnes, aes(x=Pop, y=theta, group = Pop, color = Regime)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Population")), y = expression(bold("Nucleotide diversity theta"))) + 
  scale_color_manual(values=c("blue","red")) + 
  geom_jitter(width = 0.1, pch = 21,size=0.5) + 
  theme(axis.text.x = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(vjust=0.5,size=26),
        axis.title.x = element_text(vjust=0.5,size=26))
p.carnes


#FABIAN
round(tapply(X = fabian$theta,fabian$Regime,mean,na.rm = T),5)
tapply(X = fabian$theta,fabian$Regime,std.error,na.rm = T)

fabian$Pop = factor(fabian$Pop, c("Cont_Ra", "Cont_Rb", "Sel_La", "Sel_Lb", "Sel_2La",  "Sel_2Lb"))
model = lm(theta ~ Regime/Pop + Chrom + Regime,data=fabian)
summary(model)
anova(model)

p.fabian = ggplot(data = fabian, aes(x=Pop, y=theta, group = Pop, color = Regime)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Population")), y = expression(bold("Nucleotide diversity theta"))) + 
  scale_color_manual(values=c("blue","red")) + 
  geom_jitter(width = 0.1, pch = 21,size=0.5) + 
  theme(axis.text.x = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(vjust=0.5,size=26),
        axis.title.x = element_text(vjust=0.5,size=26))
p.fabian

#HOEDJES
round(tapply(X = hoedjes$theta,hoedjes$Regime,mean,na.rm = T),5)
tapply(X = hoedjes$theta,hoedjes$Regime,std.error,na.rm = T)

#all together analyzed
model = lm(theta ~ Chrom + Diet + Diet * Regime + Regime ,data=hoedjes)
summary(model)
anova(model)

p.hoedjes = ggplot(data = hoedjes, aes(x=Pop, y=theta, group = Pop, color = Regime)) + 
  geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("Population")), y = expression(bold("Nucleotide diversity theta"))) + 
  scale_color_manual(values=c("blue","red")) + 
  geom_jitter(width = 0.1, pch = 21,size=0.5) + 
  theme(axis.text.x = element_text(size=14, colour = "black"),
        axis.text.y = element_text(size=14, colour = "black"),
        axis.title.y = element_text(vjust=0.5,size=26),
        axis.title.x = element_text(vjust=0.5,size=26))
p.hoedjes
