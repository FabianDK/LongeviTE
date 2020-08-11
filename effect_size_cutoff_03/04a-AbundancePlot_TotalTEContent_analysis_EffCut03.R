#Analysis of total genomic TE content + some plots

#For Figure 1AD, Table S6, Figure S5A and S6

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

#Characterization of single studies
library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(SuperExactTest)
library(FactoMineR)
library(patchwork)
library(grid)
library(gridExtra)
library(lemon)
library(reshape2)

dir.create("/Users/danfab/comp_all/")
output_path = "/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output/"

#TE lookup table
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T,fill=T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name)

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)

#CutOff
effcut = 0.3

#Tables with fitlered for min_coverage and only Bonferroni significant
remo.stat.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE" & abs(remo.stat.filt$Diff_SelCont) > effcut,]
carnes.stat.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE" & abs(carnes.stat.filt$Diff_SelCont) > effcut,]
fabian.stat.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE" & abs(fabian.stat.filt$Diff_SelCont) > effcut,]
hoed.stat.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == "TRUE" & abs(hoed.stat.filt$Diff_PostEarly) > effcut,]

#Background of shared TEs
TE_bg = intersect(hoed.stat.filt$TEfam,intersect(intersect(remo.stat.filt$TEfam, carnes.stat.filt$TEfam),fabian.stat.filt$TEfam)) #background of TEs
length(TE_bg) #103 still

#Total TE content
total_cont.remo = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_total_content.txt",header=T)
total_cont.carn = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_total_content.txt",header=T)
total_cont.fab = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_total_content.txt",header=T)
total_cont.hoed = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_total_content.txt",header=T)

total_cont.remo$Study = "Remolina"
total_cont.carn$Study = "Carnes"
total_cont.fab$Study = "Fabian"
total_cont.hoed$Study = "Hoedjes"

#Add missing columns in other studies, change labelling in Hoedjes, combine all tabs
total_cont.hoed$Breeding = mgsub(c("Early","Postponed"),c("Cont","Sel"),total_cont.hoed$Breeding)
names(total_cont.hoed)[3] = "Regime"
total_cont.hoed = total_cont.hoed[,c(3,2,5:11,4,1)]
total_cont.remo$Diet = NA
total_cont.carn$Diet = NA
total_cont.fab$Diet = NA
total_cont.remo$Type = NA
total_cont.carn$Type = NA
total_cont.fab$Type = NA
total_cont.all = rbind(total_cont.carn,total_cont.fab,total_cont.remo,total_cont.hoed)
total_cont.all$Study = factor(total_cont.all$Study, c('Carnes' ,'Fabian','Hoedjes','Remolina'))

total_cont.all.edit = total_cont.all[,c("Study","Regime","Pop","Diet","Type","Ins","LTR","TIR","nonLTR","RNA","DNA")]
#write.table(x = total_cont.all.edit, file = paste(output_path,"Total_TE_content.txt",sep=""), quote=F, row.names = F, sep="\t")

##############################################
#PLOT SIGNIFICANT TE INSERTIONS FOR EACH STUDY - FIGURE1
##############################################
m = matrix(c(1,2,3,4,
             5,6,7,8,
             9,9,9,9,
             10,10,10,10),nrow = 4,ncol = 4,byrow = TRUE)

pdf(file = paste(output_path,"Fig1_OLD_insertion_diff_log2_All_sign_new.pdf",sep=""),width = 17,height = 10)
par(font.axis = 1, font =1, font.lab =1,cex = 2,cex.axis = 1.9, cex.main = 3.5, font.main=1,
    mar=c(0, 5.5, 5.2, 0.5))

layout(mat = m,
       heights = c(5.5,5.5,1),
       widths = c(8,8,8,8))
carnes.stat.sign.order = carnes.stat.sign[order(carnes.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),carnes.stat.sign.order$TEsubclass)
plot(1:length(carnes.stat.sign.order$Diff_SelCont),
     carnes.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "Carnes2015",
     ylim = c(-60,60),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
# axis(2,seq(-50,30,20),seq(-50,30,20),cex.axis = 2)
title(ylab=expression(paste(delta, "Insertions (S - C)")), line=2.8, cex.lab=2.8)
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-5,35,cex=2, "", bty="n") 

fabian.stat.sign.order = fabian.stat.sign[order(fabian.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),fabian.stat.sign.order$TEsubclass)
plot(1:length(fabian.stat.sign.order$Diff_SelCont),
     fabian.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "Fabian2018",
     ylim = c(-10,10),xlim=c(0,100),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,11.2,cex=2, "", bty="n") 

hoed.stat.sign.order = hoed.stat.sign[order(hoed.stat.sign$Diff_PostEarly),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.stat.sign.order$TEsubclass)
plot(1:length(hoed.stat.sign.order$Diff_PostEarly),
     hoed.stat.sign.order$Diff_PostEarly, 
     xlab = "", 
     ylab= "", 
     main = "Hoedjes2019",
     ylim = c(-10,10),xlim=c(0,100),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,4.4,cex=2, "", bty="n") 

remo.stat.sign.order = remo.stat.sign[order(remo.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),remo.stat.sign.order$TEsubclass)
plot(1:length(remo.stat.sign.order$Diff_SelCont),
     remo.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "Remolina2012",
     ylim = c(-4,4),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,4.4,cex=2, "", bty="n") 

par(mar=c(1.2,5.5,4,0.5))
# m = matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)
carnes.stat.sign.order = carnes.stat.sign[order(carnes.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),carnes.stat.sign.order$TEsubclass)
plot(1:length(carnes.stat.sign.order$log2RatSelCont),
     carnes.stat.sign.order$log2RatSelCont, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5, yaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
title(ylab=expression(paste("log"[2], " (S / C)")), line=2.9, cex.lab=2.8,font=1)
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-2,1.5,cex=2, "", bty="n") 

fabian.stat.sign.order = fabian.stat.sign[order(fabian.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),fabian.stat.sign.order$TEsubclass)
plot(1:length(fabian.stat.sign.order$log2RatSelCont),
     fabian.stat.sign.order$log2RatSelCont, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),xlim=c(0,100),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
abline(h = 0, lty = 2, lwd =1.5)
#title(xlab="TEs sorted by effect size", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

hoed.stat.sign.order = hoed.stat.sign[order(hoed.stat.sign$log2Rat_PostEarly),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.stat.sign.order$TEsubclass)
plot(1:length(hoed.stat.sign.order$log2Rat_PostEarly),
     hoed.stat.sign.order$log2Rat_PostEarly, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),xlim=c(0,100),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

remo.stat.sign.order = remo.stat.sign[order(remo.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),remo.stat.sign.order$TEsubclass)
plot(1:length(remo.stat.sign.order$log2RatSelCont),
     remo.stat.sign.order$log2RatSelCont, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),
     xlim=c(0,80),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

par(mar=c(0,2,0,2))
plot(1, type = "n", axes=FALSE, xlab="", ylab="",bty="n")
plot_colors = c("darkorange","forestgreen","deepskyblue","grey10")
legend(x = "top",inset = 0, lty = 0,
       legend = c("TEs sorted by effect size"), cex=2.8, horiz = TRUE,
       x.intersp = -1.7,text.width=c(0.06),bty="n")

par(mar=c(0,2,0,4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="",bty="n")
plot_colors = c("darkorange","forestgreen","deepskyblue","grey10")
legend(x = "top",inset = 0, lty = 0,
       legend = c("TIR","LTR","non-LTR","Helitron"), 
       pt.bg = plot_colors, lwd=1, cex=3, pch=21, horiz = TRUE,
       x.intersp = -0.5,y.intersp = -0.4,text.width=c(0.06,0.06,0.06,0.1),bty="n")
dev.off()



#FIGURE 1 REVISED - LOG2FC
dev.off()
pdf(file = paste(output_path,"Fig1_insertion_log2_All_sign_revised.pdf",sep=""),width = 17,height = 4.2)
par(mar=c(2.4,5.5,2,0.2),mfrow=c(1,4),cex.axis = 2.2)
carnes.stat.sign.order = carnes.stat.sign[order(carnes.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),carnes.stat.sign.order$TEsubclass)
plot(1:length(carnes.stat.sign.order$log2RatSelCont),
     carnes.stat.sign.order$log2RatSelCont, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1,1),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5, yaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,round(seq(-1,1,0.5),2),round(seq(-1,1,0.5),2),cex.axis = 2.2, las = 1)
title(ylab="", line=2.9, cex.lab=2.8,font=1) #expression(paste("log"[2], " (S / C)"))
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-2,1.5,cex=2, "", bty="n") 

fabian.stat.sign.order = fabian.stat.sign[order(fabian.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),fabian.stat.sign.order$TEsubclass)
plot(1:length(fabian.stat.sign.order$log2RatSelCont),
     fabian.stat.sign.order$log2RatSelCont, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1,1),xlim=c(0,70),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1,1,0.5),seq(-1,1,0.5),cex.axis = 2.2,las = 2)
abline(h = 0, lty = 2, lwd =1.5)
#title(xlab="TEs sorted by effect size", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

hoed.stat.sign.order = hoed.stat.sign[order(hoed.stat.sign$log2Rat_PostEarly),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.stat.sign.order$TEsubclass)
plot(1:length(hoed.stat.sign.order$log2Rat_PostEarly),
     hoed.stat.sign.order$log2Rat_PostEarly, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-0.6,0.6),xlim=c(0,52),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,round(seq(-0.6,0.6,0.2),2),round(seq(-0.6,0.6,0.2),2),cex.axis = 2.2,las = 2)
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

remo.stat.sign.order = remo.stat.sign[order(remo.stat.sign$log2RatSelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),remo.stat.sign.order$TEsubclass)
plot(1:length(remo.stat.sign.order$log2RatSelCont),
     remo.stat.sign.order$log2RatSelCont, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-0.6,0.6),
     xlim=c(0,55),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,round(seq(-0.6,0.6,0.2),2),round(seq(-0.6,0.6,0.2),2),cex.axis = 2.2,las = 2)
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

dev.off()


#FIGURE S5 REVISED - dIns
pdf(file = paste(output_path,"FigS5_insertion_diff_All_sign_revised.pdf",sep=""),width = 17,height = 4.6)
par(mar=c(2.4,5.5,2,0.2),mfrow=c(1,4),cex.axis = 2.2,cex.main=2.8)
carnes.stat.sign.order = carnes.stat.sign[order(carnes.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),carnes.stat.sign.order$TEsubclass)
plot(1:length(carnes.stat.sign.order$Diff_SelCont),
     carnes.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-60,60),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="t",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
# axis(2,seq(-50,30,20),seq(-50,30,20),cex.axis = 2)
title(ylab=expression(paste(delta, "Insertions (S - C)")), line=2.8, cex.lab=2.8)
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-5,35,cex=2, "", bty="n") 

fabian.stat.sign.order = fabian.stat.sign[order(fabian.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),fabian.stat.sign.order$TEsubclass)
plot(1:length(fabian.stat.sign.order$Diff_SelCont),
     fabian.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-10,10),xlim=c(0,70),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="t",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,11.2,cex=2, "", bty="n") 

hoed.stat.sign.order = hoed.stat.sign[order(hoed.stat.sign$Diff_PostEarly),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.stat.sign.order$TEsubclass)
plot(1:length(hoed.stat.sign.order$Diff_PostEarly),
     hoed.stat.sign.order$Diff_PostEarly, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-10,10),xlim=c(0,52),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="t",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,4.4,cex=2, "", bty="n") 

remo.stat.sign.order = remo.stat.sign[order(remo.stat.sign$Diff_SelCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","Helitron"),c("darkorange","deepskyblue","forestgreen","grey10"),remo.stat.sign.order$TEsubclass)
plot(1:length(remo.stat.sign.order$Diff_SelCont),
     remo.stat.sign.order$Diff_SelCont, 
     xlab = "", 
     ylab= "", 
     main = "",
     xlim = c(0,52),
     ylim = c(-4,4),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="t",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,4.4,cex=2, "", bty="n") 
dev.off()

###############################################
#Check difference in total content
###############################################

#interaction model
intmodel_output = NULL
for (i in 3:6) {
  int = total_cont.all[,c(1,9,i)]
  int.av = matrix(tapply(int[,3],list(int$Regime),mean),ncol = 2)
  int.diff = int.av[2] - int.av[1]
  int.log2FC = log2(int.av[2] / int.av[1])
  model = lm(int[,3] ~ Study+Regime+Study*Regime,data=int) #positive means more in selected
  
  stud.stat = matrix(c(anova(model)[["Df"]][1], round(anova(model)[["F value"]][1],2), anova(model)[["Pr(>F)"]][1]),ncol=3,byrow=T)
  reg.stat = matrix(c(anova(model)[["Df"]][2], round(anova(model)[["F value"]][2],2), anova(model)[["Pr(>F)"]][2]),ncol=3,byrow=T)
  int.stat = matrix(c(anova(model)[["Df"]][3], round(anova(model)[["F value"]][3],2), anova(model)[["Pr(>F)"]][3]),ncol=3,byrow=T)
  df.resid = anova(model)[["Df"]][4]
  
  int.tab.stat = cbind(names(int)[3], 
                       int.av, int.diff, int.log2FC,
                       stud.stat,reg.stat,int.stat,
                       summary(model)$adj.r.squared,
                       df.resid)
  intmodel_output = rbind(intmodel_output,int.tab.stat)
  
}
intmodel_output = as.data.frame(intmodel_output)
names(intmodel_output) = c("Type","Cont","Sel","DiffSelCont","log2SelCont","study_Df","study_Fval","study_P","regime_Df","regime_Fval","regime_P","inter_Df","inter_Fval","inter_P","adj_R_squared","Df_resid")
for (i in 2:ncol(intmodel_output)) {
  intmodel_output[,i] = numfact(intmodel_output[,i])
}
intmodel_output #TIR and LTR significant
p.adjust(intmodel_output$study_P,method="fdr") < 0.05 #all
p.adjust(intmodel_output$regime_P[1:4],method="fdr") < 0.05 #only TIR
p.adjust(intmodel_output$inter_P[1:4],method="fdr") < 0.05 #only TIR

#Tukey tests for significant LTR and TIR
total.model = aov(LTR ~ Study+Regime+Regime*Study, data = total_cont.all)
summary(total.model)
TukeyHSD(total.model)

total.model = aov(TIR ~ Study+Regime+Regime*Study, data = total_cont.all)
summary(total.model) 
TukeyHSD(total.model)

#Save stats table
write.table(intmodel_output,paste(output_path,"Total_Content_Summary_Stats.txt",sep=""),row.names = F,quote=F,sep="\t")

#PLOT TOTAL CONTENT
#Change factors
total_cont.all.Ins = melt(total_cont.all,id.vars = c("Study","Regime"),measure.vars = "Ins")
total_cont.all.TIR = melt(total_cont.all,id.vars = c("Study","Regime"),measure.vars = "TIR") 
total_cont.all.LTR = melt(total_cont.all,id.vars = c("Study","Regime"),measure.vars = "LTR") 
total_cont.all.nonLTR = melt(total_cont.all,id.vars = c("Study","Regime"),measure.vars = "nonLTR") 

#plot
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
                      legend.key.size = unit(0.8, "cm"),
                      legend.text = element_text(size=15),
                      legend.title = element_text(size=18,face="bold"),
                      legend.position="right",
                      panel.background = element_blank(),
                      panel.border = element_rect(fill = NA, colour="black", size = 1),
                      panel.grid.major = element_line(),
                      panel.grid.minor = element_blank())

p1 = ggplot(data = total_cont.all.Ins, aes(x=Study, y=value, fill=Regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = expression(bold("")), y = expression(bold("Total Insertions"))) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")

p2 = ggplot(data = total_cont.all.TIR, aes(x=Study, y=value,fill=Regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("TIR Insertions"))) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")

p3 = ggplot(data = total_cont.all.LTR, aes(x=Study, y=value,fill=Regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("LTR Insertions"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")

p4 = ggplot(data = total_cont.all.nonLTR, aes(x=Study, y=value, fill=Regime)) + theme_bw() + geom_boxplot(outlier.shape = NA,position="dodge2") + THEME_DEF_BOX + 
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("non-LTR Insertions"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p4

p.all = grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2, position='right')
ggsave(p.all, file=paste(output_path,"Fig1_and_FigS6_insertion_Classes_new.pdf",sep=""), height=9, width=14)
