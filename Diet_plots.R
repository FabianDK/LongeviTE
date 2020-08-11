#TE abundance analysis for early vs late breeding in Hoedjes2019 using approach 1 

#For Figure S4

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(reshape2)
library(ggplot2)
library(patchwork)

setwd("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/")

#Coverage filter
propCovCut = 0.8 
output_path = "/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/"

#Effectsize cutoff
effcut = 0.3

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)

#Tables with filters and annotation
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
TE_bg = nrow(hoed.stat.filt)
bonf = 0.01 / nrow(hoed.stat.filt)
hoed.cont.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Cont == TRUE & abs(hoed.stat.filt$Diff_PostEarlyCont) > effcut,]
hoed.low.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Low == TRUE & abs(hoed.stat.filt$Diff_PostEarlyLow) > effcut,]
hoed.high.sign = hoed.stat.filt[hoed.stat.filt$Bonf_High == TRUE & abs(hoed.stat.filt$Diff_PostEarlyHigh) > effcut,]
nrow(hoed.low.sign) #86
nrow(hoed.cont.sign) #77
nrow(hoed.high.sign) #67

#All TE families that were tested
te.fam = as.data.frame(hoed.stat.filt$TEfam)
names(te.fam) = "TEfam"

##############################################
#PLOT SIGNIFICANT TE INSERTIONS FOR EACH STUDY - ABSOLUTE DIFFERENCE
##############################################
m = matrix(c(1,2,3,
             4,5,6,
             7,7,7,
             8,8,8),nrow = 4,ncol = 3,byrow = TRUE)

pdf(file = "insertion_diff_log2_Hoedjes_perDiet_sign_Effcut03.pdf",width = 13,height = 10)
par(font.axis = 1, font =1, font.lab =1,cex = 2,cex.axis = 1.9, cex.main = 3.5, font.main=1,
    mar=c(0, 5.5, 5.2, 0.5))

layout(mat = m,
       heights = c(5.5, 5.5, 1),
       widths = c(8, 8, 8))
hoed.low.sign.order = hoed.low.sign[order(hoed.low.sign$Diff_PostEarlyLow),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.low.sign.order$TEsubclass)
plot(1:length(hoed.low.sign.order$Diff_PostEarlyLow),
     hoed.low.sign.order$Diff_PostEarlyLow, 
     xlab = "", 
     ylab= "", 
     main = "Low",
     ylim = c(-10,10),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
# axis(2,seq(-50,30,20),seq(-50,30,20),cex.axis = 2)
title(ylab=expression(paste(delta, "Insertions (S - C)")), line=2.8, cex.lab=2.8)
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-5,35,cex=2, "", bty="n") 

hoed.cont.sign = hoed.cont.sign[order(hoed.cont.sign$Diff_PostEarlyCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.cont.sign$TEsubclass)
plot(1:length(hoed.cont.sign$Diff_PostEarlyCont),
     hoed.cont.sign$Diff_PostEarlyCont, 
     xlab = "", 
     ylab= "", 
     main = "Medium",
     ylim = c(-10,10),xlim=c(0,80),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,11.2,cex=2, "", bty="n") 

hoed.high.sign = hoed.high.sign[order(hoed.high.sign$Diff_PostEarlyHigh),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.high.sign$TEsubclass)
plot(1:length(hoed.high.sign$Diff_PostEarlyHigh),
     hoed.high.sign$Diff_PostEarlyHigh, 
     xlab = "", 
     ylab= "", 
     main = "High",
     ylim = c(-10,10),xlim=c(0,70),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,xaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
abline(h = 0, lty = 2, lwd =1.5)
title(xlab="", line=3, cex.lab=2.2)
legend(-5,4.4,cex=2, "", bty="n") 

par(mar=c(1.2,5.5,4,0.5))
# m = matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)
hoed.low.sign.order = hoed.low.sign[order(hoed.low.sign$log2Rat_PostEarlyLow),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.low.sign.order$TEsubclass)
plot(1:length(hoed.low.sign.order$log2Rat_PostEarlyLow),
     hoed.low.sign.order$log2Rat_PostEarlyLow, 
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5, yaxt="n",
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
title(ylab=expression(paste("log"[2], " (S / C)")), line=2.9, cex.lab=2.8,font=2)
title(xlab="", line=3, cex.lab=2.2)
abline(h = 0, lty = 2, lwd =1.5)
legend(-2,1.5,cex=2, "", bty="n") 

hoed.cont.sign.order = hoed.cont.sign[order(hoed.cont.sign$log2Rat_PostEarlyCont),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.cont.sign.order$TEsubclass)
plot(1:length(hoed.cont.sign.order$log2Rat_PostEarlyCont),
     hoed.cont.sign.order$log2Rat_PostEarlyCont, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),xlim=c(0,80),
     bg=set.col.class,col=set.col.class,
     pch=21,cex = 2.5,
     grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE))
axis(2,seq(-1.5,1.5,0.5),seq(-1.5,1.5,0.5),cex.axis = 1.8)
abline(h = 0, lty = 2, lwd =1.5)
#title(xlab="TEs sorted by effect size", line=3, cex.lab=2.2)
legend(-2,1.5,cex=2, "", bty="n") 

hoed.high.sign.order = hoed.high.sign[order(hoed.high.sign$log2Rat_PostEarlyHigh),]
set.col.class = mgsub(c("TIR","non-LTR","LTR","FB"),c("darkorange","deepskyblue","forestgreen","grey10"),hoed.high.sign.order$TEsubclass)
plot(1:length(hoed.high.sign.order$log2Rat_PostEarlyHigh),
     hoed.high.sign.order$log2Rat_PostEarlyHigh, yaxt="n",
     xlab = "", 
     ylab= "", 
     main = "",
     ylim = c(-1.5,1.5),xlim=c(0,70),
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
       legend = c("TE families sorted by effect size"), cex=2.8, horiz = TRUE,
       x.intersp = -7.5,text.width=c(0.06),bty="n")

par(mar=c(0,2,0,4))
plot(1, type = "n", axes=FALSE, xlab="", ylab="",bty="n")
plot_colors = c("darkorange","forestgreen","deepskyblue","grey10")
legend(x = "top",inset = 0, lty = 0,
       legend = c("TIR","LTR","non-LTR","Helitron"), 
       pt.bg = plot_colors, lwd=1, cex=3, pch=21, horiz = TRUE,
       x.intersp = -0.5,y.intersp = -0.4,text.width=c(0.06,0.06,0.06,0.1),bty="n")
dev.off()
