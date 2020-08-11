#Plot results from genetic drift simulations of TE frequency change
#For Figure S10

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(patchwork)

output_path = "/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output/"
setwd("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output/")

#Robert PlosGen 2015 average frequencies
robert.2015 = read.table("/Users/danfab/extra_files/Robert_Plos2015_AvFreq.txt",header=T)
names(robert.2015)[1] = "Name"
robert.2015$av_freq = numfact(robert.2015$av_freq)
robert.2015 = na.omit(robert.2015)
robert.2015$freq_rank = rank(robert.2015$av_freq)

#Tables with filters and annotation
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
hoedjes.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)

bonf.hoedjes = 0.01 / nrow(hoedjes.stat.filt)
hoedjes.stat.filt$Bonf = hoedjes.stat.filt$PvalFull_Breed < bonf.hoedjes

#Observed ratios - ONLY THOSE ALSO IN SIMULATIONS
int.fab = fabian.stat.filt[gsub('Dmel/','',fabian.stat.filt$flybase_name) %in% robert.2015$Name,]
int.remo = remo.stat.filt[gsub('Dmel/','',remo.stat.filt$flybase_name) %in% robert.2015$Name,]
int.carn = carnes.stat.filt[gsub('Dmel/','',carnes.stat.filt$flybase_name) %in% robert.2015$Name,]
int.hoed = hoedjes.stat.filt[gsub('Dmel/','',hoedjes.stat.filt$flybase_name) %in% robert.2015$Name,]

#cut off
eff_cut = 0.5

#Sign
int.fab.sign = int.fab[int.fab$Bonf == "TRUE" & abs(int.fab$Diff_SelCont) > eff_cut,]
int.remo.sign = int.remo[int.remo$Bonf == "TRUE" & abs(int.remo$Diff_SelCont) > eff_cut,]
int.carn.sign = int.carn[int.carn$Bonf == "TRUE" & abs(int.carn$Diff_SelCont) > eff_cut,]
int.hoed.sign = int.hoed[int.hoed$Bonf_Full == "TRUE" & abs(int.hoed$Diff_PostEarly) > eff_cut,]

nrow(int.carn.sign) / nrow(int.carn) #0.9339623 (96% without)
table(int.carn.sign$Diff_SelCont<0) / nrow(int.carn) #0.7358491 0.1981132 

nrow(int.fab.sign) / nrow(int.fab) #0.6380952 (77% without)
table(int.fab.sign$Diff_SelCont<0) / nrow(int.fab) #0.3428571 0.2952381

nrow(int.hoed.sign) / nrow(int.hoed) #0.4678899 (82% without)
table(int.hoed.sign$Diff_PostEarly<0) / nrow(int.hoed) #0.3669725 0.1009174

nrow(int.remo.sign) / nrow(int.remo) #0.4857143 (69% without)
table(int.remo.sign$Diff_SelCont<0) / nrow(int.remo) #0.3904762 0.0952381


fab.ratio = log2(sum(int.fab.sign$Diff_SelCont > 0) / sum(int.fab.sign$Diff_SelCont < 0))
carn.ratio = log2(sum(int.carn.sign$Diff_SelCont > 0) / sum(int.carn.sign$Diff_SelCont < 0))
hoed.ratio = log2(sum(int.hoed.sign$Diff_PostEarly > 0) / sum(int.hoed.sign$Diff_PostEarly < 0))
remo.ratio = log2(sum(int.remo.sign$Diff_SelCont > 0) / sum(int.remo.sign$Diff_SelCont < 0))

#for result tables of simulations loaded here, see: output_of_15
#Carnes2015
carn.2k = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont2000_Nsel2000_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)
carn.1k = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont1000_Nsel1000_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)
carn.1k.50 = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont1000_Nsel500_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)
carn.2k.50 = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont2000_Nsel1000_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)
carn.1k.75 = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont1000_Nsel250_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)
carn.2k.75 = read.table("/Users/danfab/effective_pop_size/Carnes/genetic_drift_sim/Ncont2000_Nsel500_C850_S170_Reps5_DirCut0.5_Runs5000_PropType.txt",header=T)

#Fabian2018
fab = read.table("/Users/danfab/effective_pop_size/Fabian/genetic_drift_sim/Fab_Ncont300_Nsel300_DirCut0.5_Runs5000_PropType.txt",header=T)
fab.50 = read.table("/Users/danfab/effective_pop_size/Fabian/genetic_drift_sim/Fab_Ncont300_Nsel150_DirCut0.5_Runs5000_PropType.txt",header=T)
fab.75 = read.table("/Users/danfab/effective_pop_size/Fabian/genetic_drift_sim/Fab_Ncont300_Nsel75_DirCut0.5_Runs5000_PropType.txt",header=T)

#Remolina2012
remo.440 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont440_Nsel220_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)
remo.440.50 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont440_Nsel110_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)
remo.440.75 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont440_Nsel55_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)
remo.640 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont640_Nsel320_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)
remo.640.50 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont640_Nsel160_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)
remo.640.75 = read.table("/Users/danfab/effective_pop_size/Remolina/genetic_drift_sim/Ncont640_Nsel80_C80_S50_Reps3_DirCut0.5_Runs5000_PropType.txt",header=T)

#Hoedjes2019
hoed.2k = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont2000_Nsel2000_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)
hoed.2k.50 = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont2000_Nsel1000_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)
hoed.2k.75 = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont2000_Nsel500_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)
hoed.4k = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont4000_Nsel4000_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)
hoed.4k.50 = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont4000_Nsel2000_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)
hoed.4k.75 = read.table("/Users/danfab/effective_pop_size/Hoedjes/genetic_drift_sim/Ncont4000_Nsel1000_C115_S58_Reps12_DirCut0.5_Runs5000_PropType.txt",header =T)

#Define x and y axis limits for plots
ylim.def = c(0,800)
xlim.def = c(-2.5,2.5)

#Fabian2018
#Pvalues
sum(log2(fab$ratio.SC.to.CS) >= fab.ratio) / length(fab$ratio.SC.to.CS) #0.7568
sum(log2(fab.50$ratio.SC.to.CS) >= fab.ratio) / length(fab.50$ratio.SC.to.CS) #.5672
sum(log2(fab.75$ratio.SC.to.CS) >= fab.ratio) / length(fab.75$ratio.SC.to.CS) #0.3254

pdf(file = "FigS10_Histogram_log2_ratio_Fabian2018.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(fab$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(fab$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=fab.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=1.2, y=830, "P = 0.757", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of S and C = 300", bty="n",cex = 1.4,text.col = "black")

hist(log2(fab.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(fab.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=fab.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=1.2, y=830, "P = 0.567", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "50% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")

hist(log2(fab.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],0.5),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(fab.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=fab.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=1.2, y=830, "P = 0.325", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "75% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")
dev.off()

#CARNES2015
#Pvalues
sum(log2(carn.1k$ratio.SC.to.CS) >= carn.ratio) / length(carn.1k$ratio.SC.to.CS) #0
sum(log2(carn.2k$ratio.SC.to.CS) >= carn.ratio) / length(carn.2k$ratio.SC.to.CS) #0
sum(log2(carn.1k.50$ratio.SC.to.CS) >= carn.ratio) / length(carn.1k.50$ratio.SC.to.CS) #0
sum(log2(carn.2k.50$ratio.SC.to.CS) >= carn.ratio) / length(carn.2k.50$ratio.SC.to.CS) #0
sum(log2(carn.1k.75$ratio.SC.to.CS) >= carn.ratio) / length(carn.1k.75$ratio.SC.to.CS) #0
sum(log2(carn.2k.75$ratio.SC.to.CS) >= carn.ratio) / length(carn.2k.75$ratio.SC.to.CS) #0

#1k population size
pdf(file = "FigS10_Histogram_log2_ratio_Carnes2015_1k.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(carn.1k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.1k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of S and C = 1000", bty="n",cex = 1.4,text.col = "black")

hist(log2(carn.1k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.1k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "50% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")

hist(log2(carn.1k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.1k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830,"P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "75% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")
dev.off()

#2k population size
pdf(file = "FigS10_Histogram_log2_ratio_Carnes2015_2k.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(carn.2k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.2k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of S and C = 2000", bty="n",cex = 1.4,text.col = "black")

hist(log2(carn.2k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.2k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "50% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")

hist(log2(carn.2k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(carn.2k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=carn.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.7, y=830,"P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "75% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")
dev.off()


#REMOLINA2012
#Pvalues
sum(log2(remo.640$ratio.SC.to.CS) >= remo.ratio) / length(remo.640$ratio.SC.to.CS) #0
sum(log2(remo.640.50$ratio.SC.to.CS) >= remo.ratio) / length(remo.640.50$ratio.SC.to.CS) #0
sum(log2(remo.640.75$ratio.SC.to.CS) >= remo.ratio) / length(remo.640.75$ratio.SC.to.CS) #0

sum(log2(remo.440$ratio.SC.to.CS) >= remo.ratio) / length(remo.440$ratio.SC.to.CS) #0
sum(log2(remo.440.50$ratio.SC.to.CS) >= remo.ratio) / length(remo.440.50$ratio.SC.to.CS) #0
sum(log2(remo.440.75$ratio.SC.to.CS) >= remo.ratio) / length(remo.440.75$ratio.SC.to.CS) #0

pdf(file = "FigS10_Histogram_log2_ratio_Remo2012_Nc640.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(remo.640$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.640$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ncont = 640; Nsel = 320", bty="n",cex = 1.4,text.col = "black")

hist(log2(remo.640.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.640.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ncont = 640; Nsel = 160", bty="n",cex = 1.4,text.col = "black")

hist(log2(remo.640.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.640.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ncont = 640; Nsel = 80", bty="n",cex = 1.4,text.col = "black")
dev.off()

#440
pdf(file = "FigS10_Histogram_log2_ratio_Remo2012_Nc440.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(remo.440$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.440$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of C = 440; Ne of S = 220", bty="n",cex = 1.4,text.col = "black")

hist(log2(remo.440.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.440.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of C = 440; Ne of S = 110", bty="n",cex = 1.4,text.col = "black")

hist(log2(remo.440.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,axes=F,
     main="",
     ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))))
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(remo.440.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=remo.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.8, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of C = 440; Ne of S = 55", bty="n",cex = 1.4,text.col = "black")

dev.off()

#HOEDJES 2019
#Pvalues
sum(log2(hoed.2k$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.2k$ratio.SC.to.CS) #0
sum(log2(hoed.4k$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.4k$ratio.SC.to.CS) #0
sum(log2(hoed.2k.50$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.2k.50$ratio.SC.to.CS) #0
sum(log2(hoed.4k.50$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.4k.50$ratio.SC.to.CS) #0
sum(log2(hoed.2k.75$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.2k.75$ratio.SC.to.CS) #0
sum(log2(hoed.4k.75$ratio.SC.to.CS) >= hoed.ratio) / length(hoed.4k.75$ratio.SC.to.CS) #0

#2k population size
pdf(file = "FigS10_Histogram_log2_ratio_Hoedjes2015_2k.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(hoed.2k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.2k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of S and C = 2000", bty="n",cex = 1.4,text.col = "black")

hist(log2(hoed.2k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.2k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "50% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")

hist(log2(hoed.2k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.2k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830,"P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "75% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")
dev.off()

#4k population size
pdf(file = "FigS10_Histogram_log2_ratio_Hoedjes2015_4k.pdf",height=10,width=6)
par(mfrow=c(3,1),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,1,2)+.4,cex.lab=1.6)

hist(log2(hoed.4k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.4k$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "Ne of S and C = 4000", bty="n",cex = 1.4,text.col = "black")

hist(log2(hoed.4k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.4k.50$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830, "P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "50% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")

hist(log2(hoed.4k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def, 
     axes=F,ylab="Simulations",
     xlab=expression(bold(paste('log'['2 '],'(S>C / C>S)'))),
     main="")
axis(1,seq(xlim.def[1],xlim.def[2],1),cex.axis=1.4)
axis(2,seq(ylim.def[1],ylim.def[2],200),cex.axis=1.4)
grid(col="grey90",lty=1,lwd=1)
hist(log2(hoed.4k.75$ratio.SC.to.CS),breaks = 25, xlim = xlim.def, ylim=ylim.def,add=T,col = "white")
abline(v=0,lty=3,lwd=2,col="black")
abline(v=hoed.ratio,lty=2,lwd=2,col="red")
box(lty = 1, col = 'black')
legend(x=0.6, y=830,"P < 0.0002", bty="n",cex = 1.4,text.col = "red")
legend("topleft", "75% reduced Ne of S", bty="n",cex = 1.4,text.col = "black")
dev.off()

