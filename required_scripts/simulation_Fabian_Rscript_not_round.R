#Simulate TE frequency change over generations in selected and control populations of Fabian et al. 2018 (Evol. Letters)

#Please change folder names and paths.
dir.create("/Users/danfab/effective_pop_size/")
setwd("/Users/danfab/effective_pop_size/")

#
library(pBrackets)
library(reshape2)
library(ggplot2)
library(Rfast)

#
#Dataset to select detected TEs
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
consistent.SC.CS.TEs = c(15,8)

#Robert PlosGen 2015 average frequencies
robert.2015 = read.table("/Users/danfab/extra_files/Robert_Plos2015_AvFreq.txt",header=T)
names(robert.2015)[1] = "Name"
robert.2015$av_freq = numfact(robert.2015$av_freq)
robert.2015 = na.omit(robert.2015)
#nrow(robert.2015) #114
robert.2015$freq_rank = rank(robert.2015$av_freq)

#####################
#CUSTOM FUNCTIONS + PLOT SETTINGS
#####################

#Transform factors to numbers function
numfact = function(x){
  transform=as.numeric(as.character(x))
  return(transform)
}

#ggplot theme
THEME_DEF = theme(axis.title.y = element_text(vjust=0.5,size=30),
                  axis.title.x = element_text(vjust=-0.2,size=28),
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
                  panel.grid.major = element_line(colour="grey85", size = (0.2)),
                  panel.grid.minor = element_blank())

#Multiple gsub function
mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}


#MANY INDEPENDENT ALLELES
drifTE = function(N, pi, gen, R){
  df = pi
  popsize = NULL
  gVarDF = NULL
  for (i in 1:gen){
    pi = vapply(pi, function(x) rbinom(1, N*2, prob = x), as.integer(1L)) / (N*2) #NEED TO TAKE x2 BECAUSE DIPLOID
    df=cbind(df,pi)
    #popsize=rbind(popsize,N)
    #gVar = q*(1-q) / N^R
    #gVarDF=rbind(gVarDF,gVar)
    N = round(N^R)
  }
  return(df)
}


#####################
#args to pass
#####################

args = commandArgs(trailingOnly=TRUE)
popsizeC = as.integer(args[1])  #e.g. 640 NEED TO TAKE x2 BECAUSE DIPLOID
popsizeS = as.integer(args[2])  #e.g. 320 NEED TO TAKE x2 BECAUSE DIPLOID
dir.cutoff = as.character(args[3]) ##needs to be >=0.5 - i.e. if a TE has a chance of >0.5 to be S>C, it is counted as S>C
dir.cutoff = as.numeric(dir.cutoff)
run.numb = as.integer(args[4]) #e.g. 100

#####################
#set output names
#####################

output_name = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,".txt",sep="")
output_name2 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_PropType.txt",sep="")
output_name_pdf.All = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_All.pdf",sep="")
output_name_pdf.Examp = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Example.pdf",sep="")
output_name_pdf2 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Boxplot.pdf",sep="")
output_name_pdf3 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Direction.pdf",sep="")
output_name_pdf4 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_PropType.pdf",sep="")
output_name_pdf5 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Ratios.pdf",sep="")
output_name_pdf6 = paste("Fab_Ncont",popsizeC,"_Nsel",popsizeS,"_DirCut",dir.cutoff,"_Runs",run.numb,"_ConsReps.pdf",sep="")

#set R = 1
R = 1 #this is not used

#####################
#RUN SCRIPT
#####################

#Select TEs also found in Fabian
int = robert.2015[robert.2015$Name %in% gsub('Dmel/','',fabian.stat.filt$flybase_name),]
af = sort(int$av_freq)
length(af)

#Observed S>C and C>S, and ns
int2 = fabian.stat.filt[gsub('Dmel/','',fabian.stat.filt$flybase_name) %in% robert.2015$Name,]
int3 = int2[int2$Bonf == "TRUE",]
observed.prop.SC = sum(int3$Mean_Sel > int3$Mean_Cont) / nrow(int2)
observed.prop.CS = sum(int3$Mean_Sel < int3$Mean_Cont) / nrow(int2)
observed.prop.Equal = (nrow(int2) - nrow(int3)) / nrow(int2)


#Show example of genetic drift through generations per TE
TEsim_S1 = drifTE(popsizeS, af, 144, R=1)
TEsim_S2 = drifTE(popsizeS, af, 144, R=1)
TEsim_S3 = drifTE(popsizeS, af, 145, R=1)
TEsim_S4 = drifTE(popsizeS, af, 148, R=1)

TEsim_C1 = drifTE(popsizeC, af, 293, R=1)
TEsim_C2 = drifTE(popsizeC, af, 310, R=1)

#EXAMPLE FOR ONE RUN
pdf(file = output_name_pdf.Examp,width = 4 * 3.4, height=6)
par(mfrow=c(2,4),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,2,1))
plot(0:144,TEsim_S1[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:144,TEsim_S1[i,],type="l",col=i)
}

plot(0:144,TEsim_S2[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:144,TEsim_S2[i,],type="l",col=i)
}

plot(0:145,TEsim_S3[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:145,TEsim_S3[i,],type="l",col=i)
}

plot(0:148,TEsim_S4[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:148,TEsim_S4[i,],type="l",col=i)
}

plot(0:293,TEsim_C1[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:293,TEsim_C1[i,],type="l",col=i)
}
plot(0:310,TEsim_C2[1,],
     type="l",ylim=c(0,1),col=1,xlab="",
     ylab="",
     main=paste("Selected N=",popsizeS),xaxt="t")
mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
grid(lty=1,col="grey95")
for (i in 2:length(af)){
  lines(0:310,TEsim_C2[i,],type="l",col=i)
}

dev.off()

#Simulate how many S>C and C>S TEs expected by chance: Proportion of TEs more in S, C, or equal
#Mean across replicates - DON'T ROUND - then check how many TEs are S>C, C>S, or equal
paste("STARTING SIMS FOR PROPORTION OF S>C AND C>S")
simtab = NULL
for (i in 1:run.numb){
  fabian.sel1_2 = replicate(2,drifTE(N = popsizeS, pi = af, gen = 144, R = 1))
  fabian.sel3 = replicate(1,drifTE(N = popsizeS, pi = af, gen = 145, R = 1))
  fabian.sel4 = replicate(1,drifTE(N = popsizeS, pi = af, gen = 148, R = 1))
  
  fabian.cont1 = replicate(1,drifTE(N = popsizeC, pi = af, gen = 293, R = 1))
  fabian.cont2 = replicate(1,drifTE(N = popsizeC, pi = af, gen = 310, R = 1))
  
  fabian.sel1_2 = lapply(seq(dim(fabian.sel1_2)[3]), function(x) fabian.sel1_2[ , , x])
  fabian.sel3 = lapply(seq(dim(fabian.sel3)[3]), function(x) fabian.sel3[ , , x])
  fabian.sel4 = lapply(seq(dim(fabian.sel4)[3]), function(x) fabian.sel4[ , , x])
  
  fabian.cont1 = lapply(seq(dim(fabian.cont1)[3]), function(x) fabian.cont1[ , , x])
  fabian.cont2 = lapply(seq(dim(fabian.cont2)[3]), function(x) fabian.cont2[ , , x])
  
  #select last generation and put regimes together
  fab.sel = c(lapply(fabian.sel1_2, function(x) x[,ncol(x)]), 
                 lapply(fabian.sel3, function(x) x[,ncol(x)]),
                 lapply(fabian.sel4, function(x) x[,ncol(x)]))
  fab.cont = c(lapply(fabian.cont1, function(x) x[,ncol(x)]), 
                 lapply(fabian.cont2, function(x) x[,ncol(x)]))
  
  #Average diff
  fab.sel.mean = Reduce("+",fab.sel) / 4
  fab.cont.mean = Reduce("+",fab.cont) / 2
  
  S.lastgen = fab.sel.mean
  C.lastgen = fab.cont.mean
  
  propSC = sum(S.lastgen > C.lastgen) / length(af)
  propCS = sum(S.lastgen < C.lastgen) / length(af)
  propEqual = sum(S.lastgen == C.lastgen) / length(af)
  ratio.SC.to.CS = propSC / propCS #larger 1 = more S>C
  ratio.Equal.to.SC.CS = propEqual / (propSC / propCS) #larger 1 = more Equal
  
  #Consistent diff across reps:
  sel.lastgen.tab = matrix(unlist(fab.sel),ncol=length(af),byrow=T)
  cont.lastgen.tab = matrix(unlist(fab.cont),ncol=length(af),byrow=T)
  
  SC.consist = sum(colMins(sel.lastgen.tab,value=T) > colMaxs(cont.lastgen.tab,value=T))
  CS.consist = sum(colMins(sel.lastgen.tab,value=T) < colMaxs(cont.lastgen.tab,value=T))
  
  #Create tabs
  int = cbind(propSC, propCS, propEqual, 
              ratio.SC.to.CS ,ratio.Equal.to.SC.CS,
              SC.consist, CS.consist)
  simtab = rbind(simtab,int)
}
simtab = as.data.frame(simtab)

paste("Average Prop. S>C TEs:", round(mean(simtab$propSC),3), "|| sd =",round(sd(simtab$propSC),2))
paste("Average Prop. C>S TEs:", round(mean(simtab$propCS),3), "|| sd =",round(sd(simtab$propCS),2))
paste("Average Prop. Equal TEs:", round(mean(simtab$propEqual),3), "|| sd =",round(sd(simtab$propEqual),2))
paste("")
paste("Average TEs Consistent S>C:", round(mean(simtab$SC.consist),3), "|| sd =",round(sd(simtab$SC.consist),2))
paste("Average TEs Consistent C>S:", round(mean(simtab$CS.consist),3), "|| sd =",round(sd(simtab$CS.consist),2))
paste("")
paste("Prop of Sims with more or equal Consistent S>C:", sum(simtab$SC.consist>=consistent.SC.CS.TEs[1]) / run.numb)
paste("Prop of Sims with more or equal Consistent C>S:", sum(simtab$CS.consist>=consistent.SC.CS.TEs[2]) / run.numb)
paste("Prop of Sims with Consistent S>C more than consistent C>S:", sum(simtab$SC.consist > simtab$CS.consist) / run.numb)
paste("")

#Write tab
write.table(simtab, output_name2,quote=F,sep="\t",row.names = F)

pval.SC = sum(simtab$propSC >= observed.prop.SC) / run.numb
pval.CS = sum(simtab$propCS >= observed.prop.CS) / run.numb
pval.Equal = sum(simtab$propEqual >= observed.prop.Equal) / run.numb

paste("Prop. of sims with larger/equal S>C than observed:", pval.SC)
paste("Prop. of sims with larger/equal C>S than observed:", pval.CS)
paste("Prop. of sims with larger/equal Equal than observed:", pval.Equal)

paste(sum(simtab$ratio.SC.to.CS>1) / run.numb*100, "% of sims with more S>C than C>S", sep="")
paste(sum(simtab$ratio.SC.to.CS<1) / run.numb*100, "% of sims with more C>S than S>C", sep="")
paste(sum(simtab$ratio.SC.to.CS==1) / run.numb*100, "% of sims with Equal numbers of S>C and C>S", sep="")

#plot expected 
pdf(file = output_name_pdf4,width = 6, height=10)
par(mfrow=c(3,1),font.axis=2,font.lab=2,lwd=1.2,cex=1,mar=c(4,4,1,1))
plot(sort(simtab$propSC),ylim=c(0,1), xlab = "", ylab = "expected Prop. S>C", pch = 21, bg = "grey")
abline(h = observed.prop.SC,col = "red")
abline(h = mean(simtab$propSC),col = "black",lty=3)
text(run.numb / 5,0.95,paste("Pval of S>C:", pval.SC), col="red")
plot(sort(simtab$propCS),ylim=c(0,1), xlab = "", ylab = "expected Prop. C>S", pch = 21, bg = "grey")
abline(h = observed.prop.CS,col = "red")
abline(h = mean(simtab$propCS),col = "black",lty=3)
text(run.numb / 5,0.95,paste("Pval of C>S:", pval.CS), col="red")
plot(sort(simtab$propEqual),ylim=c(0,1), xlab = "Simulations", ylab = "expected Prop. Equal", pch = 21, bg = "grey")
abline(h = observed.prop.Equal,col = "red")
abline(h = mean(simtab$propEqual),col = "black",lty=3)
text(run.numb / 5,0.95,paste("Pval of Equal:", pval.Equal), col="red")
dev.off()

#plot RATIOS
above0 = sum(log2(simtab$ratio.SC.to.CS) > 0) / run.numb #proportion of sims with higher number of S>C
below0 = sum(log2(simtab$ratio.SC.to.CS) < 0) / run.numb #proportion of sims with higher number of S>C
equal0 = sum(log2(simtab$ratio.SC.to.CS) == 0) / run.numb

minY = min(log2(simtab$ratio.SC.to.CS))
maxY = max(log2(simtab$ratio.SC.to.CS))

pdf(file = output_name_pdf5,width = 6, height=8)
par(mfrow=c(2,1),font.axis=2,font.lab=2,lwd=1.2,cex=1,mar=c(4,4,1,1))
plot(sort(log2(simtab$ratio.SC.to.CS)), ylim=c(floor(minY),ceiling(maxY)), xlab = "", ylab = "log2 Ratio S>C to C>S", pch = 21, bg = "grey")
#abline(h = log(observed.prop.SC / observed.prop.CS),col = "red")
abline(h =0)
abline(h = log2(mean(simtab$ratio.SC.to.CS)),col = "black",lty=3)
text(run.numb / 4, ceiling(maxY),paste("Prop of Sims with more S>C:", above0), col="red")
text(run.numb / 4, ceiling(maxY)*0.9,paste("Prop of Sims with more C>S:", below0), col="red")
text(run.numb / 4, ceiling(maxY)*0.75,paste("Prop of Sims with same numb of S>C and C>S:", equal0), col="red")
text(run.numb / 4, ceiling(maxY)*0.6,paste("Observed log2 ratio of S>C to C>S:", round(log2(observed.prop.SC / observed.prop.CS),2)), col="red")

above0.equ = sum(log(simtab$ratio.Equal.to.SC.CS) > 0) / run.numb #proportion of sims with higher number of S>C
below0.equ = sum(log(simtab$ratio.Equal.to.SC.CS) < 0) / run.numb #proportion of sims with higher number of S>C
equal0.equ = sum(log(simtab$ratio.Equal.to.SC.CS) == 0) / run.numb #proportion of sims with higher number of S>C


minY = min(log2(simtab$ratio.Equal.to.SC.CS))
maxY = max(log2(simtab$ratio.Equal.to.SC.CS))

plot(sort(log2(simtab$ratio.Equal.to.SC.CS)), ylim=c(floor(minY),ceiling(maxY)), xlab = "", ylab = "log2 Ratio Equal to S>C+C>S", pch = 21, bg = "grey")
abline(h = log2(observed.prop.Equal / (observed.prop.SC+observed.prop.CS)),col = "red")
abline(h =0)
abline(h = log(mean(simtab$ratio.Equal.to.SC.CS)),col = "black",lty=3)
text(run.numb / 3, ceiling(maxY),paste("Prop of Sims with more Equal than S>C+C>S:", above0.equ), col="red")
text(run.numb / 3, (ceiling(maxY)*1.1),paste("Prop of Sims with less Equal than S>C+C>S:", below0.equ), col="red")
text(run.numb / 3, (ceiling(maxY)*1.2),paste("Prop of Sims with Equal same as S>C+C>S:", equal0.equ), col="red")
text(run.numb / 3, ceiling(maxY)*1.4,
     paste("Observed log2 ratio of Equal to S>C+C>S:", round(log2(observed.prop.Equal / (observed.prop.SC+observed.prop.CS)),2)), 
     col="red")
dev.off()

#PLOT consistent differences
pdf(file = output_name_pdf6,width = 8, height=6)
par(mfrow=c(1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1,mar=c(4,4,1,1))
plot(simtab$SC.consist, ylim=c(0,length(af)), xlab = "", ylab = "Consistent across Replicates", pch = 21, bg = "red")
points(simtab$CS.consist,pch = 21, bg = "blue")
abline(h=consistent.SC.CS.TEs[1],col="red",lty=2,lwd=2)
abline(h=consistent.SC.CS.TEs[2],col="blue",lty=2,lwd=2)
dev.off()

paste("STARTING SIMS FOR EACH TE")


#Simulate prob ability of single TEs
fabian_prob_af = NULL
for (j in 1:length(af)){
  #print(af[j])
  fabian.sim = NULL
  for (i in 1:run.numb){
    fabian.sel1_2 = replicate(2,drifTE(N = popsizeS, pi = af[j], gen = 144, R = 1))
    fabian.sel3 = replicate(1,drifTE(N = popsizeS, pi = af[j], gen = 145, R = 1))
    fabian.sel4 = replicate(1,drifTE(N = popsizeS, pi = af[j], gen = 148, R = 1))
    
    fabian.cont1 = replicate(1,drifTE(N = popsizeC, pi = af[j], gen = 293, R = 1))
    fabian.cont2 = replicate(1,drifTE(N = popsizeC, pi = af[j], gen = 310, R = 1))
    
    fabian.sel1_2 = lapply(seq(dim(fabian.sel1_2)[3]), function(x) fabian.sel1_2[ , , x])
    fabian.sel3 = lapply(seq(dim(fabian.sel3)[3]), function(x) fabian.sel3[ , , x])
    fabian.sel4 = lapply(seq(dim(fabian.sel4)[3]), function(x) fabian.sel4[ , , x])
    
    fabian.cont1 = lapply(seq(dim(fabian.cont1)[3]), function(x) fabian.cont1[ , , x])
    fabian.cont2 = lapply(seq(dim(fabian.cont2)[3]), function(x) fabian.cont2[ , , x])
    
    fabian.sel1_2.unlist = unlist(lapply(fabian.sel1_2, function(x) x[length(x)]))
    fabian.sel3.unlist = unlist(lapply(fabian.sel3, function(x) x[length(x)]))
    fabian.sel4.unlist = unlist(lapply(fabian.sel4, function(x) x[length(x)]))
    
    fabian.cont1.unlist = unlist(lapply(fabian.cont1, function(x) x[length(x)]))
    fabian.cont2.unlist = unlist(lapply(fabian.cont2, function(x) x[length(x)]))
    
    fabian.sel.unlist = c(fabian.sel1_2.unlist,fabian.sel3.unlist,fabian.sel4.unlist)
    fabian.cont.unlist = c(fabian.cont1.unlist , fabian.cont2.unlist)
    
    fabian.sim.int = cbind(MeanSel = mean(fabian.sel.unlist),
                           MeanCont = mean(fabian.cont.unlist),
                           dFreq = mean(fabian.sel.unlist) - mean(fabian.cont.unlist)) #DO NOT ROUND
    fabian.sim = rbind(fabian.sim,fabian.sim.int)
  }
  fabian.sim = as.data.frame(fabian.sim)
  fabian_prob_af = rbind(fabian_prob_af,
                       cbind(start = af[j],
                             MeanSel = mean(fabian.sim$MeanSel),
                             MeanCont = mean(fabian.sim$MeanCont),
                             dFreq_mean = mean(fabian.sim$dFreq),
                             dFreq_med = median(fabian.sim$dFreq),
                             totalSmore = sum(fabian.sim$dFreq>0), 
                             chanceSmore = sum(fabian.sim$dFreq>0) / run.numb,
                             chanceCmore = sum(fabian.sim$dFreq<0) / run.numb,
                             chanceEqual = sum(fabian.sim$dFreq==0) / run.numb) )
  #round to 2 decimal, otherwise "too" exact and S>C, C>S grouping different - if not rounded, hardly any TEs are in group 'Same'
}
fabian_prob_af = as.data.frame(fabian_prob_af)
fabian_prob_af$chanceSmore = fabian_prob_af$chanceSmore
fabian_prob_af$chanceCmore = fabian_prob_af$chanceCmore
fabian_prob_af$chanceEqual = fabian_prob_af$chanceEqual
#fabian_prob_af

#Probability of how often a TE is more in Sel than C
pdf(file = output_name_pdf.All,width = 8, height=14)
par(mfrow=c(3,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4,mar=c(3,4,1,1))
plot(fabian_prob_af$start,fabian_prob_af$chanceSmore, 
     xlab = "", ylab = "P more in S", 
     main = paste("Runs = ",run.numb),
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = 0.5, lty =2, col = "red",lwd=2)
abline(h = mean(fabian_prob_af$chanceSmore), lty = 3, lwd=3, col ="black")
points(fabian_prob_af$start,fabian_prob_af$chanceSmore,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(fabian_prob_af$chanceSmore>0.5) / length(fabian_prob_af$chanceSmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(fabian_prob_af$chanceSmore<=0.5)) / length(fabian_prob_af$chanceSmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(fabian_prob_af$chanceSmore==0.5)+2) / length(fabian_prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(fabian_prob_af$start),"TEs from SA pop"),cex=0.7)
# dev.off()

# pdf(file = output_name_pdf.C,width = 8, height=8)
# par(mfrow=c(1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
plot(fabian_prob_af$start,fabian_prob_af$chanceCmore, 
     xlab = "", ylab = "P more in C", 
     main = "",
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = 0.5, lty =2, col = "red",lwd=2)
abline(h = mean(fabian_prob_af$chanceCmore), lty = 3, lwd=3, col ="black")
points(fabian_prob_af$start,fabian_prob_af$chanceCmore,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(fabian_prob_af$chanceCmore>0.5) / length(fabian_prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(fabian_prob_af$chanceCmore<=0.5)) / length(fabian_prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(fabian_prob_af$chanceSmore==0.5)+2) / length(fabian_prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(fabian_prob_af$start),"TEs from SA pop"),cex=0.7)
# dev.off()

# pdf(file = output_name_pdf.C,width = 8, height=8)
# par(mfrow=c(1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
plot(fabian_prob_af$start,fabian_prob_af$chanceEqual, 
     xlab = "", ylab = "P of Equal", 
     main = "",
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = 0.5, lty =2, col = "red",lwd=2)
abline(h = mean(fabian_prob_af$chanceEqual), lty = 3, lwd=3, col ="black")
points(fabian_prob_af$start,fabian_prob_af$chanceEqual,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(fabian_prob_af$chanceEqual>0.5) / length(fabian_prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(fabian_prob_af$chanceEqual<=0.5)) / length(fabian_prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(fabian_prob_af$chanceSmore==0.5)+2) / length(fabian_prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(fabian_prob_af$start),"TEs from SA pop"),cex=0.7)
mtext(line = 2,side = 1, text = "Starting Frequency (SA Pop)",cex = 1.4, font=2)
dev.off()


#Probability for each case boxplot
fabian_prob_af.chance.melt = melt(fabian_prob_af[,7:9])
fabian_prob_af.chance.melt$variable = mgsub(c("chance","Sm","Cm"),c("", "S m", "C m"),fabian_prob_af.chance.melt$variable)
fabian_prob_af.chance.melt$variable = factor(fabian_prob_af.chance.melt$variable, c( 'S more', 'C more','Equal'))

p = ggplot(data = fabian_prob_af.chance.melt, aes(x=variable, y=value, group = variable, color = variable)) + 
  geom_boxplot(color="black", outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = expression(bold("Probability"))) + 
  geom_jitter(width = 0.1, color='black', pch = 21)

ggsave(p, file=output_name_pdf2, height=7, width=6)

#Select direction and plot
Direction = NULL
for (i in 1:length(fabian_prob_af$chanceSmore)){
  valueS = fabian_prob_af$chanceSmore[i]
  valueC = fabian_prob_af$chanceCmore[i]
  valueE = fabian_prob_af$chanceEqual[i]
  
  if (valueS > dir.cutoff) {
    Direction  = c(Direction,"S>C")
  } else if (valueC > dir.cutoff) {
    Direction  = c(Direction,"C>S") 
  } else if (valueC <= dir.cutoff | valueS <= dir.cutoff) {
    Direction  = c(Direction,"Same")
  }
}
fabian_prob_af$Direction = Direction

#Write tab
write.table(fabian_prob_af, output_name,quote=F,sep="\t",row.names = F)

#Starting frequency of groups
plot.col = mgsub(c("C>S","S>C","Same"),c("blue","red","grey"),fabian_prob_af$Direction)
pdf(output_name_pdf3, width=7, height=8)
par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(fabian_prob_af$start ~ as.factor(fabian_prob_af$Direction), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(fabian_prob_af$start ~ as.factor(fabian_prob_af$Direction), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C","Same"),c(1,2,3),fabian_prob_af$Direction)),factor=0.5), 
       fabian_prob_af$start, 
       pch = 21, 
       bg= plot.col,cex=1.4)
axis(1,seq(1,3,1),labels = c("C>S","S>C","Same"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)
# axis(4,seq(0,1,0.2),labels = F, line = 0,mgp=c(3, .5, 0))
title(ylab = "Starting Frequency (SA Pop)",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
pval = round(t.test(fabian_prob_af[fabian_prob_af$Direction == "S>C",]$start,
                    fabian_prob_af[fabian_prob_af$Direction == "C>S",]$start)$p.value,2)

text(1.5,1.14,paste("P=",pval),cex=1.6)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()


fabian_prob_af
print("RUN ENDED")



