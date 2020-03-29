#Simulate TE frequency change over generations in selected and control populations

#Please change folder names and paths.
dir.create("/Users/danfab/effective_pop_size/")
setwd("/Users/danfab/effective_pop_size/")

#
library(pBrackets)
library(reshape2)
library(ggplot2)
library(Rfast)

#Datasets
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
remo.consistent.SC.CS.TEs = c(3,0)
carnes.consistent.SC.CS.TEs = c(48,2)
hoedjes.consistent.SC.CS.TEs = c(16.7,11.7)

#Sign tab for Hoedjes
#hoedjes.stat.filt.sign = read.table("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/hoedjes_all/output_stat/FullModel_Breeding_SignTEs_Bonf001.txt",header = T)
bonf.hoedjes = 0.01 / nrow(hoedjes.stat.filt)
hoedjes.stat.filt$Bonf = hoedjes.stat.filt$PvalFull_Breed < bonf.hoedjes

#Robert PlosGen 2015 average frequencies
robert.2015 = read.table("/Users/danfab/extra_files/Robert_Plos2015_AvFreq.txt",header=T)
names(robert.2015)[1] = "Name"
robert.2015$av_freq = numfact(robert.2015$av_freq)
robert.2015 = na.omit(robert.2015)
robert.2015$freq_rank = rank(robert.2015$av_freq)

#####################
#CUSTOM FUNCTIONS + PLOT SETTINGS
#####################

#Number of Runs (replicates) for experimetnal evolution studies
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

#Multiple gsub
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
    pi = vapply(pi, function(x) rbinom(1, N*2, prob = x), as.integer(1L)) / (N*2) #SET DIPLOID N * 2
    df=cbind(df,pi)
    #popsize=rbind(popsize,N)
    #gVar = q*(1-q) / N^R
    #gVarDF=rbind(gVarDF,gVar)
    N = round(N^R)
  }
  return(df)
}


#set R = 1
R = 1 #this is not used

#####################
#args to pass
#####################

args = commandArgs(trailingOnly=TRUE)
popsizeC = as.integer(args[1]) #e.g. 640 NEED TO TAKE x2 BECAUSE DIPLOID
popsizeS = as.integer(args[2]) #e.g. 320 NEED TO TAKE x2 BECAUSE DIPLOID
genC = as.integer(args[3]) #e.g. 80
genS = as.integer(args[4]) #e.g. 50
reps = as.integer(args[5]) #e.g. 3
dir.cutoff = as.character(args[6]) ##needs to be >=0.5 - i.e. if a TE has a chance of >0.5 to be S>C, it is counted as S>C
dir.cutoff = as.numeric(dir.cutoff)
run.numb = as.integer(args[7]) #e.g. 100
study = as.character(args[8])
setfreq = as.character(args[9])
ins = as.integer(args[10])

setfreq = as.numeric(setfreq)
setfreq

#To try
# genS = 100
# genC = 300
# popsizeC = 100000
# popsizeS = 100000
# af = c(0.001,0.01,0.05,0.1)
# af = c(0.25,0.5,0.75)
# af = c(0.9,0.95,0.99,0.999)
# reps = 3 
# ins = 10

#####################
#set output names
#####################

output_name = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,".txt",sep="")
output_name2 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_PropType.txt",sep="")
output_name_pdf.All = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_All.pdf",sep="")
output_name_pdf.Examp = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Example.pdf",sep="")
output_name_pdf2 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Boxplot.pdf",sep="")
output_name_pdf3 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Direction.pdf",sep="")
output_name_pdf4 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_PropType.pdf",sep="")
output_name_pdf5 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_Ratios.pdf",sep="")
output_name_pdf6 = paste("Ncont",popsizeC,"_Nsel",popsizeS,"_C",genC,"_S",genS,"_Reps",reps,"_DirCut",dir.cutoff,"_Runs",run.numb,"_ConsReps.pdf",sep="")

#####################
#RUN SCRIPT
#####################

#Only select TEs that were also detected in each study & proportion of S>C, C>S, Equal based on sharing with SA pop
if (study == "remolina") {
  int = robert.2015[robert.2015$Name %in% gsub('Dmel/','',remo.stat.filt$flybase_name),]
  af = sort(int$av_freq)
  int2 = remo.stat.filt[gsub('Dmel/','',remo.stat.filt$flybase_name) %in% robert.2015$Name,]
  int3 = int2[int2$Bonf == "TRUE",]
  observed.prop.SC = sum(int3$Mean_Sel > int3$Mean_Cont) / nrow(int2)
  observed.prop.CS = sum(int3$Mean_Sel < int3$Mean_Cont) / nrow(int2)
  observed.prop.Equal = (nrow(int2) - nrow(int3)) / nrow(int2)
  consistent.SC.CS.TEs = remo.consistent.SC.CS.TEs
} else if (study == "carnes") {
  int = robert.2015[robert.2015$Name %in% gsub('Dmel/','',carnes.stat.filt$flybase_name),]
  af = sort(int$av_freq)
  int2 = carnes.stat.filt[gsub('Dmel/','',carnes.stat.filt$flybase_name) %in% robert.2015$Name,]
  int3 = int2[int2$Bonf == "TRUE",]
  observed.prop.SC = sum(int3$Mean_Sel > int3$Mean_Cont) / nrow(int2)
  observed.prop.CS = sum(int3$Mean_Sel < int3$Mean_Cont) / nrow(int2)
  observed.prop.Equal = (nrow(int2) - nrow(int3)) / nrow(int2)
  consistent.SC.CS.TEs = carnes.consistent.SC.CS.TEs
} else if (study == "hoedjes") {
  int = robert.2015[robert.2015$Name %in% gsub('Dmel/','',hoedjes.stat.filt$flybase_name),]
  af = sort(int$av_freq)
  int2 = hoedjes.stat.filt[gsub('Dmel/','',hoedjes.stat.filt$flybase_name) %in% robert.2015$Name,]
  int3 = int2[int2$Bonf == "TRUE",]
  observed.prop.SC = sum(int3$Mean_Posponed > int3$Mean_Early) / nrow(int2)
  observed.prop.CS = sum(int3$Mean_Posponed < int3$Mean_Early) / nrow(int2)
  observed.prop.Equal = (nrow(int2) - nrow(int3)) / nrow(int2)
  consistent.SC.CS.TEs = hoedjes.consistent.SC.CS.TEs
} else if (study == "all") {
  af = sort(robert.2015$av_freq)
} else if (study == "choose") {
  af = rep(setfreq,ins)
}
#af
length(af)

#Show example of genetic drift through generations per TE
TEsim_S = replicate(reps,drifTE(popsizeS, af, genS, R))
TEsim_C = replicate(reps,drifTE(popsizeC, af, genC, R))

TEsim_S1 = lapply(seq(dim(TEsim_S)[3]), function(x) TEsim_S[ , , x])
TEsim_C1 = lapply(seq(dim(TEsim_C)[3]), function(x) TEsim_C[ , , x])

pdf(file = output_name_pdf.Examp,width = reps * 3.4, height=6)
par(mfrow=c(2,reps),font=2,font.axis=2,font.lab=2,font.main=2,mar=c(4,4,2,1))
for (j in 1:reps){
  plot(0:genS,TEsim_S[1,,j],
       type="l",ylim=c(0,1),col=1,xlab="",
       ylab="",
       main=paste("Selected N=",popsizeS),xaxt="t")
  mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
  grid(lty=1,col="grey95")
  #axis(1,seq(0,genS,5))
  for (i in 2:length(af)){
    lines(0:genS,TEsim_S[i,,j],type="l",col=i)
  }
}

for (j in 1:reps){
  plot(0:genC,TEsim_C[1,,j],
       type="l",ylim=c(0,1),col=1,xlab="",
       ylab="",
       main=paste("Control N=",popsizeC),xaxt="t")
  grid(lty=1,col="grey95")
  #axis(1,seq(0,genC,5))
  mtext(line = 2.2,side = 2, text = paste("Frequ. of",length(af),"Alleles",sep=" "),cex = 0.8, font=2)
  mtext(line = 2.2,side = 1, text = "Generations",cex = 0.8, font=2)
  for (i in 2:length(af)){
    lines(0:genC,TEsim_C[i,,j],type="l",col=i)
  }
}
dev.off()

#Simulate how many S>C and C>S TEs expected by chance: Proportion of TEs more in S, C, or equal
#Mean across replicates - DON'T ROUND - then check how many TEs are S>C, C>S, or equal
paste("STARTING SIMS FOR PROPORTION OF S>C AND C>S")
simtab = NULL
for (i in 1:run.numb){
  sel = replicate(reps,drifTE(N = popsizeS, pi = af, gen = genS, R = 1))
  cont = replicate(reps,drifTE(N = popsizeC, pi = af, gen = genC, R = 1))
  
  sel = lapply(seq(dim(sel)[3]), function(x) sel[ , , x])
  cont = lapply(seq(dim(cont)[3]), function(x) cont[ , , x])
  
  #Average diff:
  sel.mean = Reduce("+",sel) / reps
  cont.mean = Reduce("+",cont) / reps
  
  SC.lastgen = sel.mean[,ncol(sel.mean)]
  CS.lastgen = cont.mean[,ncol(cont.mean)]
  
  propSC = sum(SC.lastgen > CS.lastgen) / length(af)
  propCS = sum(SC.lastgen < CS.lastgen) / length(af)
  propEqual = sum(SC.lastgen == CS.lastgen) / length(af)
  ratio.SC.to.CS = propSC / propCS #larger 1 = more S>C
  ratio.Equal.to.SC.CS = propEqual / (propSC / propCS) #larger 1 = more Equal
 
  #Consistent diff across reps:
  sel.lastgen = lapply(sel, function(x) x[,ncol(x)])
  cont.lastgen = lapply(cont, function(x) x[,ncol(x)])
  
  sel.lastgen.tab = matrix(unlist(sel.lastgen),ncol=length(af),byrow=T)
  cont.lastgen.tab = matrix(unlist(cont.lastgen),ncol=length(af),byrow=T)
  
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
#paste(sum(simtab$ratio.Equal.to.SC.CS>1) / run.numb*100, "% of sims with more Equal than S>C+C>S", sep="")

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
text(run.numb / 4, ceiling(maxY)*0.8,paste("Prop of Sims with same numb of S>C and C>S:", equal0), col="red")
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
#Simulate probability of being more in S or C for each TE separate
prob_af = NULL
for (j in 1:length(af)){
  #print(af[j])
  sim = NULL
  for (i in 1:run.numb){
    sel = replicate(reps,drifTE(N = popsizeS, pi = af[j], gen = genS, R = 1))
    cont = replicate(reps,drifTE(N = popsizeC, pi = af[j], gen = genC, R = 1))
    
    sel = lapply(seq(dim(sel)[3]), function(x) sel[ , , x])
    cont = lapply(seq(dim(cont)[3]), function(x) cont[ , , x])
    
    sel.unlist = unlist(lapply(sel, function(x) x[length(x)]))
    cont.unlist = unlist(lapply(cont, function(x) x[length(x)]))
    
    sim.int = cbind(MeanSel = mean(sel.unlist),
                         MeanCont = mean(cont.unlist),
                         dFreq = mean(sel.unlist) - mean(cont.unlist) ) #DO NOT ROUND
    sim = rbind(sim,sim.int)
  }
  sim = as.data.frame(sim)
  prob_af = rbind(prob_af,
                       cbind(start = af[j],
                             MeanSel = mean(sim$MeanSel),
                             MeanCont = mean(sim$MeanCont),
                             dFreq_mean = mean(sim$dFreq),
                             dFreq_med = median(sim$dFreq),
                             totalSmore = sum(sim$dFreq>0), 
                             chanceSmore = sum(sim$dFreq>0) / run.numb,
                             chanceCmore = sum(sim$dFreq<0) / run.numb,
                             chanceEqual = sum(sim$dFreq==0) / run.numb) )
  #round to 2 decimal, otherwise "too" exact and S>C, C>S grouping different - if not rounded, hardly any TEs are in group 'Same'
  #DON'T ROUND
}
prob_af = as.data.frame(prob_af)
prob_af$chanceSmore = prob_af$chanceSmore
prob_af$chanceCmore = prob_af$chanceCmore
prob_af$chanceEqual = prob_af$chanceEqual
#prob_af

#Probability of how often a TE is more in Sel than C, calculated using 100 runs
pdf(file = output_name_pdf.All,width = 8, height=14)
par(mfrow=c(3,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4,mar=c(3,4,1,1))
plot(prob_af$start,prob_af$chanceSmore, 
     xlab = "", ylab = "P more in S", 
     main = paste("Runs = ",run.numb),
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = dir.cutoff, lty =2, col = "red",lwd=2)
abline(h = mean(prob_af$chanceSmore), lty = 3, lwd=3, col ="black")
points(prob_af$start,prob_af$chanceSmore,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(prob_af$chanceSmore>0.5) / length(prob_af$chanceSmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(prob_af$chanceSmore<=0.5)) / length(prob_af$chanceSmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(prob_af$chanceSmore==0.5)+2) / length(prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(prob_af$start),"TEs from SA pop"),cex=0.7)
# dev.off()

# pdf(file = output_name_pdf.C,width = 8, height=8)
# par(mfrow=c(1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
plot(prob_af$start,prob_af$chanceCmore, 
     xlab = "", ylab = "P more in C", 
     main = "",
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = dir.cutoff, lty =2, col = "red",lwd=2)
abline(h = mean(prob_af$chanceCmore), lty = 3, lwd=3, col ="black")
points(prob_af$start,prob_af$chanceCmore,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(prob_af$chanceCmore>0.5) / length(prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(prob_af$chanceCmore<=0.5)) / length(prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(prob_af$chanceSmore==0.5)+2) / length(prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(prob_af$start),"TEs from SA pop"),cex=0.7)
# dev.off()

# pdf(file = output_name_pdf.C,width = 8, height=8)
# par(mfrow=c(1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
plot(prob_af$start,prob_af$chanceEqual, 
     xlab = "", ylab = "P of Equal", 
     main = "",
     ylim = c(0,1), xlim = c(0,1))
grid(lty=1,col="grey90")
abline(h = dir.cutoff, lty =2, col = "red",lwd=2)
abline(h = mean(prob_af$chanceEqual), lty = 3, lwd=3, col ="black")
points(prob_af$start,prob_af$chanceEqual,pch = 21,bg="grey",cex = 1.2)
text(x = 0.5, y = 1, labels = paste(">0.5: ",round(sum(prob_af$chanceEqual>0.5) / length(prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
text(x = 0.5, y = 0.92, labels = paste("<=0.5: ",round( (sum(prob_af$chanceEqual<=0.5)) / length(prob_af$chanceCmore),3)*100, "% TEs", sep="" ) )
#text(x = 0.5, y = 0.84, labels = paste("=0.5 or =0 (2 fixed TEs): ",round( (sum(prob_af$chanceSmore==0.5)+2) / length(prob_af$chanceSmore),3)*100, "% TEs equal", sep="" ) )
text(0.15,0, labels=paste("n =",length(prob_af$start),"TEs from SA pop"),cex=0.7)
mtext(line = 2,side = 1, text = "Starting Frequency (SA Pop)",cex = 1.4, font=2)
dev.off()



#Probability for each case boxplot
prob_af.chance.melt = melt(prob_af[,7:9])
prob_af.chance.melt$variable = mgsub(c("chance","Sm","Cm"),c("", "S m", "C m"),prob_af.chance.melt$variable)
prob_af.chance.melt$variable = factor(prob_af.chance.melt$variable, c( 'S more', 'C more','Equal'))

#names(prob_af.chance.melt) = c("Group","Probability")

p = ggplot(data = prob_af.chance.melt, aes(x=variable, y=value, group = variable, color = variable)) + 
  geom_boxplot(color="black", outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = expression(bold("Probability"))) + 
  geom_jitter(width = 0.1, color='black', pch = 21)

ggsave(p, file=output_name_pdf2, height=7, width=6)

#Select direction and plot - dont round
Direction = NULL
for (i in 1:length(prob_af$chanceSmore)){
  valueS = prob_af$chanceSmore[i]
  valueC = prob_af$chanceCmore[i]
  valueE = prob_af$chanceEqual[i]
  
  if (valueS > dir.cutoff) {
    Direction  = c(Direction,"S>C")
  } else if (valueC > dir.cutoff) {
    Direction  = c(Direction,"C>S") 
  } else if (valueC <= dir.cutoff | valueS <= dir.cutoff) {
    Direction  = c(Direction,"Same")
  }
}
prob_af$Direction = Direction

#Write tab
write.table(prob_af, output_name,quote=F,sep="\t",row.names = F)

#Starting frequency of groups
plot.col = mgsub(c("C>S","S>C","Same"),c("blue","red","grey"),prob_af$Direction)
pdf(output_name_pdf3, width=7, height=8)
par(mar=c(4,4,1,1),font.axis=2,font.lab=2,lwd=1.2,cex=1.4)
boxplot(prob_af$start ~ as.factor(prob_af$Direction), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA)
grid(col = "lightgray", lty = 1, lwd = 0.5)
par(new=TRUE)
boxplot(prob_af$start ~ as.factor(prob_af$Direction), outpch = NA,axes = F,
        ylim = c(0,1.12),
        whisklty = 1, staplelwd = NA,col="white")
points(jitter(numfact(mgsub(c("C>S","S>C","Same"),c(1,2,3),prob_af$Direction)),factor=0.5), 
       prob_af$start, 
       pch = 21, 
       bg= plot.col,cex=1.4)
axis(1,seq(1,3,1),labels = c("C>S","S>C","Same"),cex.axis=1.2)
axis(2,seq(0,1,0.2),labels = T,cex.axis=1.2)
# axis(4,seq(0,1,0.2),labels = F, line = 0,mgp=c(3, .5, 0))
title(ylab = "Starting Frequency (SA Pop)",line = 2.6,cex.lab=1.4)
title(xlab = "All", line = 2.2,cex.lab=1.4)
pval = round(t.test(prob_af[prob_af$Direction == "S>C",]$start,
                    prob_af[prob_af$Direction == "C>S",]$start)$p.value,2)

text(1.5,1.14,paste("P=",pval),cex=1.6)
brackets(x1 = 1, y1 = 1.05, 
         x2 = 2, y2 = 1.05, 
         h = 0.05, ticks=NA, curvature = 0.5, type = 4,
         col = 1, lwd = 2, lty = 1, xpd = FALSE)
dev.off()


#prob_af
print("RUN ENDED")



