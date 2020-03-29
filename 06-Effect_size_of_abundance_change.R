#Magnitude of abundance changed in log2FC (i.e. log2 of insertions in selected divided by controls) and absolute values (insertion difference between selection and controls, i.e. delta insertions)

#For Figure 1BC and Figure S5

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(VennDiagram)
library(gridExtra)
library(ggplot2)
library(SuperExactTest)
library(FactoMineR)
library("Hmisc")
library(rstatix) 
library(patchwork)

output_path = "/Users/danfab/comp_all/"

#Plot effect size of TEs that have more/less insertions
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt.UP = carnes.stat.filt[carnes.stat.filt$Diff_SelCont>0,]
carnes.stat.filt.DOWN = carnes.stat.filt[carnes.stat.filt$Diff_SelCont<0,]

carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE",]
carnes.stat.filt.signUP = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont>0,]
carnes.stat.filt.signDOWN = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont<0,]

fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt.UP = fabian.stat.filt[fabian.stat.filt$Diff_SelCont>0,]
fabian.stat.filt.DOWN = fabian.stat.filt[fabian.stat.filt$Diff_SelCont<0,]

fabian.stat.filt.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE",]
fabian.stat.filt.signUP = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont>0,]
fabian.stat.filt.signDOWN = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont<0,]

remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
remo.stat.filt.UP = remo.stat.filt[remo.stat.filt$Diff_SelCont>0,]
remo.stat.filt.DOWN = remo.stat.filt[remo.stat.filt$Diff_SelCont<0,]

remo.stat.filt.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE",]
remo.stat.filt.signUP = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont>0,]
remo.stat.filt.signDOWN = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont<0,]

hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == TRUE,]

hoed.stat.filt.UP = hoed.stat.filt[hoed.stat.filt$Diff_PostEarly>0,]
hoed.stat.filt.DOWN = hoed.stat.filt[hoed.stat.filt$Diff_PostEarly<0,]

hoed.stat.filt.signUP = hoed.stat.filt.sign[hoed.stat.filt.sign$Diff_PostEarly>0,]
hoed.stat.filt.signDOWN = hoed.stat.filt.sign[hoed.stat.filt.sign$Diff_PostEarly<0,]

#Group TE families
carn.effect.size = rbind(cbind("Carnes","S>C",carnes.stat.filt.signUP$Diff_SelCont),
                         cbind("Carnes","C>S",abs(carnes.stat.filt.signDOWN$Diff_SelCont)))
fab.effect.size = rbind(cbind("Fabian","S>C",fabian.stat.filt.signUP$Diff_SelCont),
                        cbind("Fabian","C>S",abs(fabian.stat.filt.signDOWN$Diff_SelCont)))
hoed.effect.size = rbind(cbind("Hoedjes","S>C",hoed.stat.filt.signUP$Diff_PostEarly),
                         cbind("Hoedjes","C>S",abs(hoed.stat.filt.signDOWN$Diff_PostEarly)))
remo.effect.size = rbind(cbind("Remolina","S>C",remo.stat.filt.signUP$Diff_SelCont),
                         cbind("Remolina","C>S",abs(remo.stat.filt.signDOWN$Diff_SelCont)))

combi.effect.size = as.data.frame(rbind(carn.effect.size, fab.effect.size, hoed.effect.size, remo.effect.size))
names(combi.effect.size) = c("Study","Direction","Diff")
combi.effect.size$Diff = numfact(combi.effect.size$Diff)
sum.effect.size = tapply(combi.effect.size$Diff,list(combi.effect.size$Study,combi.effect.size$Direction),sum)

################################################################################
#T-tests on S>C vs C>S effect size (delta insertions)
################################################################################

t.test(Diff ~ Direction, data= combi.effect.size[combi.effect.size$Study == "Carnes",])
#t = 2.1833, df = 21.249, p-value = 0.04037
t.test(Diff ~ Direction, data= combi.effect.size[combi.effect.size$Study == "Fabian",])
#t = 2.9505, df = 47.605, p-value = 0.00491
t.test(Diff ~ Direction, data= combi.effect.size[combi.effect.size$Study == "Remolina",])
#t = -1.6617, df = 36.15, p-value = 0.1052
t.test(Diff ~ Direction, data= combi.effect.size[combi.effect.size$Study == "Hoedjes",])
#t = -0.30196, df = 33.397, p-value = 0.7646

#Plot
p.eff = ggplot(data = combi.effect.size, aes(x=paste(Study,Direction), y=Diff,color=Study)) + geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = expression(bold(paste(delta, "Insertions")))) + 
  scale_color_manual(values = c("steelblue1","blue","darkgreen","red")) + 
  scale_x_discrete(labels=c("C>S", "S>C", "C>S","S>C", "C>S","S>C","C>S","S>C")) + 
  geom_jitter(width=0.1) + theme(legend.position="none")
my_comp1 = list( c("Carnes C>S","Carnes S>C"),c("Fabian C>S","Fabian S>C"),c("Hoedjes C>S","Hoedjes S>C"),c("Remolina C>S","Remolina S>C"))
p.eff = p.eff + stat_compare_means(comparisons = my_comp1, 
                                   label.y = 56, 
                                   method = "t.test", 
                                   tip.length = 0.01, bracket.size = 0.6, 
                                   symnum.args = symnum.args,
                                   size = 8) + ylim(c(0,62))
combi.effect.size$Study = mgsub(c('Carnes' ,'Fabian' ,'Hoedjes' ,'Remolina'),c('Carnes2015' ,'Fabian2018' ,'Hoedjes2019' ,'Remolina2012'),combi.effect.size$Study)
p.eff.study = ggplot(data = combi.effect.size, aes(x=Study, y=Diff,color=Study)) + geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = "") + 
  scale_color_manual(values = c("steelblue1","blue","darkgreen","red")) + 
  scale_x_discrete(labels=c("Carnes", "Fabian",'Hoedjes', "Remolina")) + 
  geom_jitter(width=0.1) + theme(axis.text.y = element_blank())

stat.test = aov(Diff ~ Study, data = combi.effect.size) %>% tukey_hsd() #requires library(rstatix) 
stat.test
summary(aov(Diff ~ Study, data = combi.effect.size) )
#stat.test = stat.test[c(1,3,2),] #reaarange rows
#stat.test$Study = c("Carnes","Fabian",'Hoedjes',"Remolina")
stat.test$p.adj.lab = c("***","***","***",'ns','ns','ns')
p.eff.study = p.eff.study + 
  stat_pvalue_manual(as.data.frame(stat.test), label = "p.adj.lab",
                     y.position = c(54,60,63,52,56,50),
                     tip.length = 0.01, bracket.size = 0.8, 
                     symnum.args = symnum.args,
                     size = 8, color = "black") + ylim(c(0,63))
p.eff.study
p.eff + p.eff.study
ggsave(p.eff + p.eff.study, file=paste(output_path,"FigS5_insertion_effect_size_Regime_Study.pdf",sep=""), height=6, width=14)


#Same for log2 FC
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter_withConsistent.txt",header = T)
carnes.stat.filt.UP = carnes.stat.filt[carnes.stat.filt$Diff_SelCont>0,]
carnes.stat.filt.DOWN = carnes.stat.filt[carnes.stat.filt$Diff_SelCont<0,]

carnes.stat.filt.sign = carnes.stat.filt[carnes.stat.filt$Bonf == "TRUE",]
carnes.stat.filt.signUP = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont>0,]
carnes.stat.filt.signDOWN = carnes.stat.filt.sign[carnes.stat.filt.sign$Diff_SelCont<0,]

fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter_withConsistent.txt",header = T)
fabian.stat.filt.UP = fabian.stat.filt[fabian.stat.filt$Diff_SelCont>0,]
fabian.stat.filt.DOWN = fabian.stat.filt[fabian.stat.filt$Diff_SelCont<0,]

fabian.stat.filt.sign = fabian.stat.filt[fabian.stat.filt$Bonf == "TRUE",]
fabian.stat.filt.signUP = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont>0,]
fabian.stat.filt.signDOWN = fabian.stat.filt.sign[fabian.stat.filt.sign$Diff_SelCont<0,]

remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter_withConsistent.txt",header = T)
remo.stat.filt.UP = remo.stat.filt[remo.stat.filt$Diff_SelCont>0,]
remo.stat.filt.DOWN = remo.stat.filt[remo.stat.filt$Diff_SelCont<0,]

remo.stat.filt.sign = remo.stat.filt[remo.stat.filt$Bonf == "TRUE",]
remo.stat.filt.signUP = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont>0,]
remo.stat.filt.signDOWN = remo.stat.filt.sign[remo.stat.filt.sign$Diff_SelCont<0,]

hoed.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter_withConsistent.txt",header = T)
hoed.stat.filt.sign = hoed.stat.filt[hoed.stat.filt$Bonf_Full == TRUE,]

hoed.stat.filt.UP = hoed.stat.filt[hoed.stat.filt$Diff_PostEarly>0,]
hoed.stat.filt.DOWN = hoed.stat.filt[hoed.stat.filt$Diff_PostEarly<0,]

hoed.stat.filt.signUP = hoed.stat.filt.sign[hoed.stat.filt.sign$Diff_PostEarly>0,]
hoed.stat.filt.signDOWN = hoed.stat.filt.sign[hoed.stat.filt.sign$Diff_PostEarly<0,]

#Group
carn.effect.size = rbind(cbind("Carnes","S>C",carnes.stat.filt.signUP$log2RatSelCont),
                         cbind("Carnes","C>S",abs(carnes.stat.filt.signDOWN$log2RatSelCont)))

fab.effect.size = rbind(cbind("Fabian","S>C",fabian.stat.filt.signUP$log2RatSelCont),
                        cbind("Fabian","C>S",abs(fabian.stat.filt.signDOWN$log2RatSelCont)))

hoed.effect.size = rbind(cbind("Hoedjes","S>C",hoed.stat.filt.signUP$log2Rat_PostEarly),
                         cbind("Hoedjes","C>S",abs(hoed.stat.filt.signDOWN$log2Rat_PostEarly)))

remo.effect.size = rbind(cbind("Remolina","S>C",remo.stat.filt.signUP$log2RatSelCont),
                         cbind("Remolina","C>S",abs(remo.stat.filt.signDOWN$log2RatSelCont)))

combi.effect.size = as.data.frame(rbind(carn.effect.size, fab.effect.size, hoed.effect.size,remo.effect.size))
names(combi.effect.size) = c("Study","Direction","Diff")
combi.effect.size$Diff = numfact(combi.effect.size$Diff)
sum.effect.size = tapply(combi.effect.size$Diff,list(combi.effect.size$Study,combi.effect.size$Direction),sum)

################################################################################
#T-tests on S>C vs C>S effect size (log2 fold change insertions)
################################################################################

t.test(Diff ~ Direction, data = combi.effect.size[combi.effect.size$Study == "Carnes",])
#t = -0.739, df = 28.818, p-value = 0.4659
t.test(Diff ~ Direction, data = combi.effect.size[combi.effect.size$Study == "Fabian",])
#t = 2.6383, df = 42.995, p-value = 0.01156
t.test(Diff ~ Direction, data = combi.effect.size[combi.effect.size$Study == "Remolina",])
#t = -2.9589, df = 38.231, p-value = 0.005275
t.test(Diff ~ Direction, data = combi.effect.size[combi.effect.size$Study == "Hoedjes",])
#t = -2.8275, df = 90.718, p-value = 0.005772

#Plot
p.eff = ggplot(data = combi.effect.size, aes(x=paste(Study,Direction), y=Diff,color=Study)) + geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = expression(bold(paste("log"[2], " Insertions"))) ) + 
  scale_color_manual(values = c("steelblue1","blue",'darkgreen',"red")) + 
  scale_x_discrete(labels=c("C>S", "S>C", "C>S","S>C", "C>S","S>C","C>S","S>C")) + 
  geom_jitter(width=0.1) + theme(legend.position="none")
my_comp1 = list( c("Carnes C>S","Carnes S>C"),c("Fabian C>S","Fabian S>C"),c("Hoedjes C>S","Hoedjes S>C"),c("Remolina C>S","Remolina S>C"))
p.eff = p.eff + stat_compare_means(comparisons = my_comp1, 
                                   label.y = 1.1, 
                                   method = "t.test", 
                                   tip.length = 0.01, bracket.size = 0.6, 
                                   symnum.args = symnum.args,
                                   size = 8) + ylim(c(0,1.2))
combi.effect.size$Study = mgsub(c('Carnes' ,'Fabian' ,'Hoedjes' ,'Remolina'),c('Carnes2015' ,'Fabian2018' ,'Hoedjes2019' ,'Remolina2012'),combi.effect.size$Study)
p.eff.study = ggplot(data = combi.effect.size, aes(x=Study, y=Diff,color=Study)) + geom_boxplot(outlier.shape = NA) + theme_bw() + THEME_DEF + 
  labs(x = expression(bold("")), y = "") + 
  scale_color_manual(values = c("steelblue1","blue",'darkgreen',"red")) + 
  scale_x_discrete(labels=c("Carnes", "Fabian",'Hoedjes' ,"Remolina")) + 
  geom_jitter(width=0.1) + theme(axis.text.y = element_blank())

stat.test = aov(Diff ~ Study, data = combi.effect.size) %>% tukey_hsd() #requires library(rstatix) 
stat.test
#stat.test = stat.test[c(1,3,2),] #reaarange rows
stat.test$Study = c("Carnes","Fabian",'Hoedjes',"Remolina")
stat.test$p.adj.lab = c("***","***","***","***","***",'ns')
p.eff.study = p.eff.study + 
  stat_pvalue_manual(as.data.frame(stat.test), label = "p.adj.lab",
                     y.position = c(1.04,1.15,1.20, 1.01, 1.08, 0.98),
                     tip.length = 0.01, bracket.size = 0.8, 
                     symnum.args = symnum.args,
                     size = 8, color = "black") + ylim(c(0,1.2))
p.eff.study
p.eff + p.eff.study
ggsave(p.eff + p.eff.study, file=paste(output_path,"Fig1_insertion_effect_size_LOG2_Regime_Study.pdf",sep=""), height=6, width=14)
