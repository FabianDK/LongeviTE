source("/Users/danfab/R_functions.R")
setwd("/Users/danfab/Science/post-doc/EBI/Projects/Lifespan_TE/manuscript/submission_PLoS_Gen/revision_Code/output")

#Check different effect size cut-offs
carnes.stat.filt = read.table("/Users/danfab/Carnes/TE_maps/res_files/corrected/edit/output_stat/carnes_TE_stat_covfilter.txt",header = T)
fabian.stat.filt = read.table("/Users/danfab/Fabian/TE_maps/res_files/corrected/edit/output_stat/fabian_TE_stat_covfilter.txt",header = T)
hoedjes.stat.filt = read.table("/Users/danfab/Hoedjes/TE_maps/res_files/corrected/edit/output_stat/hoedjes_TE_stat_covfilter.txt",header = T)
hoedjes.stat.filt.sign = hoedjes.stat.filt[hoedjes.stat.filt$Bonf_Full == "TRUE",]
remo.stat.filt = read.table("/Users/danfab/Remolina/TE_maps/res_files/corrected/edit/output_stat/remolina_TE_stat_covfilter.txt",header = T)

carnes.all = nrow(carnes.stat.filt)
fabian.all = nrow(fabian.stat.filt)
hoedjes.all = nrow(hoedjes.stat.filt)
remo.all = nrow(remo.stat.filt)

#Effect size difference
eff.size.vect = seq(0,2,0.1)

carn.eff.size = NULL
fab.eff.size = NULL
remo.eff.size = NULL
hoed.eff.size = NULL

for (i in 1:length(eff.size.vect)){
  eff.size = eff.size.vect[i] #effect size difference
  print(eff.size)
  
  carn.eff = carnes.stat.filt[abs(carnes.stat.filt$Diff_SelCont) > eff.size & carnes.stat.filt$Bonf == "TRUE",]
  carn.vect = c(eff.size,
                nrow(carn.eff[carn.eff$Diff_SelCont>0,]), 
                nrow(carn.eff[carn.eff$Diff_SelCont<0,]),
                nrow(carn.eff),
                "Carnes2015")
  carn.eff.size = rbind(carn.eff.size,
                        carn.vect)
  
  fab.eff = fabian.stat.filt[abs(fabian.stat.filt$Diff_SelCont) > eff.size & fabian.stat.filt$Bonf == "TRUE",]
  fab.vect = c(eff.size,
               nrow(fab.eff[fab.eff$Diff_SelCont>0,]), 
                nrow(fab.eff[fab.eff$Diff_SelCont<0,]),
                nrow(fab.eff),
               "Fabian2018")
  fab.eff.size = rbind(fab.eff.size,
                        fab.vect)
  
  remo.eff = remo.stat.filt[abs(remo.stat.filt$Diff_SelCont) > eff.size & remo.stat.filt$Bonf == "TRUE",]
  remo.vect = c(eff.size,
                nrow(remo.eff[remo.eff$Diff_SelCont>0,]), 
               nrow(remo.eff[remo.eff$Diff_SelCont<0,]),
               nrow(remo.eff),
               "Remolina2012")
  remo.eff.size = rbind(remo.eff.size,
                        remo.vect)
  
  hoed.eff = hoedjes.stat.filt[abs(hoedjes.stat.filt$Diff_PostEarly) > eff.size & hoedjes.stat.filt$Bonf_Full == "TRUE",]
  hoed.vect = c(eff.size,
                nrow(hoed.eff[hoed.eff$Diff_PostEarly>0,]), 
                nrow(hoed.eff[hoed.eff$Diff_PostEarly<0,]),
                nrow(hoed.eff),
                "Hoedjes2019")
  hoed.eff.size = rbind(hoed.eff.size,
                        hoed.vect)
}
carn.eff.size = as.data.frame(carn.eff.size)
fab.eff.size = as.data.frame(fab.eff.size)
remo.eff.size = as.data.frame(remo.eff.size)
hoed.eff.size = as.data.frame(hoed.eff.size)

names(carn.eff.size) = c("Eff_size", "SC","CS","Total","Study")
names(fab.eff.size) = c("Eff_size", "SC","CS","Total","Study")
names(remo.eff.size) = c("Eff_size", "SC","CS","Total","Study")
names(hoed.eff.size) = c("Eff_size", "SC","CS","Total","Study")

carn.eff.size #always more
fab.eff.size #after 0.7, less
remo.eff.size #always more
hoed.eff.size #always more

carn.eff.size$Tot_all = carnes.all
fab.eff.size$Tot_all = fabian.all
remo.eff.size$Tot_all = remo.all
hoed.eff.size$Tot_all = hoedjes.all

all.eff.size = rbind(carn.eff.size,
                     fab.eff.size,
                     hoed.eff.size,
                     remo.eff.size)
all.eff.size$Tot_sign = numfact(all.eff.size$Total)
all.eff.size$CS = numfact(all.eff.size$CS)
all.eff.size$SC = numfact(all.eff.size$SC)
all.eff.size$CS = -all.eff.size$CS

all.eff.size$SCprop = all.eff.size$SC / all.eff.size$Tot_all
all.eff.size$CSprop = all.eff.size$CS / all.eff.size$Tot_all

#Plot raw numbers
all.eff.size.melt = melt(data = all.eff.size[,c(-4,-6)],id.vars = c("Eff_size","Study"),measure.vars = c("SC","CS"))

labs.bar = c( all.eff.size.melt[all.eff.size.melt$variable == "SC",]$value +3, all.eff.size.melt[all.eff.size.melt$variable == "CS",]$value -6)

barplot.eff = ggplot(data=all.eff.size.melt , aes(x=Eff_size, y=value, fill=variable)) +
  geom_bar(stat="identity", color = "black") + THEME_DEF +
  scale_y_continuous(limits=c(-40,90),breaks=seq(-40,90,20), labels= abs(seq(-40,90,20))) + 
  scale_fill_manual(values=c("red","blue")) + geom_hline(yintercept = 0) +
  theme(axis.text.y = element_text(size=18, colour = "black"),
        axis.text.x = element_text(size=18 ,angle = 45,hjust=1),
        axis.title.x = element_text(size=24, vjust=0.5),
        axis.title.y = element_text(size=24, vjust=0.5),
        strip.text = element_text(size=24)) +
  ylab("No. of significant TE families") + 
  xlab(expression(paste("Minimum insertion difference cut-off (min. |", delta, "Insertions|)")) ) + theme(legend.position = "none") + 
  facet_wrap(~Study, ncol= 2,scales="free")
barplot.eff = barplot.eff + 
  annotate("text", x = 20, y = 90, label = "S>C",size = 7,col = "red") + 
  annotate("text", x = 20, y = 80, label = "C>S",size = 7,col="blue") + 
  geom_text( aes( label = abs(value), y = labs.bar),
             vjust = 0, size = 4, color = "black" )
barplot.eff = barplot.eff + theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"))
barplot.eff
ggsave(barplot.eff, file="Effect_size_cut_off.pdf", height=12, width=14)

#Check HeT-A, G2, G-element
all.list = list(carnes.stat.filt,
                fabian.stat.filt,
                hoedjes.stat.filt,
                remo.stat.filt)

lapply(all.list, function(x) x[x$flybase_name %in% c("HeT-A","G2","G-element"),]$Diff_SelCont)
"HeT-A","G2","G-element"
0.334 -3.582 -4.161 Carnes
1.097 -0.508 -0.532 Fabian
4.433 -0.824 -0.556 Hoedjes
0.443 -0.501 -0.511 Remolina


lapply(all.list, function(x) x[x$flybase_name %in% c("HeT-A","G2","G-element"),]$Diff_PostEarly)

#Het-A 0.334 in Carnes2015
#Other 2 are >0.5

#Cut-off of 0.3
# Remolina
# > 1 - (10 + 41) / (13+63)
# [1] 0.3289474
# Hoedjes
# > 1 - (11 + 41) / (27+67)
# [1] 0.4468085
# Fabian
# > 1 - (31 + 37) / (34+51)
# [1] 0.2
# Carnes
# > 1 - (21 + 82) / (21+86)
# [1] 0.03738318
mean( c(0.3289474,0.4468085,0.2,0.03738318))
#25%

#Plot proportions
all.eff.size.melt.prop = melt(data = all.eff.size[,c(-4)],id.vars = c("Eff_size","Study"),measure.vars = c("SCprop","CSprop"))

labs.bar.prop = c( all.eff.size.melt.prop[all.eff.size.melt.prop$variable == "SCprop",]$value +0.02, all.eff.size.melt.prop[all.eff.size.melt.prop$variable == "CSprop",]$value -0.05)

barplot.eff.prop = ggplot(data=all.eff.size.melt.prop , aes(x=Eff_size, y=value, fill=variable)) +
  geom_bar(stat="identity", color = "black") + THEME_DEF +
  scale_y_continuous(limits=c(-0.4,0.8),breaks=seq(-0.4,0.8,0.2), labels= abs(seq(-0.4,0.8,0.2))) + 
  scale_fill_manual(values=c("red","blue")) + geom_hline(yintercept = 0) +
  theme(axis.text.y = element_text(size=18, colour = "black"),
        axis.text.x = element_text(size=18 ,angle = 45,hjust=1),
        axis.title.y = element_text(size=24, vjust=0.5),
        strip.text = element_text(size=15)) +
  ylab("Proportion of significant TE families") + xlab(expression(paste("Effect Size Cutoff: absolute ", delta, "Insertions (S - C)")) ) + theme(legend.position = "none") + 
  facet_wrap(~Study, ncol= 2,scales="free_y")
barplot.eff.prop = barplot.eff.prop + 
  annotate("text", x = 20, y = 0.8, label = "S>C",size = 7,col = "red") + 
  annotate("text", x = 20, y = 0.7, label = "C>S",size = 7,col="blue") + 
  geom_text( aes( label = abs(all.eff.size.melt$value), y = labs.bar.prop),
             vjust = 0, size = 4, color = "black" )
barplot.eff.prop
ggsave(barplot.eff.prop, file="Effect_size_cut_off_proportion.pdf", height=12, width=14)

