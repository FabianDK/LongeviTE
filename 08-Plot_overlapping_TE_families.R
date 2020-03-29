#Plot TE families shared between studies

#For Figure 3B and Figure S1

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

library(tidyr)
library(lemon)
source("/Users/danfab/Science/post-doc/R_custom_functions/R_functions.R")

output_path = "/Users/danfab/comp_all/"

#Annotation tab
ensembl_casey = read.table("/Users/danfab/extra_files/embl_repbase_mapping_from_Bergman_edit2.txt",header = T)
ensembl_casey$flybase_name = gsub("Dmel/","",ensembl_casey$flybase_name,fixed = T)
ensembl_casey = ensembl_casey[,c(1,4)]

mean_ins_tab = read.table("/Users/danfab/comp_all/Table_Mean_Insertions_Per_Family_With_TotalCont.txt", header=T,check.names = F)
mean_ins_tab = mean_ins_tab[,-5] #drop diet

TE.neg.complete = read.table("/Users/danfab/comp_all/shared_C>S_table.txt", header=T,check.names = F)
TE.pos.complete = read.table("/Users/danfab/comp_all/shared_S>C_table.txt", header=T,check.names = F)
head(TE.neg)
TE.neg = TE.neg.complete[rowSums(TE.neg.complete[,2:5])==4,]
TE.neg
TE.pos = TE.pos.complete[rowSums(TE.pos.complete[,2:5])==4,]
TE.pos

#Edit tables
col.sel.pos = c("study","regime","rep",as.character(TE.pos[,1]))
sel.tab.pos = mean_ins_tab[,names(mean_ins_tab) %in% col.sel.pos ]
sel.tab.pos.melt = sel.tab.pos %>% gather(key = 'TE', value = 'Insertions',-study,-regime, -rep) 
names(sel.tab.pos.melt)[4] = "TEfam"
nrow(sel.tab.pos.melt) #644
sel.tab.pos.melt = merge(ensembl_casey, sel.tab.pos.melt, by ="TEfam")

col.sel.neg = c("study","regime","rep",as.character(TE.neg[,1]))
sel.tab.neg = mean_ins_tab[,names(mean_ins_tab) %in% col.sel.neg ]
sel.tab.neg.melt = sel.tab.neg %>% gather(key = 'TE', value = 'Insertions',-study,-regime,-rep) 
names(sel.tab.neg.melt)[4] = "TEfam"
nrow(sel.tab.neg.melt) #92
sel.tab.neg.melt = merge(ensembl_casey, sel.tab.neg.melt, by ="TEfam")
nrow(sel.tab.neg.melt) #92


#Plot defs - per study
THEME_DEF_BOX = theme(axis.title.y = element_text(vjust=0.5,size=28),
                      axis.title.x = element_text(vjust=-0.2,size=18),
                      axis.text.x = element_text(size=18, colour = "black"),
                      axis.text.y = element_text(size=18, colour = "black"),
                      axis.ticks = element_line(size = 1, colour = "black"),
                      axis.ticks.length = unit(0.4, "cm"),
                      axis.line = element_line(size = 1),
                      plot.title = element_text(color="black", size=20, face="bold.italic"),
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      legend.key.size = unit(1.2, "cm"),
                      legend.text = element_text(size=20),
                      legend.title = element_text(size=20,face="bold"),
                      legend.position="right",
                      panel.background = element_blank(),
                      panel.border = element_rect(fill = NA, colour="black", size = 1),
                      panel.grid.major = element_line(),
                      panel.grid.minor = element_blank())
legend_title = "Regime"

p.pos = ggplot(data = sel.tab.pos.melt, aes(x=study, y=Insertions,fill=regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(legend_title,values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions"))) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.pos = p.pos + 
   facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

p.neg = ggplot(data = sel.tab.neg.melt, aes(x=study, y=Insertions,fill=regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(legend_title, values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions")))+
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.neg = p.neg + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))
p.pos 
pos.save.path = paste(output_path,"shared_S>C_TEs.pdf",sep="")
neg.save.path = paste(output_path,"shared_C>S_TEs.pdf",sep="")
ggsave(p.pos, file=pos.save.path, height=14, width=22)
ggsave(p.neg, file=neg.save.path, height=6, width=14)

#Per regime
p.pos = ggplot(data = sel.tab.pos.melt, aes(x=regime, y=Insertions, fill=regime)) + theme_bw() + geom_boxplot() + THEME_DEF_BOX + 
  scale_fill_manual(legend_title,values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions"))) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.pos = p.pos + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

p.neg = ggplot(data = sel.tab.neg.melt, aes(x=regime, y=Insertions, fill=regime)) + theme_bw() + geom_boxplot() + THEME_DEF_BOX + 
  scale_fill_manual(legend_title, values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions"))) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.neg = p.neg + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))
p.pos 
pos.save.path = paste(output_path,"shared_S>C_TEs_Regime.pdf",sep="")
neg.save.path = paste(output_path,"shared_C>S_TEs_Regime.pdf",sep="")
ggsave(p.pos, file=pos.save.path, height=14, width=14)
ggsave(p.neg, file=neg.save.path, height=7, width=12)


#####################
#SAME BUT NORMALIZED BY TOTAL CONTENT
#####################

mean_ins_tab.norm = cbind(mean_ins_tab[,2:4],mean_ins_tab[,5:ncol(mean_ins_tab)] / mean_ins_tab$Ins)
head(mean_ins_tab.norm)

col.sel.pos = c("study","regime","rep",as.character(TE.pos[,1]))
sel.tab.pos = mean_ins_tab.norm[,names(mean_ins_tab.norm) %in% col.sel.pos ]
sel.tab.pos.melt = sel.tab.pos %>% gather(key = 'TE', value = 'Insertions',-study,-regime, -rep) 
names(sel.tab.pos.melt)[4] = "TEfam"
nrow(sel.tab.pos.melt) #644
sel.tab.pos.melt = merge(ensembl_casey, sel.tab.pos.melt, by ="TEfam")
nrow(sel.tab.pos.melt) #644

col.sel.neg = c("study","regime","rep",as.character(TE.neg[,1]))
sel.tab.neg = mean_ins_tab.norm[,names(mean_ins_tab.norm) %in% col.sel.neg ]
sel.tab.neg.melt = sel.tab.neg %>% gather(key = 'TE', value = 'Insertions',-study,-regime,-rep) 
names(sel.tab.neg.melt)[4] = "TEfam"
nrow(sel.tab.neg.melt) #92
sel.tab.neg.melt = merge(ensembl_casey, sel.tab.neg.melt, by ="TEfam")
nrow(sel.tab.neg.melt) #92

#Plot defs - per study
THEME_DEF_BOX = theme(axis.title.y = element_text(vjust=0.5,size=28),
                      axis.title.x = element_text(vjust=-0.2,size=18),
                      axis.text.x = element_text(size=18, colour = "black"),
                      axis.text.y = element_text(size=18, colour = "black"),
                      axis.ticks = element_line(size = 1, colour = "black"),
                      axis.ticks.length = unit(0.4, "cm"),
                      axis.line = element_line(size = 1),
                      plot.title = element_text(color="black", size=20, face="bold.italic"),
                      legend.background = element_blank(),
                      legend.key = element_blank(),
                      legend.key.size = unit(1.2, "cm"),
                      legend.text = element_text(size=20),
                      legend.title = element_text(size=20,face="bold"),
                      legend.position="right",
                      panel.background = element_blank(),
                      panel.border = element_rect(fill = NA, colour="black", size = 1),
                      panel.grid.major = element_line(),
                      panel.grid.minor = element_blank())
legend_title = "Regime"

p.pos = ggplot(data = sel.tab.pos.melt, aes(x=study, y=Insertions, fill=regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(legend_title,values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Normalized Insertions"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.pos = p.pos + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))
p.pos

p.neg = ggplot(data = sel.tab.neg.melt, aes(x=study, y=Insertions,fill=regime)) + theme_bw() + geom_boxplot(outlier.shape = NA) + THEME_DEF_BOX + 
  scale_fill_manual(legend_title, values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Normalized Insertions"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.neg = p.neg + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))
p.neg
pos.save.path = paste(output_path,"shared_S>C_TEs_normalized.pdf",sep="")
neg.save.path = paste(output_path,"shared_C>S_TEs_normalized.pdf",sep="")
ggsave(p.pos, file=pos.save.path, height=14, width=22)
ggsave(p.neg, file=neg.save.path, height=6, width=14)

#Per regime
p.pos = ggplot(data = sel.tab.pos.melt, aes(x=regime, y=Insertions, fill=regime)) + theme_bw() + geom_boxplot() + THEME_DEF_BOX + 
  scale_fill_manual(legend_title,values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions / Total TE content"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.pos = p.pos + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

p.neg = ggplot(data = sel.tab.neg.melt, aes(x=regime, y=Insertions,fill=regime)) + theme_bw() + geom_boxplot() + THEME_DEF_BOX + 
  scale_fill_manual(legend_title, values = c("blue", "red")) + 
  labs(x = expression(bold("")), y = expression(bold("Insertions / Total TE content"))) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25),size=2.8,shape=21,col="white")
p.neg = p.neg + 
  facet_rep_wrap(~flybase_name, scales='free_y', repeat.tick.labels = 'left') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 20))
p.pos 
pos.save.path = paste(output_path,"shared_S>C_TEs_Regime_normalized.pdf",sep="")
neg.save.path = paste(output_path,"shared_C>S_TEs_Regime_normalized.pdf",sep="")
ggsave(p.pos, file=pos.save.path, height=14, width=14)
ggsave(p.neg, file=neg.save.path, height=7, width=12)
