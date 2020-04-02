#Plot of differentially expressed and genetically differentiated TE regulation candidate genes
#For Figure 5

#Please change folder names and edit commands accordingly.
source("/Users/danfab/R_functions.R")

setwd("/Users/danfab/comp_all/")

#Barcharts for number of differemtially expressed TE regulation genes
de_teReg = read.table("/Users/danfab/extra_files/number_prop_DE_TEregulationGENES_Main_factors.txt",header=T) #Created table in excel
de_teReg$Level = factor(de_teReg$Level, c("young","old","control","selected","f","m","f1","m1"))
de_teReg.remo = de_teReg[de_teReg$Study == "Remolina",]
de_teReg.carnes = de_teReg[de_teReg$Study == "Carnes",]

de_teReg.carnes
mat.carnes = t(rbind(c(15,4),c(6,47),c(8,NA),
                     c(22,7),c(4,6),c(11,NA)))
mat.remo = t(rbind(c(0,0),c(0,5),c(0,NA),
                   c(0,0),c(1,23),c(0,NA)))
#colnames(mat.carnes) = c("Regime (M)", "Age (M)","RxA (M)" ,"Regime (F)","Age (F)", "RxA (F)")
rownames(mat.carnes) = c("up","down")

pdf("DiffExpr_TEregulGenes_Carnes.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=1.8,mar=c(8,4.6,2,2))
xx = barplot(mat.carnes,beside = T, col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
                                            "blue","red","goldenrod2", "goldenrod4","green4","green4"),
             ylab = "",
             ylim=c(0,50), space=c(0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.carnes,beside = T, col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
                                       "blue","red","goldenrod2", "goldenrod4","green4","green4"),
        ylab = "No. of TE Regulation Genes (D.E.)",
        ylim=c(0,50), space=c(0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0),add=T)
# legend("topright",legend = c("C>S","S>C","Young>Old","Old>Young","Regime x Age"), 
#                fill= c("blue","red","goldenrod2", "goldenrod4","green4"),cex=1.5)
text(x = xx, y = mat.carnes, label = mat.carnes, pos = 3, cex = 1.2, 
     col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
             "blue","red","goldenrod2", "goldenrod4","green4","green4"))
dev.off()

pdf("DiffExpr_TEregulGenes_Remo.pdf", width=6, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=1.8,mar=c(8,4.6,2,2))
xx = barplot(mat.remo,beside = T, col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
                                          "blue","red","goldenrod2", "goldenrod4","green4","green4"),
             ylab = "",
             ylim=c(0,50), space=c(0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0,
                                   0.6,0))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.remo,beside = T, col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
                                     "blue","red","goldenrod2", "goldenrod4","green4","green4"),
        ylab = "",
        ylim=c(0,50), space=c(0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0,
                              0.6,0), add=T)
legend("topright",legend = c("C>S","S>C","Young>Old","Old>Young","Regime x Age"), 
       fill= c("blue","red","goldenrod2", "goldenrod4","green4"),bg = "white",cex=1.7)
text(x = xx, y = mat.remo, label = mat.remo, pos = 3, cex = 1.2, 
     col = c("blue","red","goldenrod2", "goldenrod4","green4","green4",
             "blue","red","goldenrod2", "goldenrod4","green4","green4"))
dev.off()


#Plot genetically differentiated candidate genes
mat.genetic_diff = c(10,3,7,6) #Numbers of candidate genes associated to TE regulation
pdf("GenetCandidate_TEregulGenes_all.pdf", width=4, height=8)
par(mfrow = c(1,1),font=1,font.lab=1,font.axis=1,cex.lab=2.2,cex.main = 2,cex.axis=1.8,mar=c(8,4.6,2,2))
xx = barplot(mat.genetic_diff,beside = T,
             ylab="No. of TE Regulation Genes (G.D.)",
             ylim=c(0,10.5))
grid(col = "lightgray", lty = 1, lwd = 0.5, equilogs = TRUE)
barplot(mat.genetic_diff,beside = T,
        ylab = "",
        ylim=c(0,10),add=T,col="gray29")
text(xx, par("usr")[3]-0.25, 
     srt = 44, adj= 1, xpd = TRUE,
     labels = c("Carnes2015","Fabian2018","Hoedjes2019", "Remolina2012"), cex=1.8)
text(x = xx, y = mat.genetic_diff, label = mat.genetic_diff, pos = 3, cex = 1.2)
dev.off()
