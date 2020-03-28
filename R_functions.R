#Libraries often used
library(reshape2)
library(ggplot2)
library(patchwork) #For ggplots, e.g. p1 + p2
library(FactoMineR)
library(ggpubr)

#Functions for R

#plot a transparent color
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#Transform a factor of numbers to numeric
numfact = function(x){
  transform=as.numeric(as.character(x))
  return(transform)
}

#Show both, minimum and maximum
minmax = function(x){
  c(min(x),max(x))
}

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

#Extract pvalues and coeff from coxme result (object)
extract_coxme <- function (mod){
  beta <- fixef(CPH)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

###layout function for assmebling multiple ggplots
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 
#lay_out(list(snphistPerc+theme(plot.margin = unit(c(0.5,0.5,0.2,0.7), "cm")) , 1, 1), list(fst_PLOT_X, 1, 2), list(fst_PLOT_2L, 2.1, 1),list(fst_PLOT_2R, 2.1, 2),list(fst_PLOT_3L, 3.2, 1),list(fst_PLOT_3R, 3.2, 2))


#FOR OVERLAPS USING SUPEREXACTTEST
#Input is a list of multiple candidate gene table - each one must have a column with "Entrez_Gene_ID" to work
superOverlap = function(candiset_list,bg_size){
  require(SuperExactTest)
  setlist = NULL
  candi_number = NULL
  for (i in 1:length(candiset_list)) {
    setlist[[i]]=as.character(candiset_list[[i]]$Entrez_Gene_ID) #needs name: Entrez_Gene_ID
    #candi_number = c(candi_number, length(candiset_list[[i]]$Entrez_Gene_ID))
  }
  names(setlist) = names(candiset_list)
  #names(candi_number) = paste(names(candiset_list),"_candi_number", sep="")
  stat = supertest(setlist,bg_size)
  #stat.edit = c(stat,overlap_size = length(stat$intersects),candi_number)
  return(summary(stat))
}

#bighead function to check big tables instead of head
bighead = function(x, rowstart,rowend,colstart,colend){
  x[rowstart:rowend,colstart:colend]
}

#check both row and column
rowcol = function(x){
  intvect = c(nrow(x),ncol(x))
  names(intvect) = c("rows","cols")
  return(intvect)
}
 
#unique and non-unique function
uniqdup = function(df,col,method){
  if (method == "uniq") {
    subdf = df[!(df[,col] %in% df[,col][duplicated(df[,col])]),]
    return(subdf)
  }
  if (method == "dupl") {
    subdf = df[df[,col] %in% df[,col][duplicated(df[,col])],]
    return(subdf)
  }
}

#Extract tested number of nodes from topGO object - input topGO object, output number of tested nodes
topGO_nodes = function(topGO_object){
  as.numeric(gsub("   - number of nodes = ","",capture.output(topGO_object)[22]))
}

#Useful:
#Add a few more colors to palette
palette("default")
palette(c(palette(),"orange","lightpink","gray47","khaki3","indianred1","tan","steelblue1","springgreen","white","brown4","orange4","burlywood"))

#get all genes in a GO term, including child terms
#1) species is selected from listDatasets(useMart("ensembl"))
#2) GO is a character vector of GO terms
#3) set.annot is a character vector with field names to be extracted, e.g. from flies: listAttributes(useMart("ensembl", dataset = 'dmelanogaster_gene_ensembl'))
get.genes.go =
  function(species, GO, annot){ 
    ensMart = useMart("ensembl", dataset = species)
    outputlist = list()
    for (i in 1:length(GO)) {
      outputlist[[i]] = getBM(annot, 
                     filters = 'go_parent_term', 
                     values = GO[i], 
                     mart = ensMart)
    }
    names(outputlist) = GO
    return(outputlist)
  }

#unique and length combined
luniq = function(x){
  return(length(unique(x)))
}

#Get CPM and RPKM values from read counts 
getRPKM = function (read_counts, gene_lengths) {
  sumtab = NULL
  scale_fact = sum(read_counts) / 1000000
  cpm = read_counts / scale_fact #counts per million
  rpkm = cpm / (gene_lengths / 1000) #divide by kilobases

  reads_per_kb = read_counts / gene_lengths #RPK
  reads_per_kb_per_mil = sum(reads_per_kb) / 1000000
  tpm = reads_per_kb / reads_per_kb_per_mil
  
  sumtab$cpm = cpm
  sumtab$rpkm = rpkm
  sumtab$tpm = tpm
  #sumtab[,1] = format(sumtab[,1],scientific = F)
  return(sumtab)
}

#Get colors for denisty plot - x and y are columns of scatter plot
# colorPal = c("cyan", "blue", "green", "yellow", "red")
density.colors = function(x,y,colorPal){
  denisty.map = densCols(x,y,colramp=colorRampPalette(c("black", "white")))
  cols.dens = col2rgb(denisty.map)[1,] + 1L
  cols =  colorRampPalette(colorPal)(256)
  cols.final = cols[cols.dens]
  return(cols.final)
}

#Create seperate ggplot legend
gglegend <- function(x){ 
  tmp <- ggplot_gtable(ggplot_build(x)) 
  leg <- which(sapply(tmp$grobs, function(y) y$name) == "guide-box") 
  tmp$grobs[[leg]]
}
# legend = gglegend(plot)
#draw it with grid.draw(legend)

#rowMins and max
rowMin = function(df,colstart,colend){ 
  rowMin.vect = NULL
  for (i in 1:nrow(df)){
    int = df[i,colstart:colend,]
    rowMin.vect = c(rowMin.vect,min(df[i,colstart:colend,]))
  }
  return(rowMin.vect)
}

rowMax = function(df,colstart,colend){ 
  rowMax.vect = NULL
  for (i in 1:nrow(df)){
    int = df[i,colstart:colend,]
    rowMax.vect = c(rowMax.vect,max(df[i,colstart:colend,]))
  }
  return(rowMax.vect)
}

#DEseq2 - get Transposable elements
deseq.results = function(DEseqobj, TEvect, cutoff, type){
  res.TE.list = NULL
  res.TE.summary = NULL
  for (i in 2:length(resultsNames(DEseqobj))){
    extr.fact = resultsNames(DEseqobj)[i]
    res.int = as.data.frame(results(DEseqobj,list( extr.fact) )) #cooksCutoff=T,independentFiltering=T
    res.int.TE = res.int[row.names(res.int) %in% TEvect,]
    res.TE.list[[extr.fact]] = res.int.TE
    res.int.TE.padj = as.vector(na.omit(res.int.TE$padj))
    sign.n = res.TE.list[[extr.fact]][(res.TE.list[[extr.fact]]$padj<cutoff),]
    sign.n.FC = sign.n[sign.n$log2FoldChange > 0,]
    sign.n.FC.count = nrow(na.omit(sign.n.FC))
    sign.n.count = nrow(na.omit(sign.n))
    ns.n = res.TE.list[[extr.fact]][(res.TE.list[[extr.fact]]$padj>=cutoff),]
    ns.n.count = nrow(na.omit(ns.n))
    res.TE.summary = as.data.frame(rbind(res.TE.summary, 
                                         cbind(extr.fact,
                                               sign.n.count,
                                               ns.n.count,
                                               sign.n.FC.count,
                                               sign.n.count-sign.n.FC.count,
                                               sign.n.count+ns.n.count)))
  }
  
  names(res.TE.list)
  names(res.TE.summary) = c("Factor","Padj005","ns","FC>0","FC<0","Total")
  res.TE.summary = as.data.frame(res.TE.summary)
  
  if (type == "summary")
    { return(res.TE.summary)
    } else (type == "table")
  { return(res.TE.list)
  }
}

#get first n and last n rows
headtail = function(x, n){
  headtail.tab =rbind(head(x,n),tail(x,n))
  return(headtail.tab)
}

#ggplot theme
THEME_DEF = theme(axis.title.y = element_text(vjust=0.5,size=20),
                  axis.title.x = element_text(vjust=-0.2,size=20),
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

THEME_DEF_small = theme(axis.title.y = element_text(vjust=0.5,size=16),
                  axis.title.x = element_text(vjust=-0.2,size=16),
                  axis.text.x = element_text(size=14, colour = "black"),
                  axis.text.y = element_text(size=14, colour = "black"),
                  axis.ticks = element_line(size = 1, colour = "black"),
                  axis.ticks.length = unit(0.4, "cm"),
                  axis.line = element_line(size = 1),
                  plot.title = element_text(color="black", size=20, face="bold.italic"),
                  legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
                  legend.key = element_blank(),
                  legend.key.size = unit(1, "cm"),
                  legend.text = element_text(size = 16),
                  legend.title = element_text(size = 18, face = "bold"),
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, colour="black", size = 1),
                  panel.grid.major = element_line(colour="grey85", size = (0.2)),
                  panel.grid.minor = element_blank())

#P-value symbols for ggpubr package
symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))


#SNPEFF OUTPUT EDITING FUNCTION
#df = "EFF=UTR_3_PRIME(MODIFIER||91||655|Nepl6|protein_coding|CODING|FBtr0299969|1|1),DOWNSTREAM(MODIFIER||30||684|Nepl5|protein_coding|CODING|FBtr0079255||1),INTRON(MODIFIER||||467|DIP-epsilon|protein_coding|CODING|FBtr0299967|7|1)"
#df = "EFF=NON_SYNONYMOUS_CODING(MODERATE||gagtacatt/gaRAL-NON|EYI26???|166|CG9240|protein_coding|CODING|FBtr0074060||1),UTR_5_PRIME(MODIFIER||212||128|CG9240|protein_coding|CODING|FBtr0308614|1|1|WARNING_REF_DOES_NOT_MATCH_GENOME),DOWNSTREAM(MODIFIER||447||224|Pis|protein_coding|CODING|FBtr0074058||1),DOWNSTREAM(MODIFIER||447||81|Pis|protein_coding|CODING|FBtr0074059||1)"
#df = "EFF=EXON(MODIFIER|||||CR18228||NON_CODING|FBtr0310082|2|1),EXON(MODIFIER|||||CR18228||NON_CODING|FBtr0300631|2|1)"
snpeff.edit = function(df){
  df[,8] = as.character(df[,8])
  df[,8] = mgsub(c("(LOW","(MODERATE","(HIGH","(MODIFIER","EFF=","|WARNING_REF_DOES_NOT_MATCH_GENOME",
                          "|WARNING_TRANSCRIPT_INCOMPLETE",
                          "|WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS",
                          "|WARNING_TRANSCRIPT_NO_START_CODON"),
                        c("","","","","","","","",""),df[,8],fixed = T)
  
  df.edit = NULL
  for (i in 1:nrow(df)){
    # split.col = unlist(strsplit(x = df[i,8],split = ","))
    #sometimes additional stuff to EFF, like LOF or NMD predicted, remove this
    split.col = unlist(strsplit(x = unlist(strsplit(x = df[i,8],split = ";"))[1],split = ","))
    split.col
    split.col.s =  strsplit(split.col,split = "|",fixed = T)
    split.col.s
    int.tab = matrix(unlist(split.col.s),
                     nrow=length(split.col.s),
                     ncol=11,byrow=T)[,c(1,3,6,7,8)]
    int.tab = unique(int.tab)
    if (length(split.col.s) == 1){
      int.tab=matrix(int.tab,nrow=1,ncol=5) 
    }
    # if (int.tab[,2] == ""){
    #   int.tab[,2] = "NA"
    # }
    
    int.tab[,2][int.tab[,2]==""] = "NA"
    int.tab[,4][int.tab[,4]==""] = "NA"
    
    int.tab
    if ("INTERGENIC" %in% c(int.tab)){
      int.tab = matrix(gsub("^$|^ $", NA, int.tab), 
                       ncol=ncol(int.tab), nrow = nrow(int.tab), byrow=F)
    }
    int.tab[,3] = gsub("INTERGENIC",NA,int.tab[,3],fixed = T)
    int.tab[,5] = gsub("INTERGENIC",NA,int.tab[,5],fixed = T)
    
    int.tab2 = cbind(df[i,c(1,2,4:7,9:ncol(df))],
                     int.tab)
    int.tab2
    df.edit = rbind(df.edit,
                    unique(int.tab2))
  }
  
  df.edit$`3` = gsub("Gene_","",df.edit$`3`)
  df.edit$`3` = sub("_","(",df.edit$`3`)
  df.edit$`3` = sub("_",")",df.edit$`3`)
  return(df.edit)
}



#MANHATTAN PLOTS
# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, hlight2, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_highlight2=ifelse(SNP %in% hlight2, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig),color = "orange", size = 1) +
    geom_hline(yintercept = -log10(sugg), color = "orange", linetype="dashed") +
    
    # Add highlighted points
    geom_point(data=subset(df.tmp, is_highlight2=="yes"), bg="blue", size=5,pch=21) +  #col="blue"
    geom_point(data=subset(df.tmp, is_highlight=="yes"), bg="red", size=5,pch=21 ) + #col="red"
    
    # Add label using ggrepel to avoid overlapping
    # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=1.5, force=1.3) +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

#List to binary vector matrix for overlaps
list_to_matrix = function(lt, universal_set = NULL) {
  if(!is.null(universal_set)) {
    lt = lapply(lt, function(x) intersect(x, universal_set))
  } else {
    universal_set = unique(unlist(lt))
  }
  
  mat = matrix(0, nrow = length(universal_set), ncol = length(lt))
  rownames(mat) = sort(universal_set)
  colnames(mat) = names(lt)
  for(i in seq_along(lt)) {
    mat[unique(lt[[i]]), i] = 1
  }
  return(mat)
}

