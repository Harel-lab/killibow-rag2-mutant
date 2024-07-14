library(edgeR)            #3.34.1
library(ggfortify)        #0.4.13
library(ggplot2)          #3.3.6
library(stats)            #4.1.1
library(factoextra)       #1.0.7
library(ggConvexHull)     #0.1.0
library(gridExtra)        #2.3
library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0
library(ComplexHeatmap)   #2.8.0
library(ggrepel)          #0.9.1    
library(reshape2)         #1.4.4
library(ggpubr)           #0.4.0
library(tidyr)            #1.2.0
library(cowplot)          #1.1.1

WTcol = '#83deb5'
RAGcol = '#00203fff' 

# sources----
path <- 'C:/Users/tehil/Dropbox/Projects/Rag2/'
expressionPath <- 'rag2_matrix_expression_RNAseq_star2206.csv' 
sampleInfoPath <- 'infoTable.csv'
humanKfConversion = read.table( "C:/Users/tehil/Dropbox/Projects/NCBI-Human-orthologs.txt", head = T, sep = "\t")

#function
createDEGobject <- function(path, expPath, infoPath){
  # take two table- 1. RNA seq with genes on the columns and samples in row 2. information table- age, genotype, feed and sex for each sample
  # create DEG object, order the samples(first - young WT male full)
  # return the DE-object
  
  # read the expression (counts) table and order the names of genes
  exMatrixT <- read.csv(paste0(path, expPath), row.names = 1)
  exMatrix <- data.frame(t(exMatrixT))
  colnames(exMatrix) <- row.names(exMatrixT)
  
  # read the samples information table, make the parameter as factor and re-factor them(wt, male, full, young as default)
  infoTable <- read.csv(paste0(path, infoPath), row.names = 1, stringsAsFactors=T)
  infoTable$Genotype <- relevel(infoTable$Genotype, ref = "WT")
  #infoTable$sex <- relevel(infoTable$sex, ref = "male")
  #infoTable$age <- relevel(infoTable$age, ref = "6.5weeks")
  #infoTable$feed <- relevel(infoTable$feed, ref = "full")
  infoTable$sample <- row.names(infoTable)
  
  #infoTable <- infoTable[, names(infoTable) %in% c("age", "genotype", "feed", "sex", "sample")]
  infoTable <- infoTable[colnames(exMatrix),] # very important line!! order the info according the RNA-seq data; the first column in count table will be corresponded to the first row in sample table
  
  # create DGE object     
  DEobj <- DGEList(exMatrix, samples = infoTable, group=infoTable$Genotype,
                   genes = row.names(exMatrix))
  
  #DEobj <- DEobj[,order(DEobj$samples$genotype)]
  #DEobj <- DEobj[,order(DEobj$samples$age)]
  #DEobj <- DEobj[,order(DEobj$samples$feed)]
  #DEobj <- DEobj[,order(DEobj$samples$sex)]
  #DEobj$samples$group <- factor(DEobj$samples$group, levels = unique(DEobj$samples$group))
  return(DEobj)
}
filterNormDEGobject <-function(DEobj){
  # filtering low expressed genes
  # normalized each sample by weight average of non-DE genes using TMM
  
  # filtering out lower expressed genes
  keep <- filterByExpr(DEobj, model.matrix(~Genotype, data=DEobj$samples), min.count = 7)
  print(table(keep))
  DEobj <- DEobj[keep,, keep.lib.sizes=FALSE]
  
  #Normalize data using TMM
  # TMM- 
  DEobj <- calcNormFactors(DEobj, method = 'TMM')
  par(mfrow = c(1, 2))
  boxplot(cpm(DEobj, normalized.lib.sizes = F, log = T), outline=F, col='white', main='Before')
  boxplot(cpm(DEobj, normalized.lib.sizes = T, log = T), outline=F, col='white', main='After')
  return(DEobj)
}
dataExplorer <- function(DEobj, titleG='', logTag = T, clusterSamples = F){
  # dataExplorer accept DGE object
  # return PCA (PC1/2 and PC1/3) with and without labels, and figure of variances against the number of dimensions
  
  #using log2(CPM)
  cpm.data <- log2(cpm(DEobj, normalized.lib.sizes = T)+1)
  cpm.data <- cpm.data[apply(cpm.data, 1, sd) != 0,]
  infoTable <- DEobj$samples  #[,c('age', 'genotype', 'feed', 'sex', 'sample', 'group')]
  
  #*****************PCA************
  #https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
  pca.de <- prcomp(t(cpm.data), scale. = T)
  fz <-fviz_eig(pca.de)
  # PCA with labels PC1/2 and PC1/3
  pcaL12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2, size = 2.5, label = F, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    scale_color_manual(name='Group', values = c(WTcol, RAGcol))+
    scale_fill_manual(name='Group', values = c(WTcol, RAGcol))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaL13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3, size = 2.5, label = F, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    scale_color_manual(name='Group', values = c(WTcol, RAGcol))+
    scale_fill_manual(name='Group', values = c(WTcol, RAGcol))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  "#PCA with dots
  pcaD12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaD13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#1D71BB', '#E30613', '#1D71BB', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13, PCAdots12 = pcaD12, PCAdots13 = pcaD13))"
  return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13, jjj=pca.de))
}
convertLoc2symbol <- function(old_names){
  names_ncbi <- read.csv('C:/Users/tehil/Dropbox/Projects/genes_names.csv')
  names(names_ncbi) <- c('NCBI', 'FinalSymbol', 'Human')
  head(names_ncbi)
  new_names <- names_ncbi[names_ncbi$NCBI %in% old_names,]
  new_names <- new_names[match(old_names, new_names$NCBI),]$FinalSymbol
  return(new_names)
}


# creating DGE object
DEobj <- createDEGobject(path, expressionPath, sampleInfoPath)

DEobjRag2 <- filterNormDEGobject(DEobj)
dev.off()

Rag2Graphs <- dataExplorer(DEobjRag2, 'Rag2')

#save the PCA graphs
pdf(paste0(path, '/figures/PCA RNAseq2.pdf'), width=10, height=10) 
grid.arrange(Rag2Graphs$PCAlabels12, Rag2Graphs$PCAlabels13, 
             Rag2Graphs$knee, nrow = 2)
dev.off()


# -------DE------
# Differential expression analysis between ages, between genotypes and in the interaction between them in the three datasets

#design <- model.matrix(~Genotype, data=DEobjRag2$samples)
DEobjRag2 <- estimateDisp(DEobjRag2)#), design)
et <- exactTest(DEobjRag2)
top <- topTags(et, n=Inf)$table
top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
top <- merge(top, humanKfConversion, by.x='genes', by.y = 'ncbi', all.x = TRUE)

write.csv(top, paste(path, 'ranking list rag2 classic.csv'))
write.csv(top[top$FDR <0.1, ], paste(path, 'ranking list rag2 classic 0.1.csv'))

dim(top[top$PValue < 0.05 & top$logFC > 0,])
dim(top[top$PValue < 0.05 & top$logFC < 0,])

dim(top[top$FDR < 0.05 & top$logFC > 0,])
dim(top[top$FDR < 0.05 & top$logFC < 0,])
head(top[order(top$FDR), ])

#GO
enrichmentClusterProfile <- function(geneListHuman, folder, fileName, backgroundGeneList = 'none') {
  
  # geneListHuman: human gene list by Symbol annotation
  # folder, fileName:  folder name and fileName for saving the files.
  # backgroundGeneList: (optional) background list. if provide use it to the universe.
  
  if (length(backgroundGeneList) == 1) {
    backgroundGeneList <- NULL
    entrezIdBG <- NULL
  } else
    entrezIdBG = bitr(as.character(backgroundGeneList), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID # Get entrez ids for annotation
  
  # go 
  ego3 <- enrichGO(gene        = geneListHuman,
                   OrgDb        = org.Hs.eg.db,
                   universe     = backgroundGeneList,
                   keyType      = 'SYMBOL',
                   ont          = c("ALL"),
                   pvalueCutoff = 1)
  
  print(length(ego3))
  if (length(ego3) > 0)
    if (nrow(ego3) > 0)
      write.table(ego3, paste0(folder, fileName, "_GO.csv"), sep = ",", quote = T, row.names = F)
  return(ego3)
}

enrichmentClusterProfile(top[top$FDR < 0.05,]$human, path, 'GO_0.1', top$human)
up = enrichmentClusterProfile(top[top$FDR < 0.05 & top$logFC > 1,]$human, path, 'GO_0.05_FC1_up', top$human)
down = enrichmentClusterProfile(top[top$FDR < 0.05 & top$logFC < -1,]$human, path, 'GO_0.05_FC1_down', top$human)


## visualization ----
cpmRag2 <- log2(cpm(DEobjRag2, normalized.lib.sizes = T)+1)

barplotExpression <- function(gene, norma = TRUE){
  expressionG <- data.frame(Exp = t(cpmRag2[gene,]), 
                     sample=DEobj$samples$sample, Genotype=DEobj$samples$Genotype)
  colnames(expressionG) = gsub('Exp.', '', colnames(expressionG))
  ylabel = 'log2CPM'
  if (norma){
    expressionG[,1:length(gene)] = 2^expressionG[,1:length(gene)] -1 # CPM without log2
    expressionG[,1:length(gene)] = expressionG[,1:length(gene)]/t(matrix(colMeans(expressionG[expressionG$Genotype == 'WT',1:length(gene)]), nrow=length(gene), ncol=6))
    ylabel = 'FC'
  }
  
  a = melt(expressionG,id.vars=c("Genotype", "sample"))
  
  bars = ggplot(a, aes(x=variable, y=value, fill=Genotype, color=Genotype)) +
          geom_bar(stat='Summary', fun = 'mean', color="black", position=position_dodge(width = 0.8), width=0.75, size=1) +
          #scale_fill_discrete(name="Genotype", breaks=c(1, 2), labels=c("WT", "KO"))+
          geom_errorbar(stat = "Summary", fun.data = mean_se, position=position_dodge(width = 0.8), width=0.3) +
          geom_point( position =position_dodge(width = 0.8),  size=2.5) +
          xlab("") + ylab(ylabel) +
          stat_compare_means(method = "t.test", label.y = max(a$value)+0.1, label = "p.format") + 
          scale_y_continuous(limits = c(0,max(a$value)+0.5), expand = c(0, 0)) +
          scale_color_manual(values=c(RAGcol, WTcol)) +
          scale_fill_manual(values=c('white', 'white')) +
          theme_classic() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black", angle=90))
  
  
  
  #pvalsdf <- data.frame(matrix(0, ncol=5, nrow = 3))
  #colnames(pvalsdf) <- c('group1', 'group2', 'p', 'p.adj', 'y.position')
  #pvalsdf[1, c('group1', 'group2')] = unique(expressionG$groupWtest)[c(2*i-1,2*i)]
  #pvalsdf[i, 'p'] <- wilcox.test(normExp ~ Genotype, data=expressionG)$p.value

  #pvalsdf[, 'p.adj'] <- round(p.adjust(pvalsdf$p, method='fdr'), 4) 
  
  #pvalsdf$y.position <- 1.2
      #stat_pvalue_manual(pvalsdf, label = 'p.adj', size = 4) +

  return(bars)
}

nameG <- read.csv(paste0(path, 'names.csv')) # the genes we highlighted
top <- merge(top, nameG[c('genes', 'name.in.the.paper')], by = 'genes', all.x = TRUE)
top$FinalSymbol <- convertLoc2symbol(top$genes)
write.csv(top, paste0(path, 'topDEG.csv'))

upG = c('LOC107379421', 'LOC107397014', 'spry2', 'spry2', 'spry2')  #top[top$FDR < 0.05 & top$logFC > 1,]$genes   #
downG =  c('LOC107372687', 'LOC107394448', 'LOC107394449', 'LOC107387651', 'LOC107372481') # top[top$FDR < 0.05 & top$logFC < -1,]$genes #

downG = c('LOC107379421','LOC107374380') #'LOC129166898')

upB = barplotExpression(upG)
downB = barplotExpression(downG)

pdf(paste0(path, '/figures/barplot3.pdf'), width=5, height=7) 
print(plot_grid(upB, downB, ncol = 1,align = "v"))
#print(downB)
dev.off()

pdf(paste0(path, '/figures/barplot2up.pdf'), width=4, height=4) 
#print(plot_grid(upB, downB, ncol = 1,align = "v"))
print(upB)
dev.off()

#volcano
top$DE = 'NO'

top[top$FDR < 0.05 & top$logFC > 1,]$DE = 'up'
top[top$FDR < 0.05 & top$logFC < -1,]$DE = 'down'
top[!top$name.in.the.paper %in% NA,]$DE = 'Immune'

top$labels <- NA
top$labels = top$name.in.the.paper

h = ggplot(data=top, aes(x=logFC, y=-log10(FDR), col=DE, label=labels)) + 
      geom_point() + 
      geom_text_repel( max.overlaps = Inf) +
      #geom_text() +
      theme_classic() +
      scale_color_manual(values=c(RAGcol,  "gray", WTcol)) + #'red',
      geom_hline(yintercept=-log10(0.05), col="gray")+
      geom_vline(xintercept=c(-1, 1), col="gray")
      

pdf(paste0(path, '/figures/volcano2.pdf'), width=6, height=4) 
print(h)
dev.off()

write.csv(top[top$DE != 'NO',], paste0(path, 'DEgenesH.csv'))

# GO visualization
GOdots <- function(GOdata, chosenPathways, title="GO"){
  #order the pathways according FDR or the original order
  #  pathwaysOrder <- chosenPathways$Description
  #
  
  GOdata <- separate(data = GOdata[chosenPathways,], col = GeneRatio, into = c("numerator", "denominator"), sep = "\\/")
  GOdata$numerator = as.numeric(GOdata$numerator)
  GOdata$denominator = as.numeric(GOdata$denominator)
  GOdata$GeneRatio = GOdata$numerator / GOdata$denominator
  GOdata$qvalue = as.numeric(GOdata$qvalue)
  
  GOdata$Description <- paste0(toupper(substr(GOdata$Description, 1, 1)), substr(GOdata$Description, 2, nchar(GOdata$Description)))
  #pathwaysOrder <- paste0(toupper(substr(pathwaysOrder, 1, 1)), substr(pathwaysOrder, 2, nchar(pathwaysOrder)))
  GOdata$Description <- factor(GOdata$Description, levels = rev(GOdata$Description))
  
  dotplotGO <- ggplot(data = GOdata, aes(x = GeneRatio, y = Description, color= qvalue)) + 
    geom_point(aes( size = Count)) + 
    scale_size(range = c(min(GOdata$Count) +0.3, 6)) +
    #  scale_color_gradient(low = "blue", high = "red") +
    scale_color_gradient(low = "red", high = "blue") +
    #scale_shape_manual(values = c(19, 21) , breaks = waiver())+, limits = c(0,0.05)
    theme_bw() +
    ylab("") + 
    ggtitle(title)
  
  return(dotplotGO)
}

chosenPathways = read.csv(paste0(path, 'GO selected pathways.csv'))
viUp = GOdots(up, chosenPathways[chosenPathways$Direction == 'up',]$ID, 'Up')
viDn = GOdots(down, chosenPathways[chosenPathways$Direction == 'down',]$ID, 'Down')

pdf(paste0(path, '/figures/GO.pdf'), width=6, height=6) 
print(plot_grid(viUp, viDn, ncol = 1,align = "v"))
dev.off()

