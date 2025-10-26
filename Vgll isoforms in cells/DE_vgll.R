library(edgeR)            #3.34.1
library(ggfortify)        #0.4.13
library(ggplot2)          #3.3.5
library(stats)            #4.1.1
library(factoextra)       #1.0.7
library(ggConvexHull)     #0.1.0
library(gridExtra)
library(ComplexHeatmap)   #2.8.0
library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0
library(xlsx)             #0.6.5
library(VennDiagram)      #1.7.0
library(reshape2)         #1.4.4
library(sva)
library(tidyr)            #1.1.4
library(stringr)
library(scales)
library(ggrepel)          #0.9.1


# sources----
path <- 'C:/Users/tehil/Dropbox/Projects/VGLL/'
expressionPath <- 'inputFiles/matrix_expression_RNAseq_star2504_VGLL.csv' 
sampleInfoPath <- 'inputFiles/infoTable VGLL.csv'

humanKfConversion = read.table(paste0(path, "inputFiles/NCBI-Human-orthologs.txt"), head = T, sep = "\t")

#function
createDEGobject <- function(path, expPath, infoPath){
  # take two table- 1. RNA seq with genes on the columns and samples in row 2. information table- age, genotype, and sex for each sample
  # create DEG object, order the samples(first - young WT male full)
  # return the DE-object
  
  # read the expression (counts) table and order the names of genes
  exMatrixT <- read.csv(paste0(path, expPath), row.names = 1)
  exMatrix <- data.frame(t(exMatrixT))
  colnames(exMatrix) <- row.names(exMatrixT)
  
  # read the samples information table, make the parameter as factor and re-factor them(wt, male, full, young as default)
  infoTable <- read.csv(paste0(path, infoPath), row.names = 1, stringsAsFactors=T)
  infoTable$genotype <- factor(infoTable$genotype, levels = c('WT', 'KD', 'null'))
  infoTable$mutation <- factor(infoTable$mutation, levels = c('WT', 'Ex1', 'Ex3'))
  #infoTable$batch <- factor(infoTable$batch)
  #infoTable$tissue <- relevel(infoTable$tissue, ref = "Brain")
  infoTable$sex <- relevel(infoTable$sex, ref = "male")
  infoTable$sample <- row.names(infoTable)
  
  infoTable <- infoTable[colnames(exMatrix),] # very important line!! order the info according the RNA-seq data; the first column in count table will be corresponded to the first row in sample table
  
  # create DGE object     infoTable$sex,
  DEobj <- DGEList(exMatrix, samples = infoTable, group=paste(infoTable$mutation),
                   genes = row.names(exMatrix))
  
  #DEobj <- DEobj[,order(DEobj$samples$genotype)]
  #DEobj <- DEobj[,order(DEobj$samples$tissue)]
  #DEobj <- DEobj[,order(DEobj$samples$feed)]
  #DEobj <- DEobj[,order(DEobj$samples$sex)]
  DEobj$samples$group <- factor(DEobj$samples$group, levels = unique(DEobj$samples$group))
  return(DEobj)
}
filterNormDEGobject <-function(DEobj){
  # filtering low expressed genes
  # normalized each sample by weight average of non-DE genes using TMM
  
  # filtering out lower expressed genes
  keep <- filterByExpr(DEobj, model.matrix(~mutation, data=DEobj$samples))# 
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
dataExplorer <- function(DEobj, titleG='', logTag = T, clusterSamples = F, colorG = 'group'){
  # dataExplorer accept DGE object
  # return PCA (PC1/2 and PC1/3) with and without labels, and figure of variances against the number of dimensions
  
  #using log2(CPM)
  cpm.data <- log2(cpm(DEobj, normalized.lib.sizes = T)+1)
  cpm.data <- cpm.data[apply(cpm.data, 1, sd) != 0,]
  infoTable <- DEobj$samples  #
  
  #*****************PCA************
  #https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
  pca.de <- prcomp(t(cpm.data), scale. = T)
  fz <-fviz_eig(pca.de)
  # PCA with labels PC1/2 and PC1/3  group batch
  pcaL12 <- autoplot(pca.de, data = infoTable, colour = colorG, x=1, y=2,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable[,colorG], color = infoTable[,colorG]), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaL13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  if (length(unique(infoTable$group))==2)
    colors_pca = c("#C0C0FF", "#0000C0") #"#DDA94B", "#bc5b3c", "#e49c66"
  else
    colors_pca = c("grey50", "#76069A", "#db81f7") #c("#22658B", "#d08435", "#B26EAB")       # "#009efa", , "#f4cac5")
  
  #PCA with dots
  pcaD12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=0.2),
                    show.legend = T)+
    scale_color_manual(name='Group', values = colors_pca)+
    scale_fill_manual(name='Group', values = colors_pca)+
    #scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaD13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=0.2),
                    show.legend = T)+
    scale_color_manual(name='Group', values = colors_pca)+
    scale_fill_manual(name='Group', values =colors_pca)+
    #scale_alpha_manual(name='Group', values = c(0.2,0.2,0,0))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)

  return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13, PCAdots12 = pcaD12, PCAdots13 = pcaD13, jjj=pca.de))
  #return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13))
}
actualHeatmap <- function(cpm.data, infoData, geneList, goID, average=F, convertNames = T, splitCol = 'none', clusterRows=T, clusterSamples = F, geneOrder=NA, scaleT = T,
                          filterGenes=F, fdr=NA, fdrMfc=NA, cutoff=0.05, hirerByAnoutherDS = F, Adataset = c(), infoDS = 'none', ...){
  # cpm.data-   dataframe with expression;    
  # infoData-   sample information: age, genotype, sex, feeding condition
  # geneList-   the genes to plot in the heatmap
  # goID-       name of the GO pathway
  # average-    use the average of each group. default False
  # splitCol-   split the columns by category. could be: 'none'(default), 'age', 'genotype', 'feed', 'sex' and any new col-name you create
  # clusterRows, clusterSamples- cluster the rows or the columns by hierarchical clustering. default- T for the rows, F for the columns
  # geneOrder-  ordering the genes by a list of genes
  # scaleT-     scale each gene by row. Boolean value; default True
  
  # filterGenes-filter the genes by FDR, FC or FDR*log(FC); Boolean value, default False
  # if filterGenes is True:
  # fdr or fdrMfc - should contain dataframe of the FDR/ FDR*log(FC) values. if not will use the FC dataframe
  # cutoff-     the threshold for the filterGenes; default- 0.05
  # 
  # hirerByAnoutherDS- use anouther dataset to order the genes by hierarchical clustering. Boolean value, default False
  # if True- Adataset should contain the dataframe and infoDS the information about the samples in the new dataset.
  # ... others parameters for the Heatmap function
  
  # took the genes are exits in the CMP data
  geneList <- geneList[geneList %in% rownames(cpm.data)]
  
  if (hirerByAnoutherDS){
    geneList <- geneList[geneList %in% rownames(Adataset)  ]
  }
  cpm.data <- cpm.data[geneList,]
  
  if (average){
    av <- averageSamples(cpm.data, infoData)
    cpm.data <- av$cpm.mean
    infoData <- av$infoTable
  }
  
  if (filterGenes){ # at least one of the column are passed the cutoff
    if (!is.na(fdr)){
      fdr <- fdr[geneList,]
      cpm.data <- cpm.data[apply(fdr,1,function(x) sum(x < cutoff) != 0),]
    }
    else if (!is.na(fdrMfc)){
      fdrMfc <- fdrMfc[geneList,]
      cpm.data <- cpm.data[apply(fdrMfc,1,function(x) sum(abs(x) > cutoff) != 0),]
    }
    else{
      cpm.data <- cpm.data[apply(cpm.data,1,function(x) sum(abs(x) > cutoff) != 0),]
    }
  }
  
  if (nrow(cpm.data) < 2) # if the matrix too small- draw a random heatmap 
    return(Heatmap(matrix = c(1,2,3,4, name='ERROR!!')))
  
  cpm.data <- na.omit(cpm.data)
  
  if (!clusterRows){
    #head(cpm.data)
    geneOrder <- geneOrder[geneOrder %in% row.names(cpm.data)]
    cpm.data <- cpm.data[geneOrder,]
    hr = F
  }
  else
    hr <- hclust(as.dist(1-cor(t(cpm.data), method="pearson"))) # Cluster rows by Pearson correlation.
  #hc <- hclust(as.dist(1-cor(cpm.data, method="pearson"))) # Clusters columns by Pearson correlation.
  
  if (scaleT)
    data.scaled <- t(scale(t(cpm.data))) 
  else{
    data.scaled <- as.matrix(cpm.data)
  }
  #print(data.scaled)
  # used hierarchical clustering from another dataset
  if(hirerByAnoutherDS){
    geneList <- geneList[geneList %in% rownames(Adataset)]
    geneList <- geneList[geneList %in% row.names(cpm.data)]
    Adataset <- Adataset[geneList,]
    
    if (average){
      av <- averageSamples(Adataset, infoDS)
      Adataset <- av$cpm.mean
    }
    hr <- hclust(as.dist(1-cor(t(Adataset), method="pearson"))) # Cluster rows by Pearson correlation.
  }
  
  if (convertNames)
    rownames(data.scaled) <- convertLoc2symbol(rownames(data.scaled))
  c.annotaion <- HeatmapAnnotation(mutation=infoData$mutation, genotype=infoData$genotype, tissue=infoData$tissue,
                                   age=infoData$age, sex=infoData$sex, feed=infoData$feed,
                                   col = list(mutation=c('WT'='black', 'Ex1'='gray60', 'Ex3'='gray'),
                                              genotype=c('WT'= 'blue','KD'='darkorchid2', 'null'='red'),
                                              sex=c('male' = '#1F8374', 'female'= '#5F2374')),
                                   simple_anno_size = unit(2, "mm"))
  if (splitCol == 'none')
    infoData['none'] = rep("", nrow(infoData))
  
  
  hm <- Heatmap(data.scaled, name = goID, top_annotation = c.annotaion, cluster_columns = clusterSamples,
                column_split  = infoData[splitCol], cluster_rows = hr, row_title_rot =0,
                row_title_gp = gpar(fontsize =8), row_names_gp = gpar(fontsize = 5), ...) 
  return(hm)
}
convertLoc2symbol <- function(oldNames){
  # convert symbol to the final symbol.
  NCBIgeneNames <- read.csv("C:/Users/tehil/Dropbox/Projects/genes_names4.csv")
  names(NCBIgeneNames) <- c('NCBI', 'FinalSymbol', 'Human')
  head(NCBIgeneNames)
  newNames <- NCBIgeneNames[NCBIgeneNames$NCBI %in% oldNames,]
  newNames <- newNames[match(oldNames, newNames$NCBI),]$FinalSymbol
  return(newNames)
}
convertsymbol2Loc <- function(oldNames){
  # convert symbol to the final symbol.
  NCBIgeneNames <- read.csv("C:/Users/tehil/Dropbox/Projects/genes_names4.csv")
  names(NCBIgeneNames) <- c('NCBI', 'FinalSymbol', 'Human')
  head(NCBIgeneNames)
  newNames <- NCBIgeneNames[NCBIgeneNames$FinalSymbol %in% oldNames,]
  newNames <- newNames[match(oldNames, newNames$FinalSymbol),]$NCBI
  return(newNames)
}
convertLoc2human <- function(oldNames){
  # convert symbol to the final symbol.
  head(humanKfConversion)
  newNames <- humanKfConversion[humanKfConversion$ncbi %in% oldNames,]
  newNames <- newNames[match(oldNames, newNames$ncbi),]$human
  return(newNames)
}


# create DGE object and separate the brains and livers.
DEobj <- createDEGobject(path, expressionPath, sampleInfoPath)

#combined technical replicates
bio_reps <- list(
  Vgll3_KD_3 = c("Vgll3_KD_3B"),  #"Vgll3_KD_3A", 
  Vgll3_KD_4 = c("Vgll3_KD_4A", "Vgll3_KD_4B"),
  Vgll3_KD_6 = c("Vgll3_KD_6A", "Vgll3_KD_6B"),
  Vgll3_null_1 = c("Vgll3_null_1A", "Vgll3_null_1B"),
  Vgll3_null_3 = c("Vgll3_null_3A", "Vgll3_null_3B"),
  Vgll3_null_4 = c("Vgll3_null_4A", "Vgll3_null_4B"),
  WT_F4 = c("WT_F4A", "WT_F4B"),
  WT_FM = c("WT_FMA", "WT_FMB"), 
  WT_F1 = c("WT_F1", "WT_F9"))

avg_counts <- sapply(names(bio_reps), function(name) {
  rowMeans(DEobj$counts[, bio_reps[[name]], drop = FALSE])
})
new_samples <- do.call(rbind, lapply(bio_reps, function(reps) {
  # Take the first one as a representative row
  rep_row <- DEobj$samples[reps[1], ]
  rownames(rep_row) <- NULL
  rep_row
}))
new_samples$sample = rownames(new_samples)

DEobj_avg <- DGEList(counts = avg_counts, samples = new_samples[, c('genotype', 'mutation', 'sex')])
DEobj_avg$samples$group = DEobj_avg$samples$mutation

DEobj <- filterNormDEGobject(DEobj)
DEobj_avg <- filterNormDEGobject(DEobj_avg)

dev.off()

vgllGraphs <- dataExplorer(DEobj, 'VGLL')#)
vgllGraphsAvg <- dataExplorer(DEobj_avg, 'VGLL')#)

cpmVgll <- log2(cpm(DEobj_avg, normalized.lib.sizes = T)+1)

#heatmaps
pdf(paste0(path, '/figures/heatmaps vgll unsupervised.pdf'), width=12, height=18) 
vgllHeatmap <- actualHeatmap(cpmVgll, DEobj$samples, rownames(DEobj), clusterSamples = T, 'Vgll') 
draw(vgllHeatmap, column_title = "vgll heatmap")
rm(vgllHeatmap)
dev.off()

#save the PCA graphs
pdf(paste0(path, '/figures/PCA VGLL avg.pdf'), width=6, height=8) 
grid.arrange(vgllGraphsAvg$PCAlabels12, vgllGraphsAvg$PCAlabels13,
             vgllGraphsAvg$PCAdots12,   vgllGraphsAvg$PCAdots13,
             vgllGraphsAvg$knee, nrow = 3)
dev.off()


# differential expression----
humanKfConversion = read.table("C:/Users/tehil/Dropbox/Projects/NCBI-Human-orthologs.txt", head = T, sep = "\t")

DErankingClasicSaving <- function(DEobj, pairs, test, pathGL, saveRanking=F){
  DEobj <- estimateDisp(DEobj)
  
  FC <- data.frame(row.names = row.names(DEobj))
  pval <- data.frame(row.names = row.names(DEobj))
  FDR <- data.frame(row.names = row.names(DEobj))
  graph <- list()
  
  for (p in 1:length(pairs)){
    et <- exactTest(DEobj, pair = pairs[[p]])
    top <- topTags(et, n=Inf)$table
    print(names(pairs[p]))
    print(paste('+', dim(top[top$FDR < 0.05 & top$logFC > 0,])))
    print(paste('-', dim(top[top$FDR < 0.05 & top$logFC < 0,])))
    print(head(top[top$FDR < 0.05 & top$logFC < 0,]))
    top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
    
    top2 = merge(top, geneLenM, by.x=0, by.y='Geneid')
    
    top2$group = with(top2, ifelse(top2$FDR < 0.01 & top2$logFC > 0, 'up',
                                   ifelse(top2$FDR < 0.01 & top2$logFC < 0, 'down', 'NS')))
    print(table(top2$group))
    print(head(top2))

    FC[[names(pairs)[p]]] <- top[rownames(FC),]$logFC
    pval[[names(pairs)[p]]] <- top[rownames(pval),]$PValue
    FDR[[names(pairs)[p]]] <- top[rownames(FDR),]$FDR
    
    #saving the ranking list
    ranking_genes <- data.frame(Gene = row.names(top), symbol= convertLoc2symbol(row.names(top)), human= convertLoc2human(row.names(top)), mlog10QvalxFC = top$mlog10QvalxFC, FC=top$logFC, pval=top$PValue, FDR=top$FDR)
    if (saveRanking)
      write.csv(ranking_genes, paste0(pathGL, 'ranking_', names(pairs[p]), test, '.csv'), row.names = F)
  }
  return(list(fc=FC, fdr=FDR, pval=pval, graph=graph, ranking_genes=top2))
}

DErankingSaving <- function(DEobj, design, con, test, pathGL, saveRanking=F){
  DEobj <- estimateDisp(DEobj, design)
  fit <- glmQLFit(DEobj, design, robust = T)
  FC <- data.frame(row.names = row.names(DEobj))
  pval <- data.frame(row.names = row.names(DEobj))
  FDR <- data.frame(row.names = row.names(DEobj))
  graph <- list()
  for (g in 1:length(con)){
    qlf <- glmQLFTest(fit, contrast = con[[g]])
    top <- topTags(qlf, n=Inf)$table
    print(names(con[g]))
    print(paste('+', dim(top[top$FDR < 0.05 & top$logFC > 0,])))
    print(paste('-', dim(top[top$FDR < 0.05 & top$logFC < 0,])))
    print(head(top[top$FDR < 0.05 & top$logFC < 0,]))
    top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
    
    top2 = merge(top, geneLen, by.x=0, by.y='Geneid')
    
    top2$group = with(top2, ifelse(top2$FDR < 0.05 & top2$logFC > 0, 'up',
                            ifelse(top2$FDR < 0.05 & top2$logFC < 0, 'down', 'NS')))
    print(table(top2$group))
    ttest_name = ifelse((table(top2$group)['up'] > 1) & (table(top2$group)['down'] > 1), 
                        paste('up:' ,round(t.test(top2[top2$group == 'up',]$Length, top2[top2$group == 'NS',]$Length)$p.value,6),
                              'down:' ,round(t.test(top2[top2$group == 'down',]$Length, top2[top2$group == 'NS',]$Length)$p.value,6),
                              '\ngenes up:', table(top2$group)['up'], 'NS:', table(top2$group)['NS'], 'down:', table(top2$group)['down']), '--')
    
    FC[[names(con)[g]]] <- top[rownames(FC),]$logFC
    pval[[names(con)[g]]] <- top[rownames(pval),]$PValue
    FDR[[names(con)[g]]] <- top[rownames(FDR),]$FDR
    
    #saving the ranking list
    ranking_genes <- data.frame(Gene = row.names(top), symbol= convertLoc2symbol(row.names(top)), human= convertLoc2human(row.names(top)), mlog10QvalxFC = top$mlog10QvalxFC, FC=top$logFC, pval=top$PValue, FDR=top$FDR)
    if (saveRanking)
      write.csv(ranking_genes, paste0(pathGL, 'ranking_', names(con[g]), test, '.csv'), row.names = F)
  }
  return(list(fc=FC, fdr=FDR, pval=pval, graph=graph))
}

geneLenW = read.csv("C:/Users/tehil/Dropbox/Projects/geneLengthWholeGene.csv", row.names=1)
geneLen = read.csv("C:/Users/tehil/Dropbox/Projects/geneLength.csv", row.names=1)
geneLenM = merge(geneLen, geneLenW, by='Geneid')

con <- list(WTvsEx1 = c(-1,1,0),
            WTvsEx3 = c(-1,0,1))

pairDE <- list(WTvsEx1 = c('WT', 'Ex1'),
               WTvsEx3 = c('WT', 'Ex3'),
               Ex1vsEx3 = c('Ex1', 'Ex3'))


tests <- c('null', 'null', 'null', 'dnull', 'dnull')

experimentGroups <- list(vgll = DEobj_avg) 

conIndex = list(classic = c(1, 2,3), ex1 = c(1), ex3 = c(2)) #, linear = c(1,6,7,8))
modeDE = c('classic','glm', 'glm','classic')

FCs <- list()
pvals <- list()
FDRs <- list()
graphs <- list()

for (i in 1:1){ #length(experimentGroups)
  y <- experimentGroups[[i]]
  if (modeDE[i] == 'classic') {
    pairs_new <- pairDE[conIndex[[i]]]
    names(pairs_new) <- paste(names(experimentGroups)[i], names(pairs_new), sep = '_')
    
    de <- DErankingClasicSaving(y, pairs_new, '', paste0(path, 'GO/GeneSets/'), F)
  } else {
    design <- model.matrix(~0+group, data=y$samples)  #
    colnames(design) <- gsub('group', '', colnames(design))
    print(colnames(design))
    
    contrast_new <- con[conIndex[[i]]]
    names(contrast_new) <- paste(names(experimentGroups)[i], names(contrast_new), sep = '_')
  
    de <- DErankingSaving(y, design, contrast_new, '', paste0(path, 'GO/GeneSets/'), F)
  }
   
  FCs[[names(experimentGroups)[i]]] = de$fc
  pvals[[names(experimentGroups)[i]]] = de$pval
  FDRs[[names(experimentGroups)[i]]] = de$fdr
  graphs = append(graphs, de$graph)
}


# linear
y <- DEobj_avg
continuousV = c(0,0,0,2,2,2,1,1,1) #dKO
groupName = 'VGLL_linear'
design <- model.matrix(~continuousV, data=y$samples) 

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = T)
for (cof in 2){
  qlf <- glmQLFTest(fit, coef=cof)
  top <- topTags(qlf, n=Inf)$table
  print(paste('+', dim(top[top$FDR < 0.05 & top$logFC > 0,])))
  print(paste('-', dim(top[top$FDR < 0.05 & top$logFC < 0,])))
  print(head(top[top$FDR < 0.05,]))
  top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
  ranking_genes <- data.frame(Gene = row.names(top),symbol= convertLoc2symbol(row.names(top)), human= convertLoc2human(row.names(top)), mlog10QvalxFC = top$mlog10QvalxFC, FC=top$logFC, pval=top$PValue, FDR=top$FDR)
  write.csv(ranking_genes, paste0(path, 'GO/GeneSets/', 'ranking_', groupName, colnames(design)[cof], '.csv'), row.names = F)
}
FCs[[groupName]][['continuous']] <- top[rownames(FCs[[groupName]]),]$logFC
FDRs[[groupName]][['continuous']] <- top[rownames(FDRs[[groupName]]),]$FDR

head(FDRs$LIVER_ATMcGAS[order(FDRs$LIVER_ATMcGAS$continuous),])
head(FDRs$BRAIN_ATMcGAS[order(FDRs$BRAIN_ATMcGAS$continuous),])

FCall <- merge(FCs$LIVER, FCs$BRAIN, by=0)
rownames(FCall) = FCall$Row.names
FCall = FCall[2:length(FCall)]

pvalall <- merge(pvals$LIVER, pvals$BRAIN, by=0)
rownames(pvalall) = pvalall$Row.names
pvalall = pvalall[2:length(pvalall)]

FDRall <- merge(FDRs$LIVER, FDRs$BRAIN, by=0)
rownames(FDRall) = FDRall$Row.names
FDRall = FDRall[2:length(FDRall)]

FDRall$FinalSymbol = convertLoc2symbol(rownames(FDRall))

write.csv(FCall, paste0(path, 'outputFiles/differential analysis_FC.csv'))
write.csv(pvalall, paste0(path, 'outputFiles/differential analysis_Pval.csv'))
write.csv(FDRall, paste0(path, 'outputFiles/differential analysis_FDR.csv'))


# volcano ####

Ex1 = data.frame(Gene = rownames(FDRs$vgll), symbol = convertLoc2symbol(rownames(FDRs$vgll)), FDR = FDRs$vgll$vgll_WTvsEx1, FC =FCs$vgll$vgll_WTvsEx1)

Ex1$DE = 'NO'
Ex1[Ex1$FDR < 0.05 & Ex1$FC > 1,]$DE = 'up'
Ex1[Ex1$FDR < 0.05 & Ex1$FC < -1,]$DE = 'down'
#top[!top$name.in.the.paper %in% NA,]$DE = 'Immune'

Ex1$labels <- NA
Ex1[!Ex1$DE == 'NO',]$labels = Ex1[!Ex1$DE == 'NO',]$symbol

h = ggplot(data=Ex1, aes(x=FC, y=-log10(FDR), col=DE)) + #, label=labels
  geom_point() + 
  #geom_text_repel(max.overlaps = Inf) +
  #geom_text() +
  theme_classic() +
  scale_color_manual(values=c('red', "gray", 'blue')) +
  geom_hline(yintercept=-log10(0.05), col="gray")+
  geom_vline(xintercept=c(-1, 1), col="gray")


pdf(paste0(path, '/figures/volcano.pdf'), width=6, height=4) 
print(h)
dev.off()

#venn
genes <- list(Ex1 = row.names(FDRs$vgll[FDRs$vgll$vgll_WTvsEx1 < 0.05, ]), Ex3 = row.names(FDRs$vgll[FDRs$vgll$vgll_WTvsEx3 < 0.05, ])) 
v <- venn.diagram(genes, category.names = names(genes), main= "venn vgll", filename = NULL)
ggsave(v, file=paste0(path, '/figures/venn.pdf'), device = "pdf", width=5.5, height = 6)


#GSEA----
library(msigdbr)
subC = c('CP:KEGG', 'CP:BIOCARTA', 'CP:REACTOME', 'GO:BP', 'GO:MF')#)
all_gene_sets = msigdbr(species = "human")
msigdbr_t2g  = all_gene_sets[all_gene_sets$gs_subcat %in% subC, c('gs_name', 'gene_symbol')]
CGP = all_gene_sets[all_gene_sets$gs_subcat %in% c('CGP'), c('gs_name', 'gene_symbol')]

GSEAparam <- function(data, savePath, path){ # thanks for Param Priya Singh
  # This input to GSEA is a ranked list of genes. Read the input ranked list.
  
  # Get human ortholog symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  dataH = merge(humanKfConversion, data, by.x = "ncbi", by.y = "Gene") 
  entrezIds = bitr(as.character(dataH[,2]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Get entrez ids for annotation
  dataHE = merge(dataH, entrezIds, by.x = "human", by.y = "SYMBOL") # Get human symbols
  head(dataHE)
  
  # There can be duplicate values because of paralogs, I take average of those for quantitative score
  unique = aggregate(dataHE[,3], list(dataHE$human), mean)
  dataHEU = merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL") 
  colnames(dataHEU) = c("human", "mlog10QvalxFC", "entrez")
  head(dataHEU)
  
  geneList = dataHEU[,2]  # gene list for GO 
  names(geneList) = as.character(dataHEU[,1]) # with entrez ids as names
  
  geneListKegg = geneList # gene list for KEGG
  names(geneListKegg) = as.character(dataHEU[,3]) #  with humna symbols as names
  
  # *** Sort the gene list based on quantitative score in decreasing order. This is critical for GSEA  
  geneList = sort(geneList, decreasing = TRUE)
  geneListKegg = sort(geneListKegg, decreasing = TRUE)
  
  print(head(dataHEU))
  print(head(geneList))
  print(tail(geneList))
  
  head(geneListKegg)
  tail(geneListKegg)
  
  # *****************  Now do different enrichment analyses *****************************
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                keyType      = 'SYMBOL',
                ont          = c("ALL"),
                pvalueCutoff = 1)
  
  mgsig <- GSEA(geneList     = geneList,
                TERM2GENE = msigdbr_t2g, 
                pvalueCutoff = 1)
  print(paste0(path, 'GO/Results/', savePath, "_GOGSEA.csv"))

  write.table(ego3, paste0(path, 'GO/Results/', savePath, "_GOGSEA.csv"), sep = ",", quote = T, row.names = F)
  write.table(mgsig, paste0(path, 'GO/Results/', savePath, "_ALL.csv"), sep = ",", quote = T, row.names = F)
}
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
  
  kegg <- enricher(gene        = geneListHuman,
                   universe     = backgroundGeneList,
                   TERM2GENE = msigdbr_t2g,
                   pvalueCutoff = 1)
  
  print(length(ego3))
  if (length(ego3) > 0)
    if (nrow(ego3) > 0)
      write.table(ego3, paste0(folder, fileName, "_GO.csv"), sep = ",", quote = T, row.names = F)
  if (length(kegg) > 0)
    if (nrow(kegg) > 0)
      write.table(kegg, paste0(folder, fileName, "_ALL.csv"), sep = ",", quote = T, row.names = F)
  
}

# run GSEA on all the ranking lists from the previous step
files <- list.files(path=paste0(path, 'GO/GeneSets/'), pattern="*.csv")
for (i in 1:length(files)){ 
  f <- gsub('.csv|ranking_| ', '', files[i])
  data = read.csv(paste0(path, 'GO/GeneSets/', files[i]))
  ego3 <- GSEAparam(data[,c('Gene', 'mlog10QvalxFC')], f, path)
}

sigGenesH = humanKfConversion[humanKfConversion$ncbi %in% sigGenesDU, ]$human
background = humanKfConversion[humanKfConversion$ncbi %in% rownames(FCs$female_LIVER_ATMcGAS), ]$human



for (i in 2:length(files)){
  f <- gsub('.csv|ranking_| ', '', files[i])
  data = read.csv(paste0(path, 'GO/GeneSets/', files[i]))
  
  sigGenesHPos = humanKfConversion[humanKfConversion$ncbi %in% data[(data$FDR < 0.05) & (data$FC > 0), ]$Gene, ]$human
  sigGenesHNeg = humanKfConversion[humanKfConversion$ncbi %in% data[(data$FDR < 0.05) & (data$FC < 0), ]$Gene, ]$human
  sigGenesH = humanKfConversion[humanKfConversion$ncbi %in% data[data$FDR < 0.05, ]$Gene, ]$human
  
  
  background = humanKfConversion[humanKfConversion$ncbi %in% data$Gene, ]$human
  
  enrichmentClusterProfile(unique(sigGenesHPos), paste0(path, 'GO/'), paste0('GO ', f, ' classic positive'), background)
  enrichmentClusterProfile(unique(sigGenesHNeg), paste0(path, 'GO/'), paste0('GO ', f, ' classic negative'), background)
  enrichmentClusterProfile(unique(sigGenesH), paste0(path, 'GO/'), paste0('GO ', f, ' classic'), background)
}

# combined GSEA result files ----
ont <- 'BP'
DBs <- 'GO'

extractSigGenes <- function(path, group, ngroup, test, ont= '', pvalCutoff=0.05) {
  # extract  the list of significant enrichment pathways and table enrichment pathways from certain database, test, ontology option (BP, CC or MF) and condition.
  #  there is option to choose p-value cutoff
  
  #file_name <-  paste0(path, 'GO/Results/', group[ngroup], test, '_GOGSEA.csv')
  #file_name <-  paste0(path, 'GO/Results/', group[ngroup], test, '_ALL.csv')
  #GO
  file_name <-  paste0(path, 'GO/', group[ngroup], '_ALL.csv')
  
  print(file_name)
  ExTable <- read.csv(file_name)
  #genes <- ExTable[ExTable$qvalues < pvalCutoff ,]$ID#   & ExTable$ONTOLOGY == ont
  genes <- ExTable[ExTable$qvalue < pvalCutoff ,]$ID#   & ExTable$ONTOLOGY == ont
  return(list(table=ExTable, genes=genes))
}

threeGroupsCompare <- function(path, group, test, DB, nameG, ont= '', pvalCutoff=0.05){
  # create table with all the significant enriched pathway in at least one of the given comparisons,
  # drawing Venn diagram with the number of significant pathways in each group (less than 5 groups, if more- without venn diagram)
  # 
  # return table with the go ID, description and for each group:
  # 1. NES (normalized enrichment score) 2.adjust p.value 3.column that contain the core genes in each pathway
  
  tables <- list()
  genes <- list()
  for (i in 1:length(group)){
    GG <- extractSigGenes(path, group, i, test, ont)
    tables[[i]] <- GG$table
    genes[[i]] <- GG$genes
  }
  
  names(genes) <- names(group)
  #https://www.rdocumentation.org/packages/VennDiagram/versions/1.6.20/topics/venn.diagram
  if (length(genes) < 6){
    v <- venn.diagram(genes, category.names = names(group), main= paste(DB, test, ont),
                      filename = NULL)
    ggsave(v, file=paste0(path, 'GO/Results/MergedAnalysis/', test, '_in_', nameG, '_venn.pdf'), device = "pdf", width=5.5, height = 6)
  }
  # one column contain which group have significant enrichment in the pathway
  vennGroups <- get.venn.partitions(genes)
  vennGroups$..set.. <- gsub("\\).*", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("\\(", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("n", '#', vennGroups$..set..)
  vennGroups$..set.. <- as.character(lapply(vennGroups$..set.., as.name))
  vennGroups$..set.. <- gsub("n", '+', vennGroups$..set..)
  vennGroups$..set.. <- gsub("#", 'n', vennGroups$..set..)
  
  tableGroups <- data.frame(group=rep('NAN', length(unique(unlist(genes)))), ID=unique(unlist(genes)))
  
  for (gr in 1:nrow(vennGroups)){
    if (length(vennGroups$..values..[[gr]]))
      tableGroups[tableGroups$ID %in% vennGroups$..values..[[gr]],]$group <- vennGroups$..set..[gr]
  }
  print(table(tableGroups$group))
  
  #col names for overlapping tables
  #colNamesGroup <- c('ID', 'NES', 'qvalues', 'core_enrichment')
  colNamesGroup <- c('ID', 'GeneRatio', 'BgRatio', 'qvalue', 'geneID')   # 'p.adjust''Count', ## GO
  N <- length(colNamesGroup) - 1
  
  #col names for describe table
  colNamesGeneral <- c('ID', 'Description')#'ONTOLOGY', 
  
  
  tablePvalNES <- merge(tables[[1]][colNamesGroup], tables[[2]][colNamesGroup], all=T, by='ID', suffixes=paste0('.', names(group)[1:2]))
  tablePathwaysNames <- merge(tables[[1]], tables[[2]], by=colNamesGeneral,all=T, suffixes=paste0('.', names(group)[1:2]))
  
  if (length(group) > 2){ # more than 2 groups in the comparisons
    for (k in 3:length(group)){
      tablePvalNES <- merge(tablePvalNES, tables[[k]][colNamesGroup], all=T, by='ID')
      
      colnames(tablePvalNES) <- c(colnames(tablePvalNES)[1:(1+(k-1)*N)], paste0(colnames(tablePvalNES)[(2+(k-1)*N):(1+k*N)], '.', names(group[k])))
      tablePathwaysNames <- merge(tablePathwaysNames, tables[[k]], by=colNamesGeneral, all=T)
    }
  }
  
  # merge the column with the groups are significant in each pathways with the pathway id and description
  tableDescriptionValues <- merge(tableGroups, tablePathwaysNames[colNamesGeneral], by='ID', all.x=T)
  
  #merge everything together
  merge(tableDescriptionValues, tablePvalNES, by='ID', all.x=T)
  
}

sendToCompareSaveXlsx <- function(path, group, DBs, nameG, test, ont=''){
  xlsxFileName <- paste0(path, 'GO/Results/MergedAnalysis/', DBs, ont, '_compare_', test, ' in ', nameG, '.xlsx')
  print(xlsxFileName)
  
  write.xlsx(t(data.frame(DB=DBs, ontology=ont, age='positive- decrease with age', genotype='positive-more in WT')),
             file = xlsxFileName, sheetName='summary', append=FALSE)
  
  comparisonsTable <- threeGroupsCompare(path, group, test, DBs, nameG, 'BP')
  if(nrow(comparisonsTable) > 0)
    write.xlsx(comparisonsTable, file = xlsxFileName, sheetName=DBs, append=TRUE)
  
}


sendToCompareSaveXlsx(path, groups$all, DBs, 'Vgll all2', '', ont)

#GO
sendToCompareSaveXlsx(path, groups$all, DBs, 'Vgll GO Ex1 and Ex3 all +kegg reactom', '', ont)


#dotplot----

GOdots <- function(GoCompareEnrichmentGroup, chosenPathways, group, title="GO", labelss=c(), cirles=c()){
  # for GSEA
  # extract the pathways and the group 
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[GoCompareEnrichmentGroup$ID %in% chosenPathways$ID, -grep('core_', colnames(GoCompareEnrichmentGroup))]
  cnames <- c(colnames(GoCompareEnrichmentGroup)[1:5], colnames(GoCompareEnrichmentGroup)[grep(paste(group, collapse = '|'), colnames(GoCompareEnrichmentGroup))])
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[,cnames]
  
  #order the pathways according FDR or the original order
  pathwaysOrder <- chosenPathways$Description
  #
  #select the pathways with significant value at least in one group
  if (length(grep('qvalues', colnames(GoCompareEnrichmentGroup))) > 1){
    checkTable <- GoCompareEnrichmentGroup[, grep('qvalues', colnames(GoCompareEnrichmentGroup))]
    row.names(checkTable)
    checkTable[is.na(checkTable)] = 1
    GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[rowSums(checkTable < 0.05) > 0,]
    #colnames(GoCompareEnrichmentGroup) <- gsub('p.adjust', 'padjust', colnames(GoCompareEnrichmentGroup))
  }
  
  meltEnrich <- melt(GoCompareEnrichmentGroup, id.vars = c('ID', 'Description'))
  meltEnrich <- meltEnrich[!meltEnrich$variable %in% c('ONTOLOGY', 'NA.'),]
  meltEnrich <- separate(data = meltEnrich, col = variable, into = c("measure", "Group"), sep = "\\.")
  meltEnrich$Group <- factor(meltEnrich$Group, levels = unique(meltEnrich$Group))
  meltEnrich$value <- as.numeric(meltEnrich$value)
  enrichment_score <- meltEnrich[meltEnrich$measure == 'NES', c('ID', 'Group', 'Description', 'value')]
  pval <- meltEnrich[meltEnrich$measure == 'qvalues', c('ID', 'Group', 'Description', 'value')]
  colnames(enrichment_score) <- gsub('value', 'enrichment_score', colnames(enrichment_score))
  colnames(pval) <- gsub('value', 'FDR', colnames(pval))
  goVisDF <- merge(pval, enrichment_score, by=c('ID', 'Group', 'Description'))
  
  goVisDF[is.na(goVisDF)] = 0

  
  #goVisDF$p_adjust[goVisDF$p_adjust > 0.05] = 0
  goVisDF$sig <- c(goVisDF$FDR < 0.05)
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==T, 'Sig')
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==F, 'NS')
  goVisDF$sig <- factor(goVisDF$sig, levels = c('Sig', 'NS'))
  
  goVisDF$Description <- paste0(toupper(substr(goVisDF$Description, 1, 1)), tolower(substr(goVisDF$Description, 2, nchar(goVisDF$Description))))
  pathwaysOrder <- paste0(toupper(substr(pathwaysOrder, 1, 1)), tolower(substr(pathwaysOrder, 2, nchar(pathwaysOrder))))
  goVisDF$Description = gsub("_", " ", goVisDF$Description)
  pathwaysOrder = gsub("_", " ", pathwaysOrder)
  goVisDF$Description <- factor(goVisDF$Description, levels = rev(pathwaysOrder))
  
  
  one <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= enrichment_score)) + 
    geom_point(aes(shape = sig, size = -log10(FDR))) +  
    scale_color_gradient(low = "blue", high = "red") + #, limits = c(-0.8,0.8)
    scale_shape_manual(values = c(19, 21), breaks = waiver())+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  if (length(labelss) == 0 | length(cirles) == 0)
    return(one)
  
  two <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= enrichment_score, shape = sig, size = -log10(FDR))) + 
    geom_point() +  
    scale_color_gradient(low = "blue", high = "red", midpoint = 0) +
    scale_size_continuous(name='-log10(FDR)', breaks=labelss, labels=labelss)+
    scale_shape_manual(values = c(19, 21)) +
    guides(size = guide_legend(override.aes = list(shape = cirles)))+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  return(two)
}
GOdotsHG <- function(GOdata, chosenPathways, title="GO"){
  # for GO enrichment using hyper geometric test
  #order the pathways according FDR or the original order
  #  pathwaysOrder <- chosenPathways$Description
  #
  
  GOdata <- separate(data = GOdata[GoCompareEnrichmentGroup$ID %in% chosenPathways$ID,], col = GeneRatio, into = c("numerator", "denominator"), sep = "\\/")
  GOdata$numerator = as.numeric(GOdata$numerator)
  GOdata$denominator = as.numeric(GOdata$denominator)
  GOdata$GeneRatio = GOdata$numerator / GOdata$denominator
  GOdata$qvalue = as.numeric(GOdata$qvalue)
  GOdata$Count = as.numeric(GOdata$Count)
  print(GOdata$Count)
  
  GOdata$Description <- paste0(toupper(substr(GOdata$Description, 1, 1)), tolower(substr(GOdata$Description, 2, nchar(GOdata$Description))))
  #pathwaysOrder <- paste0(toupper(substr(pathwaysOrder, 1, 1)), tolower(substr(pathwaysOrder, 2, nchar(pathwaysOrder))))
  GOdata$Description = gsub("_", " ", GOdata$Description)
  pathwaysOrder = GOdata[order(GOdata$GeneRatio),]$Description
  GOdata$Description <- factor(GOdata$Description, levels = pathwaysOrder) #rev(GOdata$Description)
  GOdata$Count = log10(GOdata$Count)
  
  dotplotGO <- ggplot(data = GOdata, aes(x = GeneRatio, y = Description, color= qvalue)) + 
    geom_point(aes( size = Count)) + 
    scale_color_gradient(low = "red", high = "blue") +
    #scale_shape_manual(values = c(19, 21) , breaks = waiver())+, limits = c(0,0.05)
    theme_bw() +
    ylab("") + 
    ggtitle(title)
  
  return(dotplotGO)
}


goID <- chosenPathways$ID
names(goID) <- chosenPathways$Description
godots <- GOdots(GoCompareEnrichmentGroup2, chosenPathways, c('VGLL_linear_continuousV'), "GO") #, c(1,2,3,4,5,6,7,8), c(21,19,19,19,19,19,19,19))
pdf(paste0(path, 'figures/GSEA linear vgll2.pdf'), width = 5, height = 5)
plot(godots)
dev.off()


#visualization----

DEobj_avg$samples$group = factor(DEobj_avg$samples$group, levels = c('Ex1', 'WT', 'Ex3'))


cpmVgll <- log2(cpm(DEobj_avg, normalized.lib.sizes = T)+1)


barplotExpression4 <- function(expp, infoData, gene, FCtable, FDRtable, normBar = T){
  # presenting all 4 groups
  ppar <- data.frame(Exp = expp, 
                     sample=rownames(infoData), group=infoData$group, 
                     genotype=infoData$genotype, mutation=infoData$mutation)
  
  ppar$newExp <- 2^ppar$Exp -1 #cpm without log2
  aa <- aggregate(ppar$newExp, list(ppar$group), mean)
  aa$x[1:length(unique(infoData$group))] <- aa$x[1] #normalize the males to WT young male.
  ppar = merge(ppar, aa, by.x="group", by.y="Group.1", sort=F)
  ppar$normExp = ppar$newExp / ppar$x
  ppar$gg <- gsub('WT fasted ', '', ppar$group)
  ppar$gg <- factor(ppar$gg, levels = unique(ppar$gg))
  
  'if (length(unique(infoData$group))==2)
    colors_pca = c("#0072b2", "#009e73") #"#DDA94B"
  else
    colors_pca = c("#0072b2", "#d55e00", '#9f4a96')'
  if (length(unique(infoData$group))==2)
    colors_pca = c("#C0C0FF", "#0000C0") #"#DDA94B", "#bc5b3c", "#e49c66"
  else
    colors_pca = c("#EECCF9", "#76069A", "#db81f7") #c("#22658B", "#d08435", "#B26EAB")       # "#009efa", , "#f4cac5")
  
  if (normBar){
    #tukey$y.position <- seq(max(ppar$normExp)-0.25,max(ppar$normExp)+1.5, length.out=6)
    pparPlot <- ggplot(ppar, aes(x=group, y=normExp)) +
      geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
      geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
      geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
      #scale_color_manual(values = c('black', '#008000')) +
      scale_x_discrete(guide = guide_axis(n.dodge = 4)) +
      scale_color_manual(name='group', values = colors_pca)+
      ggtitle(convertLoc2symbol(gene) , subtitle = paste(paste(colnames(FCtable), collapse = ' '),
                                                        paste(round(FCtable[gene,], 3), collapse = ' '),
                                                        paste(round(FDRtable[gene,], 3), collapse = ' '), sep="\n")) +
      xlab("") +
      #scale_y_continuous(expand = c(0, 0)) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))
  }  
  else {
    #tukey$y.position <- seq(max(ppar$Exp)+0.5,max(ppar$Exp)+3.5, length.out=6)
    pparPlot <- ggplot(ppar, aes(x=group, y=Exp)) +
      geom_bar(stat='Summary', fun = 'mean', color="black", fill='white', width=0.75, size=1) +
      geom_errorbar(stat = "Summary", fun.data = mean_se, width=0.3) +
      geom_point(position = position_jitter(width=0.15), aes(color=genotype), size=2.5) +
      #scale_color_manual(values = c('black', '#008000')) +
      scale_x_discrete(guide = guide_axis(n.dodge = 4)) +
      scale_color_manual(name='group', values = colors_pca)+
      ggtitle(convertLoc2symbol(gene) , subtitle = paste(paste(colnames(FCtable), collapse = ' '),
                                                        paste(round(FCtable[gene,], 3), collapse = ' '),
                                                        paste(round(FDRtable[gene,], 3), collapse = ' '), sep="\n")) +
      xlab("") +
      #scale_y_continuous(limits = c(0,max(ppar$Exp)+5), expand = c(0.02, 0)) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(color="black"), axis.text.x = element_text(color="black"))
    
  }
  pparPlot
}

gene = 'vgll3'
barplotExpression4(cpmVgll[gene,], DEobj_avg$samples, gene, FCs$vgll, FDRs$vgll, T)


chosenPathways = read.csv(paste0(path, '/pathways interferons.csv'))  #pathways dotplot.csv

# barplot according to list
geneSetKf = convertsymbol2Loc(read.csv(paste0(path, 'inputFiles/vgll genes.csv'))$FinalSymbol)     #purkinje cell markers human
geneSetKf = humanKfConversion[humanKfConversion$human %in% c('EMX2', 'GTF3A', 'TF3A', 'ZFP541', 'HNF4A', 'HEB', 'TCF12'),]$ncbi

geneSetKf <- geneSetKf[geneSetKf %in% row.names(cpmVgll)]
geneSetKf = c(geneSetKf, 'vgll3')

for (g in geneSetKf){
  print(paste(g, convertLoc2symbol(g), FDRs$LIVER[g,]))
}

bars = list()
for (gene in geneSetKf){ #rownames(FDRbar)
  bars[[gene]] = barplotExpression4(cpmVgll[gene,], DEobj_avg$samples, gene, FCs$vgll, FDRs$vgll, T)
}
pdf(paste0(path, '/figures/barplot vgll homer ex1 wt ex3.pdf'), width = 10, height = min(18, length(geneSetKf)/4*3.5))
if (length(geneSetKf) > 50){
  for (j in 1:(round(length(geneSetKf)/50)+1)){
    do.call("grid.arrange", c(bars[((j-1)*50+1):(j*50)], ncol = 4))
  }
} else
  do.call("grid.arrange", c(bars, ncol = 4))
dev.off()



