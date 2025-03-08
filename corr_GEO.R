
multiCorr_GEO <- function(geneName, tissueID, cutOff, tissueType){
  library(DT)
  #opening clinical table
  clinical_geo <- read_xlsx(path = "GEO_clinical.xlsx")
  clinical_geo <- as.data.frame(clinical_geo)
  NameExpTable <- paste0(tissueID,"-clean-ann.RData")
  #aim folder of expression table
  expWD <- paste0("TNMdata/GEO_Expression-Tables/", NameExpTable)
  #opening expression table
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- clinical_geo[clinical_geo$Tumor == 1 & clinical_geo$Origin == tissueID, 1]
    
  }else if (tissueType == "Normal") {
    usedTissue <- clinical_geo[clinical_geo$Normal == 1 & clinical_geo$Origin == tissueID, 1]
  }else if (tissueType == "Metastatic") {
    usedTissue <- clinical_geo[clinical_geo$Meta == 1 & clinical_geo$Origin == tissueID, 1]
  }
  
  rownames(expTable) <- expTable$gene_sym
  expTable <- expTable[,-c(1,2)]
  expsNeeded <- as.matrix(expTable[,colnames(expTable) %in% usedTissue])
  
  #make the selected gene to a numeric vector and perform correlation test
  selectedGene <- as.vector(unlist(expsNeeded[row.names(expsNeeded) == geneName,]))
  corResult <- apply(expsNeeded, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4,3)])
  
  #store the results as a dataframe and convert the list to a dataframe
  result = as.data.frame(matrix(nrow = nrow(expsNeeded), ncol = 3))
  colnames(result) <- c("Gene name" , "Spearman's Rho", "p-value")
  
  result[, 1] = names(corResult)
  result[, 2] = round(x = sapply(corResult, "[[", 1), digits = 2)
  result[, 3] = format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
  
  #order the results based on R and remove the values below cutoff
  result <- result[result$`Spearman's Rho` >= min(cutOff) & result$`Spearman's Rho` <=max(cutOff),]
  result <- result[order(result$`Spearman's Rho` ,decreasing = T),]
  
  library(ggraph)
  library(tidygraph)
  maxValue <- max(result$`Spearman's Rho`)
  plotTab <- result
  plotTab$From <- geneName
  plotTab <- plotTab[,c(4,1,2,3)]
  colnames(plotTab) <-c("From", "To", "EdgeWeights", "NodeWeights")
  plotTab$NodeWeights <- as.numeric(plotTab$NodeWeights)
  plotTab <-  as_tbl_graph(plotTab, directed = T)

  networkPlot <- ggraph(graph = plotTab,layout = "drl", circular = F) +
    geom_edge_link(aes(edge_width = EdgeWeights, color = EdgeWeights)) +
    scale_edge_width(range = c(0, maxValue*3),guide="none") +
    scale_edge_color_gradient(low = "#1AFF00",high = "#FF0000",
                              name = "Correlation \n coefficient")+
    geom_node_point(size = 8, color = "grey") +
    geom_node_text(aes(label = name),
                   repel = TRUE,
                   point.padding = unit(0.6, "lines"),
                   colour="black") +
    theme(panel.background = element_rect(fill = "white", colour = "white")) +
    theme(legend.position = "right")
  
  networkPlot
  res <- list(result = result, networkPlot = networkPlot)
  return(res)
}


singleCorr_GEO <- function(geneA, geneB, corrMethod, tissueID, tissueType) {
  #opening clinical table
  clinical_geo <- read_xlsx(path = "GEO_clinical.xlsx")
  clinical_geo <- as.data.frame(clinical_geo)
  NameExpTable <- paste0(tissueID,"-clean-ann.RData")
  #aim folder of expression table
  expWD <- paste0("TNMdata/GEO_Expression-Tables/", NameExpTable)
  #opening expression table
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- clinical_geo[clinical_geo$Tumor == 1 & clinical_geo$Origin == tissueID, 1]
    
  }else if (tissueType == "Normal") {
    usedTissue <- clinical_geo[clinical_geo$Normal == 1 & clinical_geo$Origin == tissueID, 1]
  }else if (tissueType == "Metastatic") {
    usedTissue <- clinical_geo[clinical_geo$Meta == 1 & clinical_geo$Origin == tissueID, 1]
  }
  
  rownames(expTable) <- expTable$gene_sym
  expTable <- expTable[,-c(1,2)]
  expsNeeded <- as.matrix(expTable[,colnames(expTable) %in% usedTissue])
  #make expression vectors from patients with common values 
  geneAExp <- expsNeeded[rownames(expsNeeded) == geneA,]
  geneBExp <- expsNeeded[rownames(expsNeeded) == geneB,]
  
  #select cor test
  corrTest <- cor.test(x = geneAExp, y = geneBExp, method = ifelse(test = corrMethod == "Pearson", yes = "pearson", no = "spearman"), exact = F)
    #create summary data
    summRes <- data.frame("geneANum" = length(geneAExp),
                          "geneBNum" = length(geneBExp),
                          "R" = corrTest$estimate,
                          "pValue" = corrTest$p.value,
                          "Correlation type" = ifelse(test = corrMethod == "Pearson", yes = "Pearson correlation", no = "Spearman correlation"))
    colnames(summRes) <- c(paste0(geneA, " N"), paste0(geneB, " N"), "R", "p", "Analysis type")
  
  #create df for visualization
  dat <- data.frame("geneAExp" = log(x = geneAExp,base = 10), "geneBExp" = log(x = geneBExp,base = 10))
  library(ggplot2)
  scatPlot <- ggplot(dat, aes(x = geneAExp, y = geneBExp)) +
    geom_smooth(method = "lm", size = 2) +
    geom_point(size = 4) +
    xlab(paste0("Log10 expression of ", geneA)) + 
    ylab(paste0("Log10 expression of ", geneB)) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2),
          axis.title.x = element_text(size = 26, face = "bold",  vjust = 2)
          )
  res <- list(scatPlot = scatPlot, SumStat = summRes)
  
  
  
}

polyCorr_GEO <- function(genes, corrMethod, tissueID, tissueType){
  
  #opening clinical table
  clinical_geo <- read_xlsx(path = "GEO_clinical.xlsx")
  clinical_geo <- as.data.frame(clinical_geo)
  NameExpTable <- paste0(tissueID,"-clean-ann.RData")
  #aim folder of expression table
  expWD <- paste0("TNMdata/GEO_Expression-Tables/", NameExpTable)
  #opening expression table
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- clinical_geo[clinical_geo$Tumor == 1 & clinical_geo$Origin == tissueID, 1]
    
  }else if (tissueType == "Normal") {
    usedTissue <- clinical_geo[clinical_geo$Normal == 1 & clinical_geo$Origin == tissueID, 1]
  }
  
  rownames(expTable) <- expTable$gene_sym
  expTable <- expTable[,-c(1,2)]
  #create expression matrix for the genes needed 
  expsNeeded <- as.matrix(expTable[rownames(expTable) %in% genes, colnames(expTable) %in% usedTissue])
  
  #perform correlation analysis
  corrMatrix <- cor(t(expsNeeded), method = ifelse(test = corrMethod == "Pearson", yes = "pearson", no = "spearman")) 
  #create df for visualization
  corrs <- corrMatrix
  corrs[upper.tri(x = corrMatrix, diag = T)] <- NA
  corrs <- as.data.frame(as.table(corrs))
  corrs <- corrs[!is.na(corrs$Freq),]
  colnames(corrs) <- c("Gene1", "Gene2", "R")
  corrs$R <- round(x = corrs$R, digits = 2)
  
  library(ggplot2)
  corrHeatmap <- ggplot(data = corrs, aes(Gene1, Gene2, fill = R)) +
    geom_tile(color = "white") +
    geom_text(aes(Gene1, Gene2, label = R), color = "black", size = 5,fontface = "bold") +
    scale_fill_gradient2(low = "blue",mid = "white", high = "red", limit = c(-1,1), space = "Lab", 
                         name="Correlation\nCoefficient") +
    labs(x = NULL, y = NULL) +
    
    theme_minimal() + 
    theme(axis.text.x = element_text(size = 26, face = "bold",color = "black", angle = 90, hjust = 1.1),
          axis.text.y = element_text(size = 26,face = "bold",color = "black"),
          legend.text = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold")) +
    coord_fixed()
  
  res <- list(corrHeatmap = corrHeatmap)
  
  
}
