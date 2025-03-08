multiCorr_GDC <- function(geneName, tissueID, cutOff, tissueType){
  library(DT)
  #opening clinical table
  clinTable <- read_xlsx(path = "GDC_clinical.xlsx")
  clinTable <- as.data.frame(clinTable)
  TCGAID <- unlist(strsplit(tissueID[[1]], split = ".", fixed = T))[2]
  
  clinTable <- clinTable[clinTable$`Project ID` == TCGAID,]
  
  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  #opening expression table
  load(paste0("TNMdata/GTEx_Expression-Tables/", tissueID[[1]],"-raw-clean-norm.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  gtexData <-  data.frame("SampleID" = colnames(expTable)[grep(pattern = "GTEX", x = colnames(expTable))], "type" = "Normal")
  normNames <- c(normNames, gtexData$SampleID)
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
    
  }else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  expsNeeded <- expTable[,colnames(expTable) %in% usedTissue]
  
  #make the selected gene to a numeric vector and perform correlation test
  selectedGene <- as.vector(unlist(expTable[expTable$geneid == geneName,colnames(expTable) %in% usedTissue]))
  corResult <- apply(expsNeeded, 1, function(x) cor.test(x, selectedGene, method = "spearman", exact = F)[c(4,3)])
  
  #store the results as a dataframe and convert the list to a dataframe
  result = as.data.frame(matrix(nrow = nrow(expsNeeded), ncol = 4))
  colnames(result) <- c("ID", "Gene Name" , "Spearman's Rho", "p-value")
  
  result[, 1] = names(corResult)
  result[, 2] = expTable$geneid
  result[, 3] = round(x = sapply(corResult, "[[", 1), digits = 2)
  result[, 4] = format(x = sapply(corResult, "[[", 2), digits = 3, scientific = T)
  
  result <- result[,-1]
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
  
  res <- list(result = result, networkPlot = networkPlot)
  return(res)
}


singleCorr_GDC <- function(geneA, geneB, corrMethod, tissueID, tissueType) {
  #opening clinical table
  clinTable <- read_xlsx(path = "GDC_clinical.xlsx")
  clinTable <- as.data.frame(clinTable)
  TCGAID <- unlist(strsplit(tissueID[[1]], split = ".", fixed = T))[2]
  
  clinTable <- clinTable[clinTable$`Project ID` == TCGAID,]

  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  #opening expression table
  load(paste0("TNMdata/GTEx_Expression-Tables/", tissueID[[1]],"-raw-clean-norm.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  gtexData <-  data.frame("SampleID" = colnames(expTable)[grep(pattern = "GTEX", x = colnames(expTable))], "type" = "Normal")
  normNames <- c(normNames, gtexData$SampleID)
  
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
    
  }else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  geneAExp <- as.numeric(unlist(expTable[expTable$geneid == geneA,colnames(expTable) %in% usedTissue]))
  geneBExp <- as.numeric(unlist(expTable[expTable$geneid == geneB,colnames(expTable) %in% usedTissue]))
  
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
  dat <- data.frame("geneAExp" = log(x = geneAExp, base = 10), "geneBExp" = log(x = geneBExp, base = 10))
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
          axis.title.x = element_text(size = 26, face = "bold",  vjust = 2))
  res <- list(scatPlot = scatPlot, SumStat = summRes)
 
  
}


polyCorr_GDC <- function(genes, corrMethod, tissueID, tissueType){
  
  #opening clinical table
  clinTable <- read_xlsx(path = "GDC_clinical.xlsx")
  clinTable <- as.data.frame(clinTable)
  TCGAID <- unlist(strsplit(tissueID[[1]], split = ".", fixed = T))[2]
  
  clinTable <- clinTable[clinTable$`Project ID` == TCGAID,]
  
  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  #opening expression table
  load(paste0("TNMdata/GTEx_Expression-Tables/", tissueID[[1]],"-raw-clean-norm.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  gtexData <-  data.frame("SampleID" = colnames(expTable)[grep(pattern = "GTEX", x = colnames(expTable))], "type" = "Normal")
  normNames <- c(normNames, gtexData$SampleID)
  #screening for selected tissue in clinical table
  if (tissueType == "Tumor") {
    usedTissue <- tumNames
    
  }else if (tissueType == "Normal") {
    usedTissue <- normNames
  }
  
  #create expression matrix for the genes needed
  expsNeeded <- expTable[expTable$geneid %in% genes,c(1,which(colnames(expTable) %in% usedTissue))]
  
  #remove duplicated gene names
  if (TRUE %in% duplicated(expsNeeded$geneid)) {
    dupNames <- as.character(unlist(expsNeeded[which(duplicated(expsNeeded$geneid)),1]))
    idsToRemove <- c()
    for (name in 1:length(dupNames)) {
      dupped <- dupNames[name]
      dupTab <- expsNeeded[expsNeeded$geneid == dupped,]
      idToRemove <- ifelse(test = sum(dupTab[1,2:ncol(dupTab)]) - sum(dupTab[2,2:ncol(dupTab)]) > 0, yes = rownames(dupTab)[2], no = rownames(dupTab)[1])
      idsToRemove <- c(idsToRemove, idToRemove)
      
    }
    expsNeeded <- expsNeeded[!(rownames(expsNeeded) %in% idsToRemove),]
    
  }
  
  rownames(expsNeeded) <- expsNeeded$geneid
  expsNeeded <- expsNeeded[,-1]

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
  corrHeatmap <- ggplot(data = corrs, aes(x = Gene1, y = Gene2, fill = R, width = NULL)) +
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




