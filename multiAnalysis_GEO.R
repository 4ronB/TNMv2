multipleGeneAnalysis_GEO <- function(tissueID, geneName){
  clinTable <- read_xlsx(path = "GEO_clinical.xlsx", sheet = 1)
  clinTable <- clinTable[clinTable$Origin == tissueID,]
  
  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  metNames <- as.character(unlist(clinTable[clinTable$Meta == 1, 1]))
  
  
  load(paste0("TNMdata/GEO_Expression-Tables/", tissueID,"-clean-ann.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  expsNeeded <- expTable[expTable$gene_sym %in% geneName,]
  expsNeeded <- expsNeeded[,-1]
  expsNeeded <- as.data.frame(t(expsNeeded))
  colnames(expsNeeded) <- expsNeeded[1,]
  expsNeeded <- expsNeeded[-1,]
  expsNeeded$ID <- rownames(expsNeeded)
  expsNeeded$Type <- ifelse(test = expsNeeded$ID %in% normNames, yes = "Normal", 
                            no = ifelse(test = expsNeeded$ID %in% tumNames, yes = "Tumor",
                                        no = ifelse(test = expsNeeded$ID %in% metNames, yes = "Metastatic",no = "pleaseRemoveMe")))
  
  expsNeeded <- expsNeeded[expsNeeded$Type != "pleaseRemoveMe",]
  expsNeeded <- expsNeeded[,c(ncol(expsNeeded), ncol(expsNeeded)-1, 1:(ncol(expsNeeded)-2))]
  
  #create table for ggplot
  tableForPlot <- data.frame()
  for (gene in 3:ncol(expsNeeded)) {
    selGene <- colnames(expsNeeded)[gene]
    table <- data.frame("Type" = expsNeeded$Type, 
                        "Expression" = expsNeeded[,colnames(expsNeeded) == selGene], 
                        "Gene" = selGene)
    
    tableForPlot <- rbind(tableForPlot, table)
  }
  tableForPlot$Expression <- as.numeric(tableForPlot$Expression)
  if (length(metNames) >= 1) {
    levels <-  c("Normal", "Tumor", "Metastatic")
  }else {levels <-  c("Normal","Tumor")}
  
  tableForPlot$Type <- factor(x = tableForPlot$Type, levels = levels)
  outData <- data.frame()
  for (gene in 1:length(unique(tableForPlot$Gene))) {
    selGene <- unique(tableForPlot$Gene)[gene]
    selTab <- tableForPlot[tableForPlot$Gene == selGene,]
    if (length(metNames) >= 1) {
      selOutData <- data.frame("Gene" = selGene,
                               "FC_TvsN" = mean(selTab[selTab$Type == "Tumor",2])/ mean(selTab[selTab$Type == "Normal",2]),
                               "FC_MvsT" = mean(selTab[selTab$Type == "Metastatic",2])/ mean(selTab[selTab$Type == "Tumor",2]),
                               "FC_MvsN" = mean(selTab[selTab$Type == "Metastatic",2])/ mean(selTab[selTab$Type == "Normal",2]))
      
    } else {selOutData <- data.frame("Gene" = selGene,
                                     "FC_TvsN" = mean(selTab[selTab$Type == "Tumor",2])/ mean(selTab[selTab$Type == "Normal",2]))}
    
    outData <- rbind(outData, selOutData)
    
  }
  
  
  outData <- outData[order(outData$Gene,decreasing = T),]
  tableForPlot$Expression <- ifelse(test = tableForPlot$Expression == 0, yes = 1,no = tableForPlot$Expression)
  
  twoCol <- c("green", "red")
  threeCol <- c("green", "red",  "black")
  library(ggplot2)
  library(ggridges)
  meanPlot <- ggplot(data=tableForPlot, aes(x=log2(Expression), y=Gene, col= Type)) +
    geom_boxplot() +
    theme_minimal() +
    xlab("log2(gene expressions)") +
    labs(fill = "") + 
    scale_color_manual(values= if (length(unique(tableForPlot$Type)) == 2) {twoCol} else if(length(unique(tableForPlot$Type)) == 3){threeCol}) +
    theme(legend.position="top",
          axis.text.x = element_text(size = 13, colour = "black"),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title.y = element_text(size = 13, face = "bold"),
          axis.title.x = element_text(size = 13, face = "bold")
          )
  
  densityPlot <- ggplot(data=tableForPlot, aes(x= log2(Expression), y=Gene, fill= Type)) +
    geom_density_ridges(alpha = 0.5) + 
    theme_minimal() +
    xlab("log2(gene expressions)") +
    labs(fill = "") + 
    theme(legend.position="top",
          axis.text.x = element_text(size = 13, colour = "black"),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title.y = element_text(size = 13, face = "bold"),
          axis.title.x = element_text(size = 13, face = "bold"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13))
  
  res <- list(meanPlot = meanPlot, densityPlot = densityPlot, outData = outData)
}

targetgram_GEO <- function(tissueID, geneName,logtrans){
  clinTable <- read_xlsx(path = "GEO_clinical.xlsx", sheet = 1)
  clinTable <- clinTable[clinTable$Origin == tissueID,]
  clinTable <- clinTable[!clinTable$Pair %in% "ismetlodo case de nincs par", ]
  
  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  metNames <- as.character(unlist(clinTable[clinTable$Meta == 1, 1]))
  
  load(paste0("TNMdata/GEO_Expression-Tables/", tissueID,"-clean-ann.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  expsNeeded <- expTable[expTable$gene_sym %in% geneName,]
  expsNeeded <- expsNeeded[,-1]
  expsNeeded <- as.data.frame(t(expsNeeded))
  colnames(expsNeeded) <- expsNeeded[1,]
  expsNeeded <- expsNeeded[-1,]
  expsNeeded$ID <- rownames(expsNeeded)
  expsNeeded$Type <- ifelse(test = expsNeeded$ID %in% normNames, yes = "Normal", 
                            no = ifelse(test = expsNeeded$ID %in% tumNames, yes = "Tumor",
                                        no = ifelse(test = expsNeeded$ID %in% metNames, yes = "Metastatic",no = "pleaseRemoveMe")))
  
  expsNeeded <- expsNeeded[expsNeeded$Type != "pleaseRemoveMe",]
  expsNeeded <- expsNeeded[,c(ncol(expsNeeded), ncol(expsNeeded)-1, 1:(ncol(expsNeeded)-2))]
  
  #create table for ggplot
  tableForPlot <- data.frame()
  for (type in 1:length(unique(expsNeeded$Type))) {
    selType <- unique(expsNeeded$Type)[type]
    tab <- expsNeeded[expsNeeded$Type == selType,]
    
    for (gene in 3:ncol(expsNeeded)) {
      selGene <- colnames(expsNeeded)[gene]
      geneSum <- data.frame("Gene" = selGene,
                            "Mean" = mean(as.numeric(tab[,colnames(tab) == selGene])),
                            "Median" = median(as.numeric(tab[,colnames(tab) == selGene])),
                            "Type" = selType)
      tableForPlot <- rbind(tableForPlot,geneSum)
    }
  }
  
  if (logtrans == TRUE) {
    tableForPlot$Mean <- log2(tableForPlot$Mean)
    tableForPlot$Median <- log2(tableForPlot$Median)
    
  }
  
  
  tableForPlot$Type <- factor(x = tableForPlot$Type,levels = c("Normal", "Tumor", "Metastatic"))
  plt <- ggplot(tableForPlot, aes(color = Type, fill = Type)) +
    geom_col(aes(x = Gene, y = Mean),position = "dodge2",show.legend = TRUE,alpha = .9) +
    geom_point(aes(x = Gene, y = Median),size = 3) +
    geom_segment(aes(x = Gene, xend = Gene,y = 0,yend = Median),linetype = "dashed") +
    coord_polar() +
    labs(x = NULL, y = NULL) +
    facet_wrap(~Type, scales = "fixed") +
    scale_color_manual(values = c("darkgreen", "red", "darkorange3")) +
    scale_fill_manual(values = c("darkolivegreen1", "coral1", "goldenrod1")) +
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = "grey50"),
          strip.text.x = element_text(size = 18),
          panel.grid.major.x = element_line(colour = "grey"),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.minor.y = element_line(colour = "grey"),
          axis.text.x = element_text(size = 13, colour = "black"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  res <- list(targetPlot =plt)

}


