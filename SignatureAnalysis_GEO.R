signatureAnalysis_GEO <- function(tissueID, geneName){
  clinTable <- read_xlsx(path = "GEO_clinical.xlsx", sheet = 1)
  clinTable <- clinTable[clinTable$Origin == tissueID,]
  
  tumNames <- as.character(unlist(clinTable[clinTable$Tumor == 1, 1]))
  normNames <- as.character(unlist(clinTable[clinTable$Normal == 1, 1]))
  metNames <- as.character(unlist(clinTable[clinTable$Meta == 1, 1]))
  
  load(paste0("TNMdata/GEO_Expression-Tables/", tissueID,"-clean-ann.RData"))
  expTable <- txtFile
  rm(txtFile)
  
  expsNeeded <- expTable[expTable$gene_sym %in% geneName,]
  row.names(expsNeeded) <- expsNeeded$gene_sym
  expsNeeded <- expsNeeded[,-c(1,2)]
  
  signatures <- t(rbind(colMeans(expsNeeded)))
  
  if (length(metNames) >= 5) {
    plotTable <- rbind(data.frame("Means" = signatures[rownames(signatures) %in% normNames],
                                  "Type" = "Normal", row.names = names(signatures[rownames(signatures) %in% normNames,])),
                       data.frame("Means" = signatures[rownames(signatures) %in% tumNames],
                                  "Type" = "Tumor", row.names = names(signatures[rownames(signatures) %in% tumNames,])),
                       data.frame("Means" = signatures[rownames(signatures) %in% metNames],
                                  "Type" = "Metastatic", row.names = names(signatures[rownames(signatures) %in% metNames,])))
    
    kwTab <- kruskal.test(x = list("Normal" = plotTable[plotTable$Type == "Normal", 1],
                                   "Tumor" = plotTable[plotTable$Type == "Tumor", 1],
                                   "Metastatic" = plotTable[plotTable$Type == "Metastatic", 1]))
    tableToDisplay <- data.frame("Kruskall_Wallis_p" = format(x = kwTab$p.value, digits = 3, scientific = T),
                                 "FC_T-N" = median(plotTable[plotTable$Type == "Tumor", 1])/median(plotTable[plotTable$Type == "Normal", 1]),
                                 "FC_M-T" = median(plotTable[plotTable$Type == "Metastatic", 1])/median(plotTable[plotTable$Type == "Tumor", 1]),
                                 "N_genes" = nrow(expsNeeded))
    
  }else if (length(metNames) <= 5){
    plotTable <- rbind(data.frame("Means" = signatures[rownames(signatures) %in% normNames],
                                  "Type" = "Normal", row.names = names(signatures[rownames(signatures) %in% normNames,])),
                       data.frame("Means" = signatures[rownames(signatures) %in% tumNames],
                                  "Type" = "Tumor", row.names = names(signatures[rownames(signatures) %in% tumNames,])))
    
    mwTab <- wilcox.test(x = plotTable[plotTable$Type == "Normal", 1], y = plotTable[plotTable$Type == "Tumor", 1],paired = F, exact = F, conf.int = T)
    tableToDisplay <- data.frame("Mann_Whitney_p" = format(x = mwTab$p.value, digits = 3, scientific = T),
                                 "FC_T-N" = median(plotTable[plotTable$Type == "Tumor", 1])/median(plotTable[plotTable$Type == "Normal", 1]),
                                 "N_genes" = nrow(expsNeeded))
    
  }
  
  plotTable$Type <- factor(x = plotTable$Type, levels = if (length(metNames) >= 5){c("Normal", "Tumor", "Metastatic")} else c("Normal", "Tumor"), ordered = T)
  
  twoCol <- c("darkgreen", "red")
  threeCol <- c("darkgreen", "red", "darkorange3")
  
  twoFill <- c("darkolivegreen1", "coral1")
  threeFill <- c("darkolivegreen1", "coral1", "goldenrod1")
  
  #identification of outliers
  normOuts <- boxplot.stats(plotTable[plotTable$Type == "Normal",1])$stats
  tumOuts <- boxplot.stats(plotTable[plotTable$Type == "Tumor",1])$stats
  metOuts <- boxplot.stats(plotTable[plotTable$Type == "Metastatic",1])$stats
  outs <-  if (length(metNames) >= 5){range(normOuts, tumOuts, metOuts)} else range(normOuts, tumOuts)
  library(ggplot2)
  boxPlot <- ggplot(plotTable, aes(x = Type, y = Means, color = Type, fill = Type))  + 
    geom_boxplot(outlier.size = -1, notch = F, width = 0.55) +
    geom_jitter(width = 0.25) + 
    scale_color_manual(values= if (length(unique(plotTable$Type)) == 2) {twoCol} else if(length(unique(plotTable$Type)) == 3){threeCol}) +
    scale_fill_manual(values= if (length(unique(plotTable$Type)) == 2) {twoFill} else if(length(unique(plotTable$Type)) == 3){threeFill}) +
    ylim(outs[1], outs[2]) +
    labs(x = NULL, y = "Signature expression") + 
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  vioPlot <- ggplot(plotTable, aes(x = Type, y = Means, color = Type, fill = Type))  + 
    geom_violin(lwd = 1) + 
    geom_boxplot(lwd = 1.2,  width = 0.1, color = if (length(unique(plotTable$Type)) == 2) {twoCol} else if(length(unique(plotTable$Type)) == 3){threeCol}) +
    scale_fill_manual(values= if (length(unique(plotTable$Type)) == 2) {c("white", "white")} else if(length(unique(plotTable$Type)) == 3){c("white", "white", "white")}) +
    scale_color_manual(values= if (length(unique(plotTable$Type)) == 2) {twoCol} else if(length(unique(plotTable$Type)) == 3){threeCol}) +
    ylim(outs[1], outs[2]) +
    labs(x = NULL, y = "Signature expression") + 
    theme(legend.position = "none", 
          panel.background = element_rect(fill = NA, colour = "grey50"),  
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))
  
  res <- list(boxPlot = boxPlot, vioPlot = vioPlot, tableToDisplay  = tableToDisplay)
}


