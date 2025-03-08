generalPlotter <- function(geneName,tissueType, genPlotSourceTable){
  tissueId <- read_xlsx(path = "genPlot_tissueID.xlsx", col_names =  F)
  tissueId <- as.data.frame(tissueId[tissueId$...3 %in% tissueType,])
  tissueId <- as.list(tissueId)
  plotNums <- read.table(file = "panCcNUmbers.txt",header = T,sep = "\t",row.names = 1)

  bigTable <- genPlotSourceTable[genPlotSourceTable$Gene == geneName,]
  plotTable <- data.frame(matrix(ncol= 7,nrow=0, dimnames=list(NULL, c("ExPerc", "percent", "type", "gene", "tissue", "pvalue", "signif"))))

  for (tissue in 1:length(bigTable$Tissue)){
    normTable <- data.frame("ExPerc" = unlist(bigTable[tissue,4:8]),
                            "percent" = c(0, 25, 50, 75, 100),
                            "type" = "Normal (left)",
                            "gene" = geneName,
                            "tissue" = bigTable$Tissue[tissue],
                            "pvalue" = bigTable$M.Wp[tissue],
                            "signif" = bigTable$signif[tissue])
    tumTable <- data.frame("ExPerc" = unlist(bigTable[tissue,9:13]),
                           "percent" = c(0, 25, 50, 75, 100),
                           "type" = "Tumor (right)",
                           "gene" = geneName,
                           "tissue" = bigTable$Tissue[tissue],
                           "pvalue" = bigTable$M.Wp[tissue],
                           "signif" = bigTable$signif[tissue])
    tissueTable <- rbind(normTable, tumTable)

    plotTable <- rbind(plotTable, tissueTable)

  }
  
  plotNums <- plotNums[,colnames(plotNums) %in% tissueType]
  plotTable <- plotTable[plotTable$tissue %in% tissueType,]
  signifTissue <- unique(plotTable[plotTable$signif == 1,5])
  library(ggplot2)
  expPlot <- ggplot(plotTable, aes(x = tissue, y = ExPerc, fill = tissue,line = type, size = type))  +
    geom_vline(xintercept=seq(1.5, 22-0.5, 1), colour="grey")+
    geom_boxplot(outlier.size = 1, notch = F, width = 0.75) +
    scale_x_discrete(labels = ifelse(tissueId$...3 %in% signifTissue, paste0(tissueId$...3, "*"), tissueId$...3)) +
    scale_size_manual(values = c(0.5,1.3))+
    guides(fill = "none") +
    labs(x = NULL, y = paste0(geneName, " gene expression")) +
    labs(fill = "Tumor type:", size = " ", alpha = "Significance under 0.01:") +
    theme(legend.position = "bottom",
          panel.background = element_rect(fill = NA, colour = "grey50"),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey"),
          panel.grid.minor.y = element_line(colour = "grey"),
          axis.text.x = element_text(size = 13, colour = ifelse(tissueId$...3 %in% signifTissue, "red", "black"), face = ifelse(tissueId$...3 %in% signifTissue, "bold", "plain")),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title.y = element_text(size = 13,  vjust = 2, face = "bold"),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13))
  
  res <- list(generalPlot = expPlot, plotNums = plotNums)
  
}

dotPlotter <- function(geneName, tissueType, genPlotSourceTable, tissueId){
  bigTable <- genPlotSourceTable[genPlotSourceTable$Gene %in% geneName,]
  bigTable$FC <-  bigTable$Med.Tum/bigTable$Med.Norm
  bigTable$pAdj <- p.adjust(bigTable$M.Wp,method = "fdr")
  
  plotTable <- bigTable[,colnames(bigTable) %in% c("Gene","Tissue","pAdj","FC" )]
  plotTable <- plotTable[plotTable$Tissue %in% tissueType,]
  plotTable$FC <- ifelse(test = plotTable$FC %in% c("Inf", "NaN"), yes = 0,no = plotTable$FC) 
  plotTable$pAdj <- ifelse(plotTable$FC == 0, yes = 1,no = plotTable$pAdj)
  plotTable$signif <- ifelse(plotTable$pAdj < 0.01, yes = 1,no = 0)
  
  dotPlot <- ggplot(data = plotTable,aes(x = Tissue,y = Gene, color = log2(FC), size = -log10(pAdj))) +
    geom_point() +
    geom_vline(xintercept=seq(from = 1.5, to = length(unique(plotTable$Tissue)) - 0.5, by = 1), colour="grey")+
    geom_hline(yintercept=seq(from = 1.5, to = length(unique(plotTable$Gene)) - 0.5, by = 1), colour="grey")+
    
    scale_color_gradient2(low = "#28F903",mid = "white", high = "#F93333", space = "Lab") +
    geom_point(data = plotTable[plotTable$signif == 1,],colour = "black", pch = 21) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      panel.background = element_rect(fill = NA, colour = "white"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 13, colour = "black",face = "bold", angle = 67.5,hjust=1),
      axis.text.y = element_text(size = 13, colour = "black",face = "bold"),
      axis.title = element_blank(),
      legend.text = element_text(size = 10,face = "bold"),
      legend.title = element_text(size = 10,face = "bold")
      )
  
  res <- list(dotPlotPanCan = dotPlot)
  
} 