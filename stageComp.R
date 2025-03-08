nextComparer <- function(tissueID, selGene, clinTab){
  selSubtype <- "Stage"
  clinTab <- as.data.frame(read_xlsx("GEO_clinical_mergedStages.xlsx", sheet = 1))
  clinTab <- clinTab[clinTab$Origin == tissueID,]
  NameExpTable <- paste0(tissueID,"-clean-ann.RData")
  #aim folder of expression table
  expWD <- paste0("TNMdata/GEO_Expression-Tables/", NameExpTable)
  #opening expression table
  load(expWD)
  expTable <- txtFile
  rm(txtFile)
  rownames(expTable) <- expTable$gene_sym
  expTable <- expTable[,-c(1,2)]

  normNames <- clinTab[clinTab$Normal == 1, colnames(clinTab) == "GSM"]
  norms <- clinTab[clinTab$Normal == 1, colnames(clinTab) %in% c("GSM", selSubtype)]
  norms$Stage <- "Normal"
  tums <- clinTab[clinTab$Tumor == 1, colnames(clinTab) %in% c("GSM", selSubtype)]
  tums <- tums[!is.na(tums$Stage),]
  tums <- rbind(tums, norms)
  redExp <- as.data.frame(t(expTable[rownames(expTable) == selGene,colnames(expTable) %in% c(normNames,tums$GSM)]))
  redExp$GSM <- rownames(redExp)
  redExp <- merge(x = redExp,y = tums,by = "GSM")
  colnames(redExp) <- c("GSM", "Exp", "Stage")

  kw <- kruskal.test(Exp~Stage,data = redExp)
  outs <- range(boxplot.stats(redExp$Exp)$stats)
  pvalue <-format(kw$p.value, digits = 3, scientific = T)
  dunnTest <- dunn.test(x = redExp$Exp,g = redExp$Stage,label = T, kw = T)
  statSum <- data.frame("Comparisons" = dunnTest$comparisons,
                        "Dunn_Test_pAdj" = dunnTest$P.adjusted)
  sumTab <- data.frame()
  for (var in 1:length(unique(redExp$Stage))) {
   selVar <-  unique(redExp$Stage)[var]
   tab <- c(selVar,fivenum(redExp[redExp$Stage == selVar,2]),boxplot.stats(redExp[redExp$Stage == selVar,2])$stats[1],boxplot.stats(redExp[redExp$Stage == selVar,2])$stats[5], length(redExp[redExp$Stage == selVar,2]))
   tab <- t(as.data.frame(x = tab,row.names = c("Stage","Min", "Q1", "Med", "Q3", "Max","Lower whisker", "Upper whisker", "N")))
   rownames(tab) <- selVar
   sumTab <- rbind(sumTab, tab)
  }
  
  sumTab <- sumTab[order(rownames(sumTab)),]
  stages <- unique(redExp[redExp$Stage != "Normal",3])
  levels <- c("Normal", sort(stages))
  stage4Color <- c("darkgreen", "goldenrod1", "blue", "purple", "red")
  stage3Color <- c("darkgreen", "goldenrod1", "purple", "red")
  stage4fill <- c("darkolivegreen1", "goldenrod1", "blue", "purple", "red")
  stage3fill <- c("darkolivegreen1", "goldenrod1", "purple", "red")
  
  
  redExp$Stage <- factor(x = redExp$Stage, levels = levels)
  boxPlot <- ggplot(redExp, aes(x = Stage, y = Exp, color = Stage, fill = Stage))  +
    geom_boxplot(outlier.size = -1, notch = F, width = 0.55,alpha = 0.3)+ ylim(outs[1], outs[2]) +
    geom_jitter(width = 0.25, aes(fill = Stage)) +
    labs(x = NULL, y = paste0(selGene, " gene expression")) +
    geom_label(label = paste0("KWp = ", pvalue), fill = "white", color = "black",x =(length(unique(redExp$Stage))),
               y = outs[2] , label.size = 0, size = 7) +
    scale_color_manual(values = if(length(stages) == 4){stage4Color}else if(length(stages) == 3){stage3Color})+
    scale_fill_manual(values = if(length(stages) == 4){stage4fill}else if(length(stages) == 3){stage3fill})+
    theme(legend.position = "none",
          panel.background = element_rect(fill = NA, colour = "grey50"),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 26, colour = "black"),
          axis.text.y = element_text(size = 26, colour = "black"),
          axis.title.y = element_text(size = 26, face = "bold",  vjust = 2))


  res <- list(boxPlot = boxPlot, kwp = pvalue,sumTab = sumTab, expTab = redExp[,-1], statSum = statSum)
}