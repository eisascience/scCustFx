#' @title PlotFeatThrTab
#' @description This function from a seurat object, takes the feature selected and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SeuratObj, A Seurat object.
#' @param FeatName, the Feature name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotFeatThrTab <- function(SeuratObj, FeatName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(FeatName)) stop("FeatName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(FeatName) != 1) stop ("FeatName Length != 1")
  
  tempDF <- FetchData(SeuratObj, vars = c(FeatName, MetaDataName))
  colnames(tempDF) <- c("var1", "var2")
  tempTab <- table(tempDF$var1>CutThresh, tempDF$var2)
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    
    
    
    vp <- VlnPlot(SeuratObj, features = FeatName, 
                  group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=var1)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .03) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
  }
  
  
  
}

#' @title PlotExprThrTab
#' @description This function from a seurat object, takes the gene expression of a single gene and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SeuratObj, A Seurat object.
#' @param GeneName, the Gene name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotExprThrTab <- function(SeuratObj, GeneName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(GeneName)) stop("GeneName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(GeneName) != 1) stop ("GeneName Length != 1")
  
  
  # tempTab <- table(SeuratObj@assays$RNA@data[GeneName, ]>CutThresh, SeuratObj@meta.data[,MetaDataName])
  tempTab <- table(GetAssayData(object = SeuratObj, slot = 'data')[GeneName, ]>CutThresh, 
                   SeuratObj@meta.data[,MetaDataName])
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    
    tempDF  <- as.data.frame(cbind(Expr = SeuratObj@assays$RNA@data[GeneName, ], CellN = colnames(SeuratObj@assays$RNA@data)), stringsAsFactors = F)
    tempDF$Expr <- as.numeric(tempDF$Expr)
    
    vp <- VlnPlot(SeuratObj, features = GeneName, group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=Expr)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .3) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
    
    
    
  }
  
  
  
}

#' @title QuantileADTcounter
#'
#' @description Plots a histogram of ADT/protein/CITE-seq markers optional: identify percent of cells cut by a threshold
#' @param so A Seurat Object 
#' @param min.cells a horizontal line is plotted at this value
#' @param features Features
#' @param slot which ADT slot data or count
#' @param thresh a vertical line is plotted at this value and perc. of cells below it is computed and plotted
#' @return A ggplot object
#' @export
QuantileADTcounter <- function(so, features=NULL, min.cells=10, 
                               slot="data", thresh = NULL){
  print("Setting default assay to ADT")
  DefaultAssay(so) = "ADT"
  
  ADT_range <- (GetAssayData(so, slot=slot))
  # max(ADT_range)
  # ADT_range[1:10, 1:10]
  
  ADT_range.melt <- reshape2::melt(quantile(ADT_range, 
                                            sort(c(1:9.9/100, 
                                                   1:9/10, .75, .85, .995, .996, .997, 
                                                   .998, .999, .9993, .9995, .9998, .9999,
                                                   round(seq(from=90, to=100, length.out = 50), 2)/100)) ))
  
  ADT_range.melt$value <- round(ADT_range.melt$value, 3)
  ADT_range.melt <- unique(ADT_range.melt)
  
  # ADT_range.melt$Quant  <- as.numeric(gsub("%", "", rownames(ADT_range.melt)))/100
  
  if(is.null(features)) {
    warning("features was NULL, using all features of ADT")
    features = rownames(so)
  }
  
  ggls <- lapply(features, function(adt_f, ADT.DF){
    
    # adt_f = "CD8a"
    
    ADT_range.melt.sub <- ADT_range.melt
    
    ADT_range.melt.sub$Bin  <- as.numeric(table(cut(GetAssayData(so, slot = slot)[adt_f, ], 
                                                    breaks = c(-5, ADT_range.melt.sub$value), 
                                                    include.lowest = T)))
    
    
    ADT_range.melt.sub$Percent = ADT_range.melt.sub$Bin/sum(ADT_range.melt.sub$Bin)
    
    ADT_range.melt.sub$BinFac = factor(paste0("bin-", ADT_range.melt.sub$value), labels  = ADT_range.melt.sub$value)
    
    if(is.null(thresh)) thresh = as.numeric(quantile(ADT_range.melt.sub$value, .7))
    
    
    
    # sum(ADT_range.melt.sub$Percent[ADT_range.melt.sub$value < thresh])
    
    ThreshPerc  = sum(ADT_range.melt.sub[ADT_range.melt.sub$value < thresh,]$Bin)/sum(ADT_range.melt.sub$Bin)
    
    BinThresh <- which(ADT_range.melt.sub$BinFac ==
                         min(as.numeric(as.character(ADT_range.melt.sub$BinFac[as.numeric(as.character(ADT_range.melt.sub$BinFac))>thresh]))))
    
    ggplot(ADT_range.melt.sub) + 
      geom_bar(aes(x=BinFac, y=Bin, fill=Percent), position = 'dodge',stat='identity') + theme_bw()  +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      labs(title = paste0(adt_f, '\n No. of cells binned per global ADT quantiles\nPercent below Threshold (',round(thresh,2),') = %' , round(ThreshPerc*100, 3))) + 
      geom_hline(yintercept = min.cells) + geom_vline(xintercept = BinThresh, col="red", lty=2)
    
  })
  return(ggls)
}



#' @title plot_libSize
#' 
#' @description Plots a histogram of ADT/protein/CITE-seq markers (optional: identify percent of cells cut by a threshold)
#' @param so A Seurat Object 
#' @param assay Assay within the seurat object
#' @param slot which ADT slot data or count
#' @param featureSplit a feature of metadata; if NULL no split id done. 
#' @return A ggplot object
#' @export
plot_libSize <- function(so, assay="ADT", slot = "counts", 
                         featureSplit=NULL, col_vector = NULL){
  
  if(is.null) col_vector = scCustFx:::ColorTheme()$col_vector
  
  if(is.null(featureSplit)) stop("featureSplit is null")
    
  if(!is.null(featureSplit)) {
    
    if(!featureSplit %in% colnames(so@meta.data)) {
      stop("featureSplit not in meta.data")
      
    } else {
      
      so$tempLibSize = colSums(GetAssayData(so, assay = assay, slot=slot))
      
      so$tempGrouping = factor(so@meta.data[,featureSplit])
      
      tempOut = so@meta.data %>% group_by(tempGrouping)  %>% 
        summarize(mean_LibSize = sum(tempLibSize, na.rm = TRUE))
      
      
    }} 
  
  
  tempOut %>% 
    reshape2::melt() %>% 
    ggplot(aes(y=log10(value), x=tempGrouping, fill=tempGrouping)) + 
    geom_bar( position = 'dodge',stat='identity') +
    theme_classic() +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90)) + 
    ggtitle(paste0(Title, "\nSum log10 library size per barcode Prefix")) +
    scale_fill_manual(values=col_vector)
  
}
                               
                               
