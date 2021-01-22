

#' @title QuantileADTcounter
#'
#' @description Plots a histogram of ADT/protein/CITE-seq markers (optional: identify percent of cells cut by a threshold)
#' @param so A Seurat Object 
#' @param min.cells a horizontal line is plotted at this value
#' @param slot which ADT slot data or count
#' @param thresh a vertical line is plotted at this value and % of cells below it is computed and plotted
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
                               
                               
