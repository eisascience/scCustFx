
#' Plot Explained Variance of Principal Components
#'
#' This function generates a bar plot illustrating the explained variance of principal components (PCs) in a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param pcs Numeric vector specifying which PCs to plot. Default is 1:10.
#'
#' @return A ggplot object representing the explained variance of the specified principal components.
#'
#' @export
plot_pca_variance <- function(SerObj, pcs = 1:10) {
  
  # Get the explained variance by the PCs
  variance_explained <- SerObj@reductions$pca@stdev^2 / sum(SerObj@reductions$pca@stdev^2)
  
  # Subset data for the specified PCs
  variance_data <- data.frame(PC = pcs, Variance_Explained = variance_explained[pcs])
  
  # Plot the explained variance
  gg <- ggplot(variance_data, aes(x = factor(PC), y = Variance_Explained, fill = ifelse(PC > 3, "Pink", "Lightgray"))) + 
    geom_col() + 
    xlab("Principal Component") + 
    ylab("Explained Variance (%)") + 
    theme_classic(base_size = 14)+
    ggtitle(paste("Explained Variance of PCs", min(pcs), "to", max(pcs))) +
    scale_x_discrete(labels = paste("PC", pcs)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Pink", "Lightgray")) +
    geom_hline(yintercept = 0.05, linetype = "dashed") + 
    theme(legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, vjust = .8, hjust=.9))
  
  return(gg)
}


#' Create a grouped plot with percentage labels in the legend.
#'
#' This function generates a grouped plot for the specified feature in a Seurat object, 
#' adding percentage labels to the legend to represent the proportion of total points 
#' in each category.
#'
#' @param SerObj A Seurat object containing the data.
#' @param group_feature The feature based on which the data will be grouped and plotted.
#' @param color_values A vector of color values for the groups.
#' @param alpha The alpha (transparency) value for the points. Default is 0.3.
#' @param raster Logical, indicating whether raster graphics should be used. Default is FALSE.
#' @param raster.dpi A numeric vector of length 2 indicating the reSerObjlution of the raster 
#'   graphics. Default is c(2000, 2000).
#' @param pt.size The size of the points in the plot. Default is 0.5.
#' @param base_size The base font size for the plot. Default is 10.
#' @param theme_b The theme to be applied to the plot. Default is theme_bw with base_size parameter.
#' @return A ggplot object with grouped points and legend labels representing percentages.
#' @examples
#' \dontrun{
#' # Example usage:
#' # gg <- create_grouped_plot(SerObj, group_feature = "IsAlphaBeta", 
#' #                           color_values = c("dodgerblue", "maroon"), alpha = 0.3,
#' #                           raster = FALSE, raster.dpi = c(2000, 2000), pt.size = 0.5,
#' #                           base_size = 10, theme_b = theme_classic(base_size = 10))
#' }
#'
#' @export
create_grouped_plot <- function(SerObj, group_feature, color_values = c("dodgerblue", "maroon"), 
                                alpha = .3, raster=F, raster.dpi = c(2000, 2000), 
                                pt.size = .5,base_size = 10, shuffle = F,
                                theme_b = theme_classic(base_size = base_size)) {
  # Extract the specified group_feature column from the Seurat object
  group_column <- SerObj[[group_feature]]
  
  # Calculate percentage of total points for the specified group_feature
  percentages <- round(table(group_column) / nrow(group_column) * 100 , 2)
  
  gg <- DimPlot(SerObj, shuffle = shuffle, 
                group.by = group_feature, 
                raster = raster, raster.dpi = raster.dpi, 
                cols = alpha(color_values, alpha),
                pt.size = pt.size, order = T) +
    scale_color_manual(values = color_values, 
                       breaks = names(percentages), 
                       labels = paste(names(percentages), 
                                      " (", round(percentages, 1), "%)")) +
    theme_b +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    labs(color = group_feature)
  
  return(gg)
}


#' Plot 3D visualization of Seurat object
#'
#' @param SerObj A Seurat object
#' @param feature1 The name of the first feature to plot along the x-axis
#' @param feature2 The name of the second feature to plot along the y-axis
#' @param feature3 The name of the third feature to plot along the z-axis
#' @param color_feature The name of the feature to use for coloring the points (optional)
#' @param cols The color palette to use for the points (optional)
#' @param plot_title The title of the plot (optional)
#'
#' @return A 3D plotly visualization of the Seurat object
#' @export
plot3d_seurat <- function(SerObj, 
                          feature1, feature2, feature3, 
                          color_feature=NULL, 
                          cols = colorRamp(c("navy", "red")), 
                          plot_title="", slot = "data", markerSize=2) {
  library(plotly)
  
  
  
  # Check if feature1 is a gene name or metadata name
  # if (feature1 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   x <- SerObj@assays$RNA@data[feature1, ]
  # } else if (feature1 %in% colnames(SerObj@meta.data)) {
  #   x <- SerObj@meta.data[, feature1]
  # } else {
  #   stop(paste("Feature", feature1, "not found."))
  # }
  x = FetchData(SerObj, feature1, slot = slot)[,1]
  
  # Check if feature2 is a gene name or metadata name
  # if (feature2 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   y <- SerObj@assays$RNA@data[feature2, ]
  # } else if (feature2 %in% colnames(SerObj@meta.data)) {
  #   y <- SerObj@meta.data[, feature2]
  # } else {
  #   stop(paste("Feature", feature2, "not found."))
  # }
  y = FetchData(SerObj, feature2, slot = slot)[,1]
  
  # Check if feature3 is a gene name or metadata name
  # if (feature3 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   z <- SerObj@assays$RNA@data[feature3, ]
  # } else if (feature3 %in% colnames(SerObj@meta.data)) {
  #   z <- SerObj@meta.data[, feature3]
  # } else {
  #   stop(paste("Feature", feature3, "not found."))
  # }
  z = FetchData(SerObj, feature3, slot = slot)[,1]
  
  
  tempDF = data.frame(x=x, y=y, z=z)
  
  # Set point colors if color_feature is specified
  if (!is.null(color_feature)) {
    
    if (color_feature %in% colnames(SerObj@meta.data)) {
      color_vals <- SerObj@meta.data[, color_feature]
      
      if(!is.factor(color_vals)) tempDF$Cols = naturalsort::naturalfactor(color_vals) else tempDF$Cols = color_vals
      
      cols_red = cols[1:length(levels(tempDF$Cols))]
      names(cols_red) = levels(tempDF$Cols)
      
      p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, color = ~Cols, 
                   colors = cols_red, type = "scatter3d", mode = "markers", 
                   marker = list(size = markerSize))
    } else {
      GeneExpr = FetchData(SerObj, color_feature, slot = slot)[,1]
      p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, 
                   color = ~GeneExpr, colors = cols, 
                   type = "scatter3d", mode = "markers", marker = list(size = markerSize))
      
    }
  } else {
    # Set default color
    p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", marker = list(size = 2))
  }
  
  # Set plot title and axis labels
  # plot_title <- paste(feature1, "vs", feature2, "vs", feature3)
  x_axis_label <- feature1
  y_axis_label <- feature2
  z_axis_label <- feature3
  
  p <- layout(p, title = plot_title, scene = list(xaxis = list(title = x_axis_label), 
                                                  yaxis = list(title = y_axis_label), 
                                                  zaxis = list(title = z_axis_label)))
  
  return(p)
}


#' Calculate the percent expressed cells for each group in a Seurat object
#'
#' @param SerObj A Seurat object
#' @param group.by A meta feature used to group the cells
#' @param features A vector of features (genes) to calculate the percent expressed cells for
#' @param plot A logical indicating whether to plot 
#' @param topN select top N genes  
#' @param cols color vector
#' @return A data frame with the percent expressed cells for each group and feature
#' @export
VlnPlot_subset <- function(SerObj, group.by, features, plot=F, topN = 10, xlab = "",
                           cols=c("#8DD3C7", "#33A02C", "#F4CAE4", "#A6CEE3", "#FDCDAC",
                                  "#B2DF8A","#E78AC3", "#1F78B4", "#FB9A99", "#E31A1C", 
                                  "#66C2A5", "#FF7F00"), addDimplot = T, returnLs = T){
  
  
  #fix factor
  SerObj@meta.data[,group.by] = naturalsort::naturalfactor(as.character(SerObj@meta.data[,group.by]))
  
  features = unique(features)
  
  
  PctExprDF =  scCustFx:::PercentExpressed(SerObj, 
                                           group.by=group.by, 
                                           features=features, 
                                           plot_perHM = plot)
  
  PctExprDF$max = apply(PctExprDF, 1, max)
  PctExprDF$min = apply(PctExprDF, 1, min)
  
  PctExprDF$maxminDelta = PctExprDF$max-PctExprDF$min
  
  
  Idents(SerObj) = group.by
  
  vln <- VlnPlot(SerObj, 
                 features = rownames( PctExprDF[order(PctExprDF$maxminDelta, decreasing = T)[1:topN],]), 
                 stack = TRUE, flip = TRUE, fill.by = "ident",
                 cols = cols) + NoLegend() + xlab(xlab) #+
  # theme_classic(base_size = 14) 
  
  
  if(addDimplot){
    dim <- DimPlot(SerObj, 
                   label = T, label.size = 10, label.box = T, repel = T,
                   cols=cols,
                   group.by=group.by) + 
      coord_flip()  + scale_y_reverse() +
      theme_classic(base_size = 14) + 
      NoLegend() +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank()
      )
    if(returnLs) {
      list(dim=dim, vln=vln)
    } else {
      dim / vln
    }
  } else  {
    vln
  }
  
  
}

#' @title PlotFeatThrTab
#' @description This function from a seurat object, takes the feature selected and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SerObj, A Seurat object.
#' @param FeatName, the Feature name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotFeatThrTab <- function(SerObj, FeatName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(FeatName)) stop("FeatName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(FeatName) != 1) stop ("FeatName Length != 1")
  
  tempDF <- FetchData(SerObj, vars = c(FeatName, MetaDataName))
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
    
    
    
    vp <- VlnPlot(SerObj, features = FeatName, 
                  group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=var1)) + 
        geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", binwidth = .03) +
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
#' @param SerObj, A Seurat object.
#' @param GeneName, the Gene name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotExprThrTab <- function(SerObj, GeneName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(GeneName)) stop("GeneName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(GeneName) != 1) stop ("GeneName Length != 1")
  
  
  # tempTab <- table(SerObj@assays$RNA@data[GeneName, ]>CutThresh, SerObj@meta.data[,MetaDataName])
  tempTab <- table(GetAssayData(object = SerObj, slot = 'data')[GeneName, ]>CutThresh, 
                   SerObj@meta.data[,MetaDataName])
  
  
  
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
    
    tempDF  <- as.data.frame(cbind(Expr = SerObj@assays$RNA@data[GeneName, ], CellN = colnames(SerObj@assays$RNA@data)), stringsAsFactors = F)
    tempDF$Expr <- as.numeric(tempDF$Expr)
    
    vp <- VlnPlot(SerObj, features = GeneName, group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
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
#' @param SerObj A Seurat Object 
#' @param min.cells a horizontal line is plotted at this value
#' @param features Features
#' @param slot which ADT slot data or count
#' @param thresh a vertical line is plotted at this value and perc. of cells below it is computed and plotted
#' @return A ggplot object
#' @export
QuantileADTcounter <- function(SerObj, features=NULL, min.cells=10, 
                               slot="data", thresh = NULL){
  print("Setting default assay to ADT")
  DefaultAssay(SerObj) = "ADT"
  
  ADT_range <- (GetAssayData(SerObj, slot=slot))
  # max(ADT_range)
  # ADT_range[1:10, 1:10]
  
  ADT_range.melt <- reshape2::melt(quantile(ADT_range, 
                                            SerObjrt(c(1:9.9/100, 
                                                   1:9/10, .75, .85, .995, .996, .997, 
                                                   .998, .999, .9993, .9995, .9998, .9999,
                                                   round(seq(from=90, to=100, length.out = 50), 2)/100)) ))
  
  ADT_range.melt$value <- round(ADT_range.melt$value, 3)
  ADT_range.melt <- unique(ADT_range.melt)
  
  # ADT_range.melt$Quant  <- as.numeric(gsub("%", "", rownames(ADT_range.melt)))/100
  
  if(is.null(features)) {
    warning("features was NULL, using all features of ADT")
    features = rownames(SerObj)
  }
  
  ggls <- lapply(features, function(adt_f, ADT.DF){
    
    # adt_f = "CD8a"
    
    ADT_range.melt.sub <- ADT_range.melt
    
    ADT_range.melt.sub$Bin  <- as.numeric(table(cut(GetAssayData(SerObj, slot = slot)[adt_f, ], 
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
#' @param SerObj A Seurat Object 
#' @param assay Assay within the seurat object
#' @param slot which ADT slot data or count
#' @param featureSplit a feature of metadata; if NULL no split id done. 
#' @return A ggplot object
#' @export
plot_libSize <- function(SerObj, assay="ADT", slot = "counts", 
                         featureSplit=NULL, col_vector = NULL){
  
  if(is.null) col_vector = scCustFx:::ColorTheme()$col_vector
  
  if(is.null(featureSplit)) stop("featureSplit is null")
    
  if(!is.null(featureSplit)) {
    
    if(!featureSplit %in% colnames(SerObj@meta.data)) {
      stop("featureSplit not in meta.data")
      
    } else {
      
      SerObj$tempLibSize = colSums(GetAssayData(SerObj, assay = assay, slot=slot))
      
      SerObj$tempGrouping = factor(SerObj@meta.data[,featureSplit])
      
      tempOut = SerObj@meta.data %>% group_by(tempGrouping)  %>% 
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
                               
                               
