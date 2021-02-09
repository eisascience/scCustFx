

#' @title heatmap_table
#'
#' @description Takes a table and plots it. Either as a table or a pheatmap heatmap with clustering
#' @param TBL A table (base R table()) 
#' @param clust.meth default is "ward.D2", see pheatmap for more options
#' @param clust.dist default is "euclidean", see pheatmap for more options
#' @param Title A title else ''
#' @param plotTable default is F, if T, it will plot an actual table
#' @return A pheatmap heatmap or a plotted table
#' @export
heatmap_table <- function(TBL, clust.meth="ward.D2", clust.dist="euclidean", Title="", plotTable=F){
  if(plotTable){
    grid.table(TBL)
  } else {
    pheatmap::pheatmap(as.matrix(TBL),
                       cluster_rows = T, cluster_cols = T,
                       main = paste0(Title, "\n", clust.dist, " dist w. ", clust.meth, "clustering"),
                       color = viridis(20), cellwidth = 10, show_colnames = T,
                       clustering_distance_rows = clust.dist,
                       clustering_distance_cols = clust.dist, #euclidean
                       clustering_method = clust.meth,
                       # annotation_col = ColAnnRow,
                       # annotation_colors = my_colour
                       #, labels_row = IDgenes
    )
  }
}


#' @title gg_barplot_table
#'
#' @description Takes a table and plots a gg barplot 
#' @param TBL A table (base R table()) 
#' @param Title A title else ''
#' @param gg_theme a ggplot theme or NULL
#' @return A ggplot 
#' @export
gg_barplot_table <- function(TBL, Title = "", gg_theme=NULL){
  
  tbl.m = reshape2::melt(TBL)
  if("Var1" %in% colnames(tbl.m)) {
    if(class(tbl.m$Var1)=="integer") tbl.m$Var1 = as.character(tbl.m$Var1)
    ggp = ggplot(tbl.m, aes(x=Var1, y=value)) 
  }
  
  if("variable" %in% colnames(tbl.m)) {
    if(class(tbl.m$variable)=="integer") tbl.m$variable = as.character(tbl.m$variable)
    ggp = ggplot(tbl.m, aes(x=variable, y=value)) 
  }
  
  if(!is.null(gg_theme)){
    ggp = ggp + gg_theme()
  } else {
    ggp = ggp + ggplot2::theme_bw()
  }
  
  
  ggp   +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggplot2::ggtitle(Title)
}



#' @title gg_tSNE_df
#'
#' @description Takes a data.frame with tSNE1 and tSNE2 precomputed and plots it.
#' @param DF A data.frame 
#' @param SavePath A data.frame 
#' @param Markers A char. vec. that that all must exist in the data.frame.
#' @param splitby a single character which data.frame must contain this value in the cols; the resulting plot is faceted by this.
#' @param verbose if T verbose output is printed 
#' @param Title A title else ''
#' @param gg_theme a ggplot theme or NULL
#' @return A ggplot 
#' @export
gg_tSNE_df <- function(DF, 
                       SavePath=NULL, 
                       Markers=NULL, 
                       splitby=NULL,
                       verbose = T,
                       size=0.2, alpha=0.6,
                       Title="",
                       width = 2000, height = 2300, pointsize = 20, ncol=4, 
                       limits = c(0, 1),
                       gg_theme=NULL, printPlot = T, vridis_col_opt = "A"){
  
  if(! ("tSNE1" %in% colnames(DF)) & "tSNE2" %in% colnames(DF)) stop("tSNE1 and/or tSNE2 not found in df")
  if(is.null(Markers)) stop("Markers is NULL")
  
  ExprPlots <- list()
  
  if(is.null(splitby)){
    
    for(xN in Markers){
      if(verbose) print(xN)
      

      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", xN)], id.vars=c("tSNE1", "tSNE2")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        ggthemes::scale_colour_gradient2_tableau('pi0', low = "blue", mid = "white", high = "red", midpoint = 0)
        scale_color_viridis(option=vridis_col_opt,limits = limits, oob = scales::squish) +
        ggtitle(paste0(Title, "\ntSNE - ", xN))
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  } else {
    
    if(! splitby %in% colnames(DF)) stop("splitby not in colnames of DF")
    
    DF$VarSplit = DF[,splitby]
    
    for(xN in Markers){
      # xN = Markers[1]
      if(verbose) print(xN)

      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", "VarSplit", xN)], id.vars=c("tSNE1", "tSNE2", "VarSplit")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        scale_color_viridis(option="A") +
        ggtitle(paste0(Title, "\ntSNE - ", xN, "\n per", splitby))+ facet_wrap(~VarSplit)
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  }
  
  
  if(!is.null(SavePath)) {
    
    png(SavePath, width = width, height = height, pointsize = pointsize)
    
    print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol))
    
    dev.off()
    
    
  } else { 
    if(printPlot) {
      return(print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol)))
    } else {
      return(ExprPlots)
    } 
    
    }
  
}