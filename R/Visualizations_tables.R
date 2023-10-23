

#' gg_pie_table
#'
#' This function creates a pie chart with labels and frequency values using ggplot2.
#'
#' @param TBL A table or data frame containing the data for the pie chart.
#' @param cols A vector of colors to be used for filling the pie slices.
#' @param legend.position The position of the legend in the plot.
#'
#' @export
gg_pie_table <- function(TBL = NULL, cols = col_vector, legend.position = "none"){
  FreqtempDF = as.data.frame.table(TBL) #%>% tibble::column_to_rownames("Var1") 
  FreqtempDF$fraction = round(FreqtempDF$Freq/sum(FreqtempDF$Freq), 4)
  FreqtempDF$ymax = cumsum(FreqtempDF$fraction)
  FreqtempDF$ymin = c(0, head(FreqtempDF$ymax, n=-1))
  
  # Compute label position
  FreqtempDF$labelPosition <- (FreqtempDF$ymax + FreqtempDF$ymin) / 2
  
  # Compute a good label
  FreqtempDF$label <- paste0(FreqtempDF$Var1, "\n value: ", FreqtempDF$Freq)
  
  ggplot(FreqtempDF, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) +
    geom_rect() +
    geom_label( x=4.7, aes(y=labelPosition, label=label), size=3) +
    scale_fill_manual(values = cols)+
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    theme_void() +
    theme(legend.position =legend.position)
  
  
}



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

#' @title gg_barplot_table, aka OOSAP::PlotMyTable
#'
#' @description Takes a table and plots a gg barplot 
#' @param MyTable A 2way table (base R table()) 
#' @param Title A title else ''
#' @param legend.position legend position defautl "bottom"
#' @param PlotCombo if T plots % and count plots together
#' @param xlab x axis label
#' @param xtext_angle x axis label
#' @param hjust if angle is 90 set to 1 
#' @param vjust if angle is 90 set to 0.5
#' @param returnlist logical
#' @param col_vector a vector of colors
#' @param theme a ggplot theme, default theme_bw(base_size = base_size)
#' @param base_size numerical
#' @return A ggplot 
#' @export
gg_barplot_2waytable = function(MyTable, Title="", legend.position="bottom", PlotCombo = F, xlab="", 
                                xtext_angle = 45, hjust = 0.9,vjust = 0.8,base_size = 10,
                                returnlist = T, col_vector=col_vector, 
                                theme_b = theme_bw(base_size = base_size)){
  tempMeltDF <- reshape2::melt(MyTable)
  
  
  p1 <- (ggplot(tempMeltDF) +
           geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
           theme_b + scale_fill_manual(values=col_vector) +
           theme(legend.position=legend.position,
                 legend.direction="horizontal",
                 axis.text.x = ggplot2::element_text(angle = xtext_angle, vjust = vjust, hjust=hjust),
                 legend.title = element_blank()) +
           ggtitle(paste0(Title, "\n Total Contribution")) + ylab("Total No.")) + xlab(xlab)
  
  
  p2<- (ggplot(tempMeltDF) +
          geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
          theme_b  + scale_fill_manual(values=col_vector) +
          theme(legend.position=legend.position,
                legend.direction="horizontal",
                axis.text.x = ggplot2::element_text(angle = xtext_angle, vjust = vjust, hjust=hjust),
                legend.title = element_blank())+
          scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                             labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
          ggtitle(paste0(Title, "\n Relative % Contribution")) + ylab("(%)")) + xlab(xlab)
  
  
  
  
  if(PlotCombo) {
    print(cowplot::plot_grid(p1, p2, ncol = 1)) 
    p1 + theme(legend.position="none") / p2
  } else {
    if(returnlist) {
      list(p1 = p1, p2=p2)
    } else {
      p1
      p2
    }
    
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



