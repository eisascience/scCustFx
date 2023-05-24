
#' Perc2Radar
#'
#' This function creates a radar chart based on the percentage data provided in a data frame.
#'
#' @param PercDF A data frame containing percentage data.
#' @param Features A character vector specifying the features (columns) to be included in the radar chart.
#' @param col_vector A vector of colors to be used for plotting the radar chart.
#' @param plotLegend A logical value indicating whether to plot the legend.
#'
#' @export
Perc2Radar <- function(PercDF= NULL, Features=NULL, col_vector = col_vector, plotLegend=T){
  if(is.null(PercDF) | is.null(Features)) stop("bad input")
  
  tempDF_radar <- PercDF[, Features]
  
  tempDF_radar <- rbind(rep(100,ncol(tempDF_radar)) , rep(0,ncol(tempDF_radar)) , tempDF_radar)
  
  fmsb::radarchart(tempDF_radar, axistype=1, 
                   pcol = col_vector , 
                   plwd=4 , plty=1,
                   #custom the grid
                   cglcol="grey", cglty=1, axislabcol="grey", 
                   #caxislabels=seq(0,100,10), 
                   cglwd=0.8,
                   #custom labels
                   vlcex=0.8, 
                   calcex= 0.8
  )
  
  if(plotLegend){
    par(xpd=TRUE)
    legend(x="topright", legend = rownames(tempDF_radar[-c(1,2),]), 
           bty = "n", pch=20 , col=col_vector , 
           text.col = "black", cex=1.2, pt.cex=3)
    par(xpd=F)
  }
}

