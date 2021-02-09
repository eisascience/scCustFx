#' @title ColorTheme
#'
#' @description a vector of colors that are in sequence pretty and distinguishable
#' @return A vector of colors.
#' @export
ColorTheme <- function(){
  scaleyellowred <- grDevices::colorRampPalette(c("dodgerblue", "lightyellow", "red"), space = "rgb")(30)
  
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[-4]
  
  col_vector2 <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52",
                   "#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                   '#e6194b', '#3cb44b', '#ffe119', '#4363d8',
                   '#f58231', '#911eb4', '#46f0f0', '#f032e6',
                   '#bcf60c', '#fabebe', '#008080', '#e6beff',
                   '#9a6324', '#fffac8', '#800000', '#aaffc3',
                   '#808000', '#ffd8b1', '#000075', '#808080',
                   '#ffffff', '#000000')
  # c("#00AFBB", "#E7B800", "#FC4E07")
  
  col_vector2 <- as.character(sapply(1:20, function(xN){
    c(col_vector[xN], col_vector2[xN])
  }))
  col_vector  <- c(col_vector2, col_vector[21:length(col_vector)])
  
  col_vector <- col_vector[c(1:4, 6:length(col_vector))]
  
  return(list(col_vector=col_vector, scaleyellowred=scaleyellowred))
}



#' @title Cold2Hot
#'
#' @description a vector of colors that are in sequence well for a heatmap
#' @return A vector of colors.
#' @export
Cold2Hot <- function(plotDemo = F){
  library(wesanderson)
  
  Rbrew10_Spec = RColorBrewer::brewer.pal(10, "Spectral")
  Wes10_Ziss = as.character(wes_palette("Zissou1", 10, type = "continuous"))
  Rbrew10_PRGn = RColorBrewer::brewer.pal(10, "PRGn")
  
  col_vec = c(Wes10_Ziss, Rbrew10_Spec, Rbrew10_PRGn, c("black", "navy", "gold", "firebrick", "dodgerblue"))
  
  # barplot(rep(1, 35), col =col_vec,
  #         las=2, names.arg = 1:35)
  
  if(plotDemo) barplot(rep(1, 14), col =col_vec[c(32, 35, 19, 4, 18, 28, 15, 6, 13, 10, 11, 22, 21)])
  
  return(col_vec[c( 32, 35, 19, 4, 18, 28, 15, 6, 13, 10, 11, 22, 21)])
  
  
}

##TODO: make a function that takes parameters like color set and returns the scale_color_gradientn 

#' @title scale_color_gradient_heat1
#'
#' @description a vector of colors that are in sequence well for a heatmap
#' @return A vector of colors.
#' @export
scale_color_gradient_heat1 <- ggplot2::scale_color_gradientn(colours = Cold2Hot() )

#' @title scale_color_gradient_heat2
#'
#' @description a vector of colors that are in sequence well for a heatmap
#' @return A vector of colors.
#' @export
scale_color_gradient_heat2 <- ggplot2::scale_color_gradientn(colours = c("black", "navy", "blue", "dodgerblue", "gold", "red"))
