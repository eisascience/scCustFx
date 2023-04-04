

#' MDS Plot of Data Frame
#'
#' Performs multidimensional scaling (MDS) on a data frame and generates a plot of the
#' results. The function uses the ggplot2 package to create a scatterplot of the MDS
#' coordinates with colored points based on a factor variable. The ggrepel package is used
#' to add labels to the points without overlap.
#'
#' @param dfx A data frame with samples as rows and variables as columns.
#' @param labelsDF A data frame with the labels for each sample in the same order as dfx.
#' @param factorV A vector containing the factor variable to use for coloring the points.
#' @param title A character string for the title of the plot.
#' @param col_vector A vector of colors to use for coloring the points.
#' @param returnMDS A logical value indicating whether to return the MDS coordinates as a data frame.
#'   Default is FALSE.
#' @return If returnMDS is TRUE, returns a data frame of the MDS coordinates. Otherwise, returns
#'   a ggplot object.
#' @examples
#' data(iris)
#' MDSmyDF(dfx = iris[, 1:4], labelsDF = iris[, 5], factorV = iris[, 5], title = "Iris Data MDS Plot", 
#'   col_vector = c("#FF0000", "#00FF00", "#0000FF"))
#' @import ggplot2
#' @import ggrepel
#' @importFrom stats cmdscale dist
#' @export
MDSmyDF <- function(dfx, labelsDF, factorV, title = "MDS Plot", col_vector, returnMDS = F){
  
  
  if(length(factorV) != nrow(dfx)) stop("rotate t() your dataX")
  
  mds <- cmdscale(as.matrix(dist(dfx)))
  colnames(mds) <- c("MDS1", "MDS2")
  
  mds <- cbind(mds, labelsDF) #samples as rows add a cols of labels
  
  
  p1 <- ggplot(mds, aes(x = MDS1, y = MDS2)) +
    theme_bw() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
    coord_cartesian(xlim = c(min(mds[,1])-5, max(mds[,1])+5)) +
    scale_color_manual(values = col_vector, name="Samples")
  
  # the graphic with ggrepel
  p1 <- p1 + geom_text_repel(aes(y = MDS2 + 0.25), label = factorV) +
    ggtitle(paste("MDS of:",title ,sep=" "))+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  if(returnMDS) return(mds) else p1
  
  
}