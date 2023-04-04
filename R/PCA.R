
#' Principal Component Analysis (PCA) Plot
#'
#' This function generates a PCA plot using the prcomp function in R.
#'
#' @param dfx a numeric matrix containing the data to be plotted
#' @param labels a character vector containing the labels for each sample in \code{dfx}
#' @param factorV a factor or character vector indicating the group each sample belongs to
#' @param title a character string specifying the title of the plot
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance
#' @param center a logical value indicating whether the variables should be centered
#' @param col_vector a vector of color names to be used for plotting the samples
#' @param namePointBL a logical value indicating whether to include labels for each point
#'
#' @return If \code{returnPCA} is TRUE, the function returns a list containing the PCA results. Otherwise, the function generates a PCA plot.
#'
#' @examples
#' data("iris")
#' PCAmyDF(iris[,1:4], iris$Species, iris$Species, "Iris dataset", TRUE, TRUE, c("red", "blue", "green"), TRUE)
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @export
PCAmyDF <- function (dfx, labels, factorV, title = "PCA Plot", scale, center, col_vector, namePointBL = F) {
  if(class(labels) == "function") {
    print("no labels, using factor as names")
    labels = as.character(factorV)
  }
  if(length(as.character(factorV)) != length(labels)) {
    print("labels length != factor length, using factor as names")
    labels = as.character(factorV)
  }
  
  dfx.pca <- prcomp(t(dfx), scale.=scale, center = center)
  
  MinXY <- min(c(round(min(dfx.pca$rotation[,1:2]) - abs(min(dfx.pca$rotation[,1:2])*0.5),1), -1) )
  MaxXY <- max(c(round(max(dfx.pca$rotation[,1:2]) + abs(max(dfx.pca$rotation[,1:2])*0.5),1),  1) )
  
  
  if(namePointBL){
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      scale_color_manual(values = col_vector, name="Samples") +
      geom_text_repel(aes(y = PC2, label = labels))  +
      ggtitle(paste("PCA of:",title ,sep=" ")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  } else {
    autoplot(dfx.pca) +
      theme_bw() +
      geom_point(aes(color = factorV), alpha = 0.8, size = 3) +
      scale_color_manual(values = col_vector, name="Samples")  +
      ggtitle(paste("PCA of:",title ,sep=" "))+
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(MinXY,MaxXY)) + ylim(c(MinXY,MaxXY))
  }
  
  
  
  
}