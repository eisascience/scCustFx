#' Analyze Enrichment in a Contingency Table
#'
#' This function performs a chi-square test on a given contingency table and
#' generates various diagnostic plots to assess the association between the table's
#' dimensions. Optional plots include a balloon plot, mosaic plot, association plot,
#' heatmap of Pearson's residuals, correlation plot of Pearson's residuals, and a plot
#' displaying the relative contribution of Pearson's residuals.
#'
#' @param tempTbl A contingency table (e.g., produced by \code{\link{table}}) containing counts.
#' @param plotBalloon Logical; if \code{TRUE}, a balloon plot is generated using the
#'   \code{gplots} package.
#' @param plotMosaic Logical; if \code{TRUE}, a mosaic plot is generated.
#' @param plotAssoc Logical; if \code{TRUE}, an association plot is generated using the
#'   \code{vcd} package.
#' @param plotResHeatmap Logical; if \code{TRUE}, a heatmap of Pearson's residuals is produced
#'   using the \code{pheatmap} package. Default is \code{TRUE}.
#' @param plotCorrplot Logical; if \code{TRUE}, a correlation plot of Pearson's residuals is generated
#'   using the \code{corrplot} package.
#' @param plotContrib Logical; if \code{TRUE}, a plot showing the relative contribution of Pearson's residuals
#'   is generated using the \code{corrplot} package. Default is \code{TRUE}.
#' @param plotContrib_dir Logical; if \code{TRUE}, a plot showing the directional contribution of Pearson's residuals
#'   is generated using the \code{corrplot} package. Default is \code{TRUE}. The directional contribution retains the sign of the residuals.
#' @param Title Character; a title string to be used in the generated plots.
#'
#' @return A list containing:
#' \item{test}{The chi-square test object.}
#' \item{observed}{Observed counts from the contingency table.}
#' \item{expected}{Expected counts, rounded to two decimal places.}
#' \item{residuals}{Pearson's residuals, rounded to three decimal places.}
#' \item{p.value}{The p-value of the chi-square test.}
#' \item{contrib}{Relative contribution of Pearson's residuals (in percent), rounded to three decimal places.}
#'
#' @details This function requires the following packages: \code{gplots}, \code{pheatmap},
#'   \code{vcd}, and \code{corrplot}. Ensure these packages are installed prior to using the function.
#'
#' @examples
#' \dontrun{
#' # Create an example contingency table
#' tempTbl <- table(sample(letters[1:3], 100, replace = TRUE),
#'                  sample(c("Low", "High"), 100, replace = TRUE))
#'
#' # Analyze the table with selected plots
#' result <- analyzeTable_Enrichment(tempTbl,
#'                                   plotBalloon = TRUE,
#'                                   plotMosaic = TRUE,
#'                                   plotResHeatmap = TRUE,
#'                                   plotContrib = TRUE,
#'                                   Title = "Example Plot Title")
#'
#' # Print chi-square test results and associated metrics
#' print(result)
#' }
#'
#' @export
analyzeTable_Enrichment <- function(tempTbl,
                         plotBalloon = FALSE,
                         plotMosaic = FALSE,
                         plotAssoc = FALSE,
                         plotResHeatmap = T,
                         plotCorrplot = FALSE,
                         plotContrib = F, 
                         plotContrib_dir = T, 
                         Title="") {
  
  # Load required packages
  if (!requireNamespace("gplots", quietly = TRUE)) {
    stop("Please install the 'gplots' package.")
  }
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install the 'pheatmap' package.")
  }
  if (!requireNamespace("vcd", quietly = TRUE)) {
    stop("Please install the 'vcd' package.")
  }
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("Please install the 'corrplot' package.")
  }
  
  # Compute Chi-square test and associated values
  tempTbl.chisq <- chisq.test(tempTbl)
  
   
  
  chisqOutput <- list(
    test = tempTbl.chisq,
    observed = tempTbl.chisq$observed,
    expected = round(tempTbl.chisq$expected, 2),
    residuals = round(tempTbl.chisq$residuals, 3),
    p.value = tempTbl.chisq$p.value,
    contrib = round(100 * (tempTbl.chisq$residuals^2) / tempTbl.chisq$statistic, 3),
    contrib_dir = round(100 * tempTbl.chisq$residuals * abs(tempTbl.chisq$residuals) / tempTbl.chisq$statistic, 3)
  )
  
  # Plot a balloon plot if specified
  if (plotBalloon) {
    gplots::balloonplot(t(tempTbl), 
                        main = Title, 
                        xlab = "", ylab = "",
                        label = TRUE, show.margins = TRUE)
  }
  
  # Plot a mosaic plot if specified
  if (plotMosaic) {
    mosaicplot(tempTbl, shade = TRUE, las = 2,
               main = paste0(Title))
  }
  
  # Plot an association plot if specified
  if (plotAssoc) {
    vcd::assoc(tempTbl, shade = TRUE, las = 3, rot_labels = 45)
  }
  
  # Plot a heatmap of Pearson's residuals if specified
  if (plotResHeatmap) {
    pheatmap::pheatmap(round(tempTbl.chisq$residuals, 3), 
                       main = paste0("Pearson's chisqr residual \n", Title), 
                       cluster_cols = FALSE)
  }
  
  # Plot a correlation plot of Pearson's residuals if specified
  if (plotCorrplot) {
    corrplot::corrplot(tempTbl.chisq$residuals, is.cor = FALSE, 
                       title = paste0("Pearson's chisqr residual \n", Title))
  }
  
  # Plot the relative contribution of Pearson's residuals if specified
  if (plotContrib) {
    
    corrplot::corrplot(chisqOutput$contrib, is.cor = FALSE, 
                       title = paste0("Relative Contribution of residual \n", Title))
  }
  
  # Plot the relative contribution of Pearson's residuals if specified, but with direction
  if (plotContrib_dir) {
    
    corrplot::corrplot(chisqOutput$contrib_dir, is.cor = FALSE, 
                       title =  paste0("Directional Relative Contribution of residual \n", Title))
  }
  
  # Return a list with key Chi-square test results
  return(chisqOutput)
}


#' Get Top Correlation Pairs from a Correlation Matrix
#'
#' Extracts the top positively and negatively correlated variable pairs from a correlation matrix.
#'
#' @param corMat A square numeric matrix representing correlations between variables. Should have row and column names.
#' @param nTop An integer specifying the number of top positive and negative correlation pairs to return. Defaults to 5.
#'
#' @return A list with two data frames:
#' \describe{
#'   \item{topPositive}{Data frame of the top \code{nTop} most positively correlated variable pairs.}
#'   \item{topNegative}{Data frame of the top \code{nTop} most negatively correlated variable pairs.}
#' }
#'
#' @examples
#' # Generate a random correlation matrix
#' set.seed(123)
#' mat <- matrix(rnorm(100), ncol = 10)
#' corMat <- cor(mat)
#' result <- getTopCorrelationPairs(corMat, nTop = 3)
#' result$topPositive
#' result$topNegative
#'
#' @export
getTopCorrelationPairs <- function(corMat, nTop = 5) {
  # Check that corMat is a square matrix
  if (!is.matrix(corMat) || nrow(corMat) != ncol(corMat)) {
    stop("Input must be a square correlation matrix.")
  }
  # Get row and column names
  rn <- rownames(corMat)
  cn <- colnames(corMat)
  if (is.null(rn) || is.null(cn)) {
    rn <- paste0("V", seq_len(nrow(corMat)))
    cn <- rn
  }
  # Create a data frame of all pairs (i, j) with i < j
  pairs <- data.frame()
  for (i in 1:(nrow(corMat) - 1)) {
    for (j in (i + 1):ncol(corMat)) {
      pairs <- rbind(pairs, data.frame(
        Var1 = rn[i],
        Var2 = cn[j],
        Correlation = corMat[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }
  # Order by correlation value for positive and negative correlations
  topPositive <- pairs[order(-pairs$Correlation), ]
  topNegative <- pairs[order(pairs$Correlation), ]
  # Select top n pairs from each group
  topPositive <- head(topPositive, nTop)
  topNegative <- head(topNegative, nTop)
  return(list(
    topPositive = topPositive,
    topNegative = topNegative
  ))
}