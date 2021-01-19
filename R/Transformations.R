#' @title Range01
#'
#' @description Scales the range of input to be between 0 and 1
#' @param x a numerical vector
#' @param MaxN set Max value; NULL default find max(x)
#' @param MinN set Min value; NULL default find min(x)
#' @return A numeric vector
#' @keywords range
range01 <- function(x, MaxN = NULL, MinN = NULL){
  if(is.null(MaxN)) MaxN = max(x)
  if(is.null(MinN)) MinN = min(x)
  
  (x - MinN)/(MaxN - MinN)
}

#' @title QuantileEdgeClean
#'
#' @description Scales the range of input to be between 0 and 1
#' @param mat a numerical matrix, rows are features cols are cells/measured units
#' @param LowSet a numerical value for lower send quantiles; default 0
#' @param HighSet a numerical value for higher send quantiles; default NULL sets to high quantile
#' @param lowQuant a value to set the low quantiles
#' @param highQuant a value to set the high quantiles
#' @return A numeric mat
QuantileEdgeClean <- function(mat, LowSet=0, HighSet=NULL, lowQuant = 0.1, highQuant=0.9){
  
  mat = apply(mat, 1, function(x) {
    x[x<quantile(asinh(x), lowQuant)] = LowSet
    if(is.null(HighSet)) HS = quantile(asinh(x), highQuant)
    x[x>quantile(asinh(x), highQuant)] = HS
    x
  }) %>% t()
  
}