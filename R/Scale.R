
#' Convert Vector to 0-1 Range
#'
#' This function converts a numeric vector to a range of 0 to 1 by linearly scaling the values
#' between the minimum and maximum of the input vector.
#'
#' @param x Numeric vector to be scaled.
#' @param MaxN Numeric value representing the maximum limit of scaling. If NULL, defaults to the maximum value of x.
#' @param MinN Numeric value representing the minimum limit of scaling. If NULL, defaults to the minimum value of x.
#' @return Numeric vector scaled to range of 0 to 1.
#' @examples
#' x <- c(2, 5, 9, 11, 14)
#' range01(x)
#' @export
range01 <- function(x, MaxN = NULL, MinN = NULL){
  if(is.null(MaxN)) MaxN = max(x)
  if(is.null(MinN)) MinN = min(x)
  
  (x - MinN)/(MaxN - MinN)
}


