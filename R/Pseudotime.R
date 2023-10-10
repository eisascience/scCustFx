
#' Principal Curve Bootstrap
#'
#' This function performs bootstrap estimation of principal curves based on the given input data.
#'
#' @param input A matrix or data frame representing the input data.
#' @param df_range A numeric vector specifying the range of degrees of freedom (df) for principal curve fitting.
#' If any element is set to NA, the range will be determined by df_start and df_end parameters.
#' @param df_start A numeric value indicating the starting degree of freedom for principal curve fitting.
#' Ignored if df_range is specified.
#' @param df_end A numeric value indicating the ending degree of freedom for principal curve fitting.
#' Ignored if df_range is specified.
#' @param plot A logical value indicating whether to plot the principal curves.
#' @return A list containing the principal curves fitted for different degrees of freedom.
#'
#' @import princurve
#' @export
princurve_bootstrap <- function (input, df_range = NA, df_start = NA, df_end = NA, plot = TRUE) 
{

  principal_curves <- list()
  if(anyNA(df_range)) df_range <- df_start:df_end else df_range <- sort(df_range)
  if(is.na(df_start)) df_start <- df_range[1]
  
  
  for (i in seq_along(df_range)) {
    # i = 1
    cat(paste("Fitting princurve with df", df_range[i], "\n"))
    
    if (df_range[i] == df_start) {
      principal_curves[[i]] <- princurve::principal_curve(input, df = df_range[i])
    } else {
      principal_curves[[i]] <- princurve::principal_curve(input, df = df_range[i], 
                                                          start = principal_curves[[i - 1]])
    }
    
    if (plot) plot(principal_curves[[i]])
  }
  names(principal_curves) <- paste0("df_", df_range)
  return(principal_curves)
}



#' Find Pseudotime
#'
#' This function finds the pseudotime corresponding to a given sample point.
#'
#' @param sample.point A numeric vector representing the sample point.
#' @param pseudotime A matrix or data frame containing pseudotime values.
#'
#' @return The index of the pseudotime corresponding to the closest match to the given sample point.
#'
#' @examples
#' sample.point <- c(0.5, 0.3, 0.7)
#' pseudotime <- matrix(c(0.1, 0.4, 0.8, 0.2, 0.3, 0.6), ncol = 3)
#' find_pseudotime(sample.point, pseudotime)
#'
#' @export
find_pseudotime <- function(sample.point, pseudotime){ which.min(colSums((t(pseudotime) - c(sample.point))^2)) }

