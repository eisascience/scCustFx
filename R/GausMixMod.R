
#' Compute the intersection threshold for a two‐component Gaussian mixture
#'
#' Given the parameters of a two‐component Gaussian mixture model (means and variances),
#' this function solves for the point where the two Gaussian density curves intersect.
#'
#' @param params A list of mixture parameters, typically `model$parameters` from an \code{Mclust} object,
#'   containing:
#'   \itemize{
#'     \item \code{mean}: Numeric vector of component means.
#'     \item \code{variance$sigmasq}: Numeric vector of component variances.
#'   }
#' @param dataV A vector of the data used to make gmm
#' @return A single numeric value representing the threshold at which the two Gaussian
#'   components have equal density (the intersection point).
#' @examples
#' \dontrun{
#'   fit <- Mclust(my_data, G = 2)
#'   thresh <- compute_gmm_threshold(fit$parameters, dataV=my_data)
#' }
#' @export
compute_gmm_threshold <- function(params, dataV) {
  mu    <- params$mean
  sigma <- sqrt(params$variance$sigmasq)
  # Order by mean
  ord    <- order(mu)
  mu1    <- mu[ord[1]]; mu2 <- mu[ord[2]]
  s1     <- sigma[ord[1]]; s2  <- sigma[ord[2]]
  # Solve for x: dnorm(x,mu1,s1) = dnorm(x,mu2,s2)
  a <- 1/(2*s1^2) - 1/(2*s2^2)
  b <- mu2/(s2^2) - mu1/(s1^2)
  c <- mu1^2/(2*s1^2) - mu2^2/(2*s2^2) - log(s2/s1)
  roots <- polyroot(c(c, b, a))
  # Choose real root within data range
  real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
  thresh <- real_roots[real_roots > min(dataV, na.rm=TRUE) &
                         real_roots < max(dataV, na.rm=TRUE)]
  return(thresh[1])
}