


#' Simulate 3D Data
#'
#' The function generates 3D data with N individuals, L variables and P replicates.
#'
#' @param N integer, the number of individuals.
#' @param L integer, the number of variables.
#' @param P integer, the number of replicates.
#' @param ncomps integer, the number of components.
#' @return data list, a list containing the generated data including Y, A, B, X, and noise.
simulate_3D_data <- function(N, L, P, ncomps) {
  A <- matrix(rnorm(N*ncomps), N, ncomps) # individual scores
  B <- matrix(sample(c(-1,1), P*ncomps, replace=TRUE), 
              P, ncomps) # individual scores
  
  X <- matrix(rnorm(ncomps*L), ncomps, L) # loadings
  X[sample(1:length(X), 0.9*length(X))] <- 0  # make loadings sparse
  
  noise <- array(rnorm(N*L*P, 0, 1), c(N, L, P))
  
  Y <- array(NA, c(N, L, P))
  for (p in 1:P) {
    Btmp <- matrix(B[p,], N, ncomps, byrow=TRUE)
    Y[,,p] <- (A * Btmp) %*% X + noise[,,p]
  }
  
  data <- list()
  data$Y <- Y
  data$A <- A
  data$B <- B
  data$X <- X
  data$noise <- noise
  
  return(data)
}

