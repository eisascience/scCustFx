
#' @title quantile_norm
#'
#' @description quantile normalize a matrix
#' @param mat A matrix
#' @return A quantile-based normalized matrix
#' @export
quantile_norm <- function(mat){
  #ref/cite online source:
  mat_rank <- apply(mat,2,rank,ties.method="min")
  mat_sorted <- data.frame(apply(mat, 2, sort))
  mat_mean <- apply(mat_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  mat_final <- apply(mat_rank, 2, index_to_mean, my_mean=mat_mean)
  rownames(mat_final) <- rownames(mat)
  return(t(mat_final))
}

#' @title avgExpr_norm
#'
#' @description Normalize a matrix by average expression (total and or per feature)
#' @param mat A matrix where the features are the rows, cells/measured units are columns
#' @return A normalized matrix
#' @export
avgExpr_norm <- function(mat, doGlobal=T, doPerFeat = T, rm.NA = T){

  if(doGlobal) mat = mat/mean(rowSums(mat))
  if(doPerFeat) mat = (mat)/rowMeans(mat) 
  if(length(is.na(mat))>0) warning("NAs are found in matrix")
  if(rm.NA) mat[is.na(mat)]=0
  
}

#' Normalize multiple matrices using quantile normalization
#'
#' @param matrices A list of matrices
#' @return A list of normalized sparse matrices
normalizeMultipleMatrices <- function(matrices) {
  
  require(preprocessCore)
  
  print("preProcessing")
  # Concatenate matrices along columns
  matrices_concatenated <- do.call(cbind, matrices)
  
  print("starting normalization")
  # Perform quantile normalization
  normalized_matrices_concatenated <- normalize.quantiles(matrices_concatenated)
  
  print("postProcessing")
  # Get the column index for each matrix
  indexes <- cumsum(sapply(matrices, ncol))
  
  # Initialize a list to store the normalized matrices
  normalized_matrices <- vector("list", length(matrices))
  
  
  # Loop over the index and create matrices
  for (i in 1:(length(matrices))) {
    if (i==1) {
      normalized_matrices[[i]] <- as(normalized_matrices_concatenated[,1:indexes[i]], "sparseMatrix")
    } else {
      normalized_matrices[[i]] <- as(normalized_matrices_concatenated[,(indexes[i-1]+1):indexes[i]],"sparseMatrix")
    }
  }
  
  names(normalized_matrices) = names(matrices)
  return(normalized_matrices)
}

