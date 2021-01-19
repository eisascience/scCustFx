
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