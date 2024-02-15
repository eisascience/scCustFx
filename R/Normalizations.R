
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






#' @title dropsimNorm
#'
#' @description passes to dropsim::normaliseDGE() fx
#' @param GeneExprMat A matrix rows of genes and cols of cells
#' @return A quantile-based normalized matrix
#' @export
dropsimNorm <- function(GeneExprMat, featuresToUse, minLibrarySize = 0, maxLibrarySize = 10, center = FALSE, scale = TRUE, gene_subset = 1){
  
  
  normedDGE <- dropsim::normaliseDGE(Matrix::as.matrix(GeneExprMat[featuresToUse, 
  ]), center = center, scale = scale, threshold = maxLibrarySize, min_library_size = minLibrarySize, 
  gene_subset = gene_subset)
  
  
  # normedDGE <- as.matrix(normedDGE)
  
  
  
  return(normedDGE)
  
}

#' @title dropsimNormV2 speeds up using sparse matrix fx as much as possible
#'
#' @description  calls normaliseDGEv2 instead of dropsim::normaliseDGE()
#' @param GeneExprMat A matrix rows of genes and cols of cells
#' @return A quantile-based normalized matrix
#' @export
dropsimNormV2 <- function(GeneExprMat, featuresToUse, minLibrarySize = 0, maxLibrarySize = 10, center = FALSE, scale = TRUE, gene_subset = 1){
  
  
  normedDGE <- normaliseDGEv2(Matrix::as.matrix(GeneExprMat[featuresToUse, 
  ]), center = center, scale = scale, threshold = maxLibrarySize, min_library_size = minLibrarySize, 
  gene_subset = gene_subset)
  
  
  # normedDGE <- as.matrix(normedDGE)
  
  
  
  return((normedDGE))
  
}


#' @title normaliseDGEv2
#'
#' @description a sparse implementation of dropsim::normaliseDGE()
#' @param dge A matrix rows of genes and cols of cells
#' @return A quantile-based normalized matrix
#' @export
normaliseDGEv2 <-function (dge, verbose = FALSE, center = TRUE, scale = TRUE, 
          outliers = FALSE, threshold = 10, min_library_size = 0, 
          gene_subset = 0.8, transformation = "sqrt") 
{
  cells = colnames(dge)
  
  library_size <- Matrix::colSums(dge)
  lib_size_correction_factor <- library_size/median(sqrt(library_size))
  
  
  # Normalize each column by its library size correction factor using diagonal matrix multiplication
  #instead of # dge <- Matrix::t(Matrix::t(dge)/lib_size_correction_factor)
  correction_matrix <- Diagonal(x = 1 / lib_size_correction_factor)
  dge <- dge %*% correction_matrix
  
  colnames(dge) = cells
  
  cell_subset <- library_size >= min_library_size
  gene_mean <-  Matrix::rowMeans(dge[, cell_subset])
  gene_subset <- data.table(gene_mean, names(gene_mean))[order(-gene_mean)][1:round(length(gene_mean) * 
                                                                                      gene_subset)]$V2
  dge <- dge[gene_subset, cell_subset]
  
  
  if (transformation == "asinh") {
    dge <- asinh(10000 * dge)
  }
  else {
    dge <- sqrt(10000 * dge)
  }
  dge <-  Matrix::t(dge)
  
  
 
  
  if (verbose) {
    
    sample_dge <- function(dge) {
      if (prod(dim(dge)) > 1e+06) {
        tmp <- sample(dge, 1e+06)
      }
      else if (prod(dim(dge)) < 1e+06 & class(dge) == "dgCMatrix") {
        tmp <- as.matrix(dge)
      }
      else {
        tmp <- dge
      }
      return(tmp)
    }
    
    
    str(dge)
    
    hist(sample_dge(dge), breaks = 2000, ylim = c(0, 1000), 
         xlim = c(0, 50), xlab = "Expression", main = "Histogram of data post cell-wise (& sqrt) normalisation")
    
  }
  
  if (center == TRUE) {
    dge <- scale(dge, center = TRUE, scale = scale)
  } else {
    if(scale){
      dge <- Matrix::t(Matrix::t(dge)/colSdColMeans(dge))
    }
  }
  
  dge <- dge[, !is.na(colSums(dge))]
  
  if (verbose) {
    hist(sample_dge(dge), breaks = 2000, ylim = c(0, 10000), 
         main = "Histogram of data post gene-wise normalisation")
  }
  
  stopifnot(sum(is.na(dge)) == 0)
  
  if (outliers) {
    max_by_cell <- apply(dge, 1, max)
    max_by_gene <- apply(dge, 2, max)
    plot(max_by_gene, main = "Maximum expression value per gene", 
         ylab = "Maximum Expression", xlab = "Gene Index (by mean expression)")
    abline(h = threshold, col = "red")
    plot(max_by_cell[], main = "Maximum expression value per cell", 
         ylab = "Maximum Expression", xlab = "Cell Index (by library size)")
    abline(h = threshold, col = "red")
  }
  
  
  sum(dge > threshold)
  dge[dge > threshold] <- threshold
  dge[dge < (-threshold)] <- (-threshold)
  
  
  if (verbose) {
    hist(sample_dge(dge), breaks = 500, ylim = c(0, 1000), 
         main = "Histogram of final data", xlab = "Expression")
    str(dge)
  }
  
  gc()
  
  return(dge)
}

#' @title colSdColMeans
#' 
#' Calculate Column Standard Deviations
#'
#' Computes the standard deviation of each column in a matrix, optionally removing NA values.
#' This function is optimized for matrices, including sparse matrices from the Matrix package,
#' ensuring efficient computation without densifying the sparse matrix. The standard deviation
#' is computed in a way that is mindful of potential NA values, based on the user's choice.
#'
#' @description a sparse implementation of dropsim::normaliseDGE()
#' @param dge A matrix rows of genes and cols of cells
#' @return A quantile-based normalized matrix
#' @export
colSdColMeans <- function(x, na.rm = TRUE) {
  if (na.rm) {
    n <-  Matrix::colSums(!is.na(x))
  }
  else {
    n <- nrow(x)
  }
  colVar <-  Matrix::colMeans(x * x, na.rm = na.rm) - ( Matrix::colMeans(x, 
                                                                         na.rm = na.rm))^2
  return(sqrt(colVar * n/(n - 1)))
}
