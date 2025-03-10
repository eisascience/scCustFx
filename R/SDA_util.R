
#' Preprocess Raw Count Matrix for SDA
#'
#' This function preprocesses a raw gene expression count matrix by performing basic quality control checks,
#' computing library sizes, and visualizing feature distributions per cell/sample. It returns a list containing 
#' the filtered expression matrix and selected features.
#'
#' @param rawCountMat A numeric matrix or sparse matrix where rows represent features (genes) and columns represent cells/samples.
#' @param minLibrarySize Integer. Minimum library size threshold for filtering (default: 1).
#'
#' @return A list containing:
#' \item{GEX}{Preprocessed count matrix (same as input for now).}
#' \item{featuresToUse}{Vector of feature (gene) names included in preprocessing.}
#' \item{minLibrarySize}{The minimum library size threshold used.}
#'
#' @details The function checks if the number of cells/samples exceeds recommended limits for SDA (Sparse Decomposition of Arrays).
#' It then calculates total features per cell, generates a density plot with peak estimation, and returns the input data along with feature metadata.
#' 
#' @importFrom Matrix colSums
#' @importFrom stats density
#' @importFrom ggplot2 ggplot aes geom_density scale_x_sqrt geom_vline ggtitle labs theme_bw
#'
#' @export
#'
#' @examples
#' \dontrun{
#' if(requireNamespace("Matrix", quietly = TRUE)) {
#'     raw_counts <- Matrix::Matrix(rpois(1000, lambda = 10), nrow = 100, ncol = 10, sparse = TRUE)
#'     processed_data <- SDA_counts_preproc(raw_counts)
#'     str(processed_data)
#' }
#' }
SDA_counts_preproc <- function (rawCountMat, minLibrarySize = 1) 
{
  #counter part with Seurat object is SDA_ser_preproc
  
  # rawCountMat = count_data_EM;
  
  n_cells <- ncol(rawCountMat) #for bulk this is n_samples
  if (n_cells > 2e+05) {
    stop("SDA has shown to max handle ~200K cells ")
  } else if (n_cells > 150000) {
    warning("SDA has shown to max handle ~150K cells ")
  }
  print(paste0("Initial features: ", nrow(rawCountMat)))
  print(paste0("Initial cells: ", ncol(rawCountMat)))
  
  featuresToUse <- rownames(rawCountMat)
  
  df <- data.frame(x = Matrix::colSums(rawCountMat[featuresToUse, ]))
  
  dens <- stats::density(df$x)
  mx <- dens$x[which.max(dens$y)]
  
  ggplot(df, aes(x = x)) + scale_x_sqrt() + geom_density() + 
    ggtitle("Total features per cell") + labs(x = "Features/Cell", 
                                              y = "Density") + geom_vline(xintercept = mx, color = "red") + 
    ggtitle(paste0("Library Size: Peak = ", mx)) + theme_bw()
  
  print("starting dropsim normaliseDGE")
  return(list(GEX = rawCountMat, 
              featuresToUse = featuresToUse, 
              minLibrarySize = minLibrarySize))
}



#' Gene overlap across components
#'
#' Given a set of genes and a dataframe of gene expression values across multiple components, this function computes the number of overlapping genes in the top \code{Ngenes} of each component with the input gene set. 
#'
#' @param GeneSet a character vector of gene names to compare across the components
#' @param compDF a data.frame object containing gene expression values across components, where each column corresponds to a component and each row to a gene
#' @param Ngenes the number of genes to consider in the top overlap
#' @param verbose a boolean indicating whether to print the top found genes for each component
#'
#' @return A named list of length equal to the length of the input \code{GeneSet}, where each element of the list is a numeric vector of length equal to the number of components in \code{compDF} representing the number of overlapping genes in the top \code{Ngenes} of each component with the input gene set.
#' 
#'
#' @export
GeneOverlapAcrossComps  = function(GeneSet, compDF, Ngenes = 10, verbose = T){

  tempLS = lapply(GeneSet, function(xG){
    # xG = GeneSet[6]
    
    
    TopFound = lapply(compDF, function(gS){
      sum(any(xG %in% as.character(gS[1:Ngenes])))
    }) %>% unlist()
    
    
    if(verbose){
      
      TopFoundDF = as.data.frame(cbind(Top_Pos, Top_Neg)[1:10,names(TopFound[TopFound>0])]); TopFoundDF
      
      print(TopFoundDF)
    }
    
 
    comp2show = names(TopFound[TopFound>0]);
    
    
  })
  names(tempLS) = GeneSet
  
  return(tempLS)
}



#' Rotate SDA factorisation
#'
#' @param X SDA results list to be rotated
#' @param reference SDA results list to for X to be rotated towards
#' 
#' @return SDA results list, but with the gene loadings and scores matrices rotated by procrustes
#' 
#' @export
#' @import vegan
#' 
rotate_SDA <- function(X=results9_WT, reference=results1){
  
  if(ncol(X$loadings[[1]])>ncol(reference$loadings[[1]])){
    common_genes <- colnames(X$loadings[[1]])[colnames(X$loadings[[1]]) %in% colnames(reference$loadings[[1]])]
  }else{
    common_genes <- colnames(reference$loadings[[1]])[colnames(reference$loadings[[1]]) %in% colnames(X$loadings[[1]])]
  }
  
  rot_9WT <- vegan::procrustes(t(reference$loadings[[1]][,common_genes]), 
                               t(X$loadings[[1]][,common_genes]))
  
  colnames(rot_9WT$Yrot) <- paste0(rownames(X$loadings[[1]]),"rot")
  
  colnames(rot_9WT$rotation) <- rownames(reference$loadings[[1]])
  rownames(rot_9WT$rotation) <- rownames(X$loadings[[1]])
  
  # rotate scores too
  rot_9WT_scoresrot <- X$scores %*% rot_9WT$rotation
  colnames(rot_9WT_scoresrot) <- paste0(rownames(X$loadings[[1]]),"rot")
  str(rot_9WT_scoresrot)
  
  row_9WT_list <- list(loadings=list(t(rot_9WT$Yrot)), scores=rot_9WT_scoresrot, rotation=rot_9WT$rotation)  
  
  return(row_9WT_list)
}


#' @title Predict using SDA
#' @description Predict using sparse discriminant analysis
#' @param new.data a matrix containing new data to be predicted
#' @param sda.model a fitted sparse discriminant analysis model
#' @param return_cell_score a logical indicating whether to return the cell score or the predicted new data
#' @return If return_cell_score is T, returns cell score, else returns the predicted new data.
#' @export
predict.SDA = function(new.data, sda.model, return_cell_score=T) {
  # new.data = DGEnew ; sda.model = SDAres
  
  if(colnames(sda.model$scores)[1] != rownames(sda.model$loadings[[1]])[1]){
    colnames(sda.model$scores) = rownames(sda.model$loadings[[1]])
  }
  
  L = sda.model$loadings[[1]]
  
  
  commonGenes = intersect(colnames(new.data), colnames(L))
  
  if(length(colnames(new.data))/length(commonGenes)<1) print("genes needed to be removed")
  
  # CSnew = tcrossprod(Matrix::as.matrix(new.data[,commonGenes]), L[,commonGenes])
  
  CSnew = Matrix::as.matrix(new.data[,commonGenes]) %*% t( (L[,commonGenes]))
  
  
  if(return_cell_score) {
    return(CSnew) 
  } else {
    new.predicted = CSnew %*% L[,commonGenes]
    return(new.predicted)
  }
}
