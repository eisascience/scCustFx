#' @title Run UMAP dimensionality reduction on a matrix using umap-learn which needs python
#' @param DGEmat matrix of gene expression data
#' @param assay character name for the assay
#' @param n.neighbors integer, number of neighbors to use
#' @param n.components integer, number of dimensions to keep
#' @param metric character, metric to use for distance calculation
#' @param n.epochs integer, number of epochs for optimization
#' @param learning.rate numeric, learning rate for optimization
#' @param min.dist numeric, minimum distance between points
#' @param spread numeric, spread of clusters
#' @param set.op.mix.ratio numeric, mix ratio for set operation
#' @param local.connectivity integer, controls local connectivity
#' @param repulsion.strength numeric, repulsion strength
#' @param negative.sample.rate integer, number of negative samples
#' @param a numeric, parameter for umap
#' @param b numeric, parameter for umap
#' @param seed.use integer, seed for reproducibility
#' @param metric.kwds list, additional arguments for metric
#' @param angular.rp.forest logical, use angular random projection forest
#' @param reduction.key character, prefix for reduction result names
#' @param verbose logical, verbosity flag
#' @return matrix of reduced dimensions
#' @export
RunUMAP.Matrix.py <- function(
  #originally from Seurat pacakge,
  DGEmat,
  assay = NULL,
  n.neighbors = 30L,
  n.components = 2L,
  metric = "correlation",
  n.epochs = NULL,
  learning.rate = 1.0,
  min.dist = 0.3,
  spread = 1.0,
  set.op.mix.ratio = 1.0,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5,
  a = NULL,
  b = NULL,
  seed.use = 42,
  metric.kwds = NULL,
  angular.rp.forest = FALSE,
  reduction.key = 'UMAP_',
  verbose = TRUE,
  ...
) {
  require(reticulate)

  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }

  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  if (typeof(x = n.epochs) == "double") {
    n.epochs <- as.integer(x = n.epochs)
  }
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    metric_kwds = metric.kwds,
    angular_rp_forest = angular.rp.forest,
    verbose = verbose
  )
  umap_output <- umap$fit_transform(as.matrix(x = DGEmat))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- rownames(DGEmat)

  return(umap_output)
}




#' Perform UMAP Dimensionality Reduction on a Matrix
#'
#' This function applies the Uniform Manifold Approximation and Projection (UMAP) algorithm 
#' to reduce the dimensionality of a given data matrix using the native R `uwot` package.
#' 
#' @param DGEmat A numeric matrix where rows represent samples and columns represent features.
#' @param assay Optional assay information (not used in this implementation).
#' @param n.neighbors Integer. Number of nearest neighbors to consider (default: 30).
#' @param n.components Integer. Number of dimensions for the reduced space (default: 2).
#' @param metric Character. Distance metric to use (default: "correlation"). Options include `"euclidean"`, `"cosine"`, `"manhattan"`, etc.
#' @param min.dist Numeric. Minimum distance between embedded points (default: 0.3).
#' @param spread Numeric. Spread of the UMAP embedding (default: 1).
#' @param seed.use Integer. Seed for reproducibility (default: 42).
#' @param reduction.key Character. Prefix for naming UMAP dimensions in the output (default: "UMAP_").
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#' @param ... Additional arguments passed to `uwot::umap()`.
#'
#' @return A matrix of reduced UMAP coordinates with row names corresponding to `DGEmat`.
#'
#' @importFrom uwot umap
#'
#' @examples
#' \dontrun{
#' if(requireNamespace("uwot", quietly = TRUE)) {
#'     data_matrix <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#'     umap_results <- RunUMAP.Matrix(DGEmat = data_matrix)
#' }
#' }
#' @export
RunUMAP.Matrix <- function(DGEmat, assay = NULL, n.neighbors = 30L, n.components = 2L, 
                           metric = "correlation", min.dist = 0.3, spread = 1, 
                           seed.use = 42, reduction.key = "UMAP_", verbose = TRUE, ...) {
  
  # Load required package
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required but not installed. Install it using install.packages('uwot').")
  }
  
  # Set seed for reproducibility
  if (!is.null(seed.use)) {
    set.seed(seed.use)
  }
  
  # Run UMAP using the 'uwot' package
  umap_output <- uwot::umap(as.matrix(DGEmat),
                            n_neighbors = as.integer(n.neighbors),
                            n_components = as.integer(n.components),
                            metric = metric,
                            min_dist = min.dist,
                            spread = spread,
                            verbose = verbose)
  
  # Assign column and row names
  colnames(umap_output) <- paste0(reduction.key, 1:ncol(umap_output))
  rownames(umap_output) <- rownames(DGEmat)
  
  return(umap_output)
}
