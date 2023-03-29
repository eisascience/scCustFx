#' @title Run UMAP dimensionality reduction on a matrix
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
RunUMAP.Matrix <- function(
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
