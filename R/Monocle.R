#' PseudoMySeur
#'
#' A function to create a pseudo-time trajectory using Monocle on a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param features A character vector specifying the features (genes) to include in the analysis. If NULL, all features are used.
#' @param removeCells A character vector specifying the cell names to remove from the analysis. If NULL, no cells are removed.
#' @param removeNoisyFeats A logical value indicating whether to remove noisy features. Default is TRUE.
#'
#' @return A Monocle CellDataSet object representing the pseudo-time trajectory.
#'
#' @export
PseudoMySeur <- function(SerObj, features=NULL, removeCells=NULL, removeNoisyFeats = T){
  
  if(is.null(features)) {
    features = rownames(SerObj)
  }
  
  if(removeNoisyFeats){
    features <- features[!grepl("RPL", features)]
    features <- features[!grepl("RPS", features)]
    features <- features[!grepl("^MT", features)]
    
  }
  
  features = features[features %in% rownames(SerObj)]  
  
  
  DGE_mat <- GetAssayData(SerObj, slot = "data")[features, ]
  
  metadata = SerObj@meta.data
  
  if(!is.null(removeCells)){
    DGE_mat = DGE_mat[,setdiff(colnames(DGE_mat), c(removeCells))]
  }
  
  
  
  metadata <- metadata[colnames(DGE_mat), ]
  
  library(monocle)
  
  pd <- new("AnnotatedDataFrame", data = metadata[colnames(DGE_mat), ])
  fd <- new("AnnotatedDataFrame", data = data.frame(varpos = rownames(DGE_mat), gene_short_name = rownames(DGE_mat)))
  rownames(fd) <- fd$varpos
  
  monocle <- newCellDataSet(scale(as.matrix(DGE_mat), center = F, scale = F),
                            phenoData = pd,
                            featureData = fd,
                            lowerDetectionLimit = 0,
                            expressionFamily = uninormal())
  fData(monocle)$varpos <- TRUE
  
  
  print("Reducing dimension by DDRTree...")
  monocle <- reduceDimension(monocle,
                             max_components = 15,  # attempted: 50, 20, 10, 15, 14
                             reduction_method = 'DDRTree',
                             norm_method = "none",
                             pseudo_expr = 0,
                             verbose = T )
  
  print("Ordering cells...")
  
  monocle <- orderCells(monocle)
  
  
  monocle$Pseudotime_cut <- cut(monocle$Pseudotime, 5, labels = paste0("cut", 1:5))
  monocle$State_Pseudo <- paste0("State", monocle$State, "_", monocle$Pseudotime_cut)
  
  return(monocle)
  
  
  
  
}



#' Ser2Monocle_MakeNProcess
#'
#' Performs the preprocessing steps for single-cell RNA-seq data using Monocle 2 package
#'
#' @param SeurObj A Seurat object with a count matrix and metadata
#' @param retunMon Logical value indicating whether to return the Monocle 2 object (default = TRUE)
#' @param PCAnDim Number of principal components to retain for dimensionality reduction (default = 20)
#' @param doUMAP Logical value indicating whether to perform UMAP reduction (default = TRUE)
#' @param min_dist The effective minimum distance between embedded points in the UMAP plot (default = 0.3)
#' @param n_neighbors The size of local neighborhood (default = 40)
#' @param dotSNE Logical value indicating whether to perform tSNE reduction (default = FALSE)
#' @param doClust Logical value indicating whether to perform clustering (default = TRUE)
#' @param ClusLouvRes Louvain resolution parameter for clustering (default = 0.00005)
#' @param louvain_iter The number of iterations for the Louvain algorithm (default = 3)
#' @param KeepTopNgenes The number of top variable genes to use in the analysis (default = 3000)
#' @param minExprGeneDet The minimum expression threshold for genes to be considered expressed (default = 0.1)
#' @param upperScale The scaling factor for the upper bound for filtering cells by total mRNA expression (default = 2)
#' @param lowerScale The scaling factor for the lower bound for filtering cells by total mRNA expression (default = 2.5)
#' @param assay The assay type to use (default = NULL)
#' @import monocle
#' 
#' @return A Monocle object
#'
#' @export
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- as.Seurat(pbmc_small)
#' SeurObj_cds <- Ser2Monocle_MakeNProcess(SeurObj = pbmc_small, retunMon = T,
#' PCAnDim = 20, doUMAP = T, min_dist=.3,
#' n_neighbors = 40, dotSNE = F, doClust = T,
#' ClusLouvRes = 0.00005, louvain_iter = 3,
#' KeepTopNgenes = 3000, minExprGeneDet = 0.1,
#' upperScale = 2, lowerScale = 2.5, assay = NULL)
Ser2Monocle_MakeNProcess <- function(SeurObj = NULL, retunMon = T, PCAnDim = 20,
                                     doUMAP = T, min_dist=.3, n_neighbors = 40,
                                     dotSNE = F,
                                     doClust = T, ClusLouvRes = 0.00005, louvain_iter = 3,
                                     KeepTopNgenes = 3000, 
                                     minExprGeneDet = 0.1, upperScale = 2, lowerScale = 2.5, assay = NULL){
  
  if (is.null(assay)){
    assay = DefaultAssay(SeurObj)
  }
  
  SeurObj.counts <- Seurat::GetAssayData(SeurObj, slot='counts', assay = assay)
  
  SeurObj.meta.data <- SeurObj@meta.data
  
  gene_ann <- data.frame(gene_short_name = row.names(SeurObj.counts), row.names = row.names(SeurObj.counts))
  
  pd <- new("AnnotatedDataFrame",data=SeurObj.meta.data)
  fd <- new("AnnotatedDataFrame",data=gene_ann)
  
  SeurObj_cds <- newCellDataSet(SeurObj.counts,
                                phenoData = pd, featureData =fd,
                                expressionFamily = negbinomial.size(),
                                lowerDetectionLimit=1)
  
  
  SeurObj_cds <- detectGenes(SeurObj_cds, min_expr = minExprGeneDet)
  
  
  # plot(density(SeurObj_cds$num_genes_expressed))
  
  
  
  SeurObj_cds <- estimateSizeFactors(SeurObj_cds)
  SeurObj_cds <- estimateDispersions(SeurObj_cds)
  
  disp_table = dispersionTable(SeurObj_cds)
  
  disp_table = disp_table %>% mutate(excess_disp =
                                       (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(dplyr::desc(excess_disp))
  
  top_subset_genes = as.character(head(disp_table, KeepTopNgenes)$gene_id)
  
  SeurObj_cds = setOrderingFilter(SeurObj_cds, top_subset_genes)
  
  
  pData(SeurObj_cds)$Total_mRNAs <- Matrix::colSums(exprs(SeurObj_cds))
  
  
  upper_bound <- 10^(mean(log10(pData(SeurObj_cds)$Total_mRNAs)) +
                       upperScale*sd(log10(pData(SeurObj_cds)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(SeurObj_cds)$Total_mRNAs)) -
                       lowerScale*sd(log10(pData(SeurObj_cds)$Total_mRNAs)))
  
  
  
  
  SeurObj_cds <- SeurObj_cds[,pData(SeurObj_cds)$Total_mRNAs > lower_bound &
                               pData(SeurObj_cds)$Total_mRNAs < upper_bound]
  
  SeurObj_cds <- detectGenes(SeurObj_cds, min_expr = minExprGeneDet)
  
  
  
  SeurObj_cds <- preprocessCDS(SeurObj_cds,
                               method = 'PCA',
                               norm_method = 'log',
                               num_dim = PCAnDim,
                               verbose = T)
  
  if(doUMAP) SeurObj_cds <- reduceDimension(SeurObj_cds, max_components = 2,
                                            reduction_method = 'UMAP',
                                            metric="correlation",
                                            min_dist = min_dist,
                                            n_neighbors = n_neighbors,
                                            verbose = T)
  
  # plot_pc_variance_explained(SeurObj_cds, return_all = F) # norm_method='log'
  
  if(dotSNE) SeurObj_cds <- reduceDimension(SeurObj_cds, max_components = 2, num_dim = 15,
                                            reduction_method = 'tSNE', verbose = T)
  
  if(doClust) SeurObj_cds <- clusterCells(SeurObj_cds,
                                          method = 'louvain',
                                          res = ClusLouvRes,
                                          louvain_iter = louvain_iter,
                                          verbose = T)
  
  
  
  
  
  return(SeurObj_cds)
  
  
}