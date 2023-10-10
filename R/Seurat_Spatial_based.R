

#' process Seurat Spatial object via SCT
#'
#' process Seurat Spatial object via SCT
#' 
#' @param serObj a seurat obj 
#' @return A Processed Seurat Obj
#' @export
process_SCTbase = function(serObj = NULL, ncells = 3000, assay = "Spatial", 
                           vst.flavor = "v2",
                           dims=1:30, verbose = T, 
                           do_tSNE = F,
                           clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2)){
  
  
  print("SCT transformation start")
  serObj <- SCTransform(serObj, assay = assay, 
                        ncells = ncells, # of cells used to build NB regression def is 5000
                        verbose = verbose,
                        vst.flavor = vst.flavor)
  print("Running PCA start")
  serObj <- RunPCA(serObj)
  if(do_tSNE) serObj <- RunTSNE(serObj, dims = dims, 
                                perplexity = CellMembrane:::.InferPerplexityFromSeuratObj(serObj, perplexity = 30),  
                                check_duplicates = FALSE)
  print("Running UMAP start")
  serObj <- RunUMAP(serObj, dims = dims)
  print("Running Clustering start")
  serObj <- FindNeighbors(serObj, dims = dims)
  for (resolution in clusterResolutions) {
    serObj <- FindClusters(object = serObj, resolution = resolution, verbose = verbose)
    serObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = serObj)
  }
  return(serObj)
}



#' Create a Seurat object from Spatial SlideSeq data
#'
#' @param dge.path path to the expression data
#' @param pos.path path to the position data
#' @param celltype.path path to the cell type data
#' @param celltype.header header for the cell type column in celltype.path
#' @param radius radius for filtering cells
#' @param do.plot logical indicating whether to plot the filtered cells
#' @param doFilter logical indicating whether to filter cells
#' @param log_nCount.bool logical indicating whether to log transform nCount
#'
#' @return a Seurat object for Spatial SlideSeq data
#'
#' @export
SpatialSer.SlideSeq.local = function(dge.path=NULL,
                                     pos.path = NULL, celltype.path = NULL, celltype.header = "max_cell_type",
                                     radius = 2450, do.plot = T, doFilter = T, log_nCount.bool = T){
  
  if(is.null(dge.path) | is.null(pos.path)) stop("paths are null")
  if(!file.exists(dge.path) | !file.exists(pos.path)) stop("files dont exist")
  
  
  expr <- data.table::fread(
    input = dge.path,
    sep = ",",
    data.table = FALSE
  )
  if (any(c(0, 1) %in% expr[, 1])) expr = expr[,-1]
  rownames(x = expr) <- expr[, 1] #toupper(x = expr[, 1])
  # colnames(x = expr) <- toupper(x = colnames(x = expr))
  expr <- t(x = expr[, -1])
  
  SpatialSer <- CreateSeuratObject(
    counts = expr,
    project = 'SlideSeq',
    assay = 'Spatial'
  )
  
  positions <- read.csv(
    file = pos.path
  )
  if (any(c(0, 1) %in% positions[, 1])) positions = positions[,-1]
  
  rownames(x = positions) <- positions$barcode
  positions <- positions[, 2:3]
  
  #  SpatialSer[['image']] <- new(
  #    Class = 'SlideSeq',
  #    assay = "Spatial",
  #    coordinates = positions
  # )
  
  SpatialSer@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = positions
  )
  
  if(!is.null(celltype.path)){
    
    celltypes <- read.csv(
      file = celltype.path
    )
    celltypesV = celltypes[,celltype.header]
    names(celltypesV) = celltypes$barcode
    
    SpatialSer = AddMetaData(SpatialSer, metadata = celltypesV, col.name = "celltype.input")
  }
  
  if(doFilter) SpatialSer <- FilterSlideSeq(object = SpatialSer, radius = radius, do.plot = T)
  
  if(log_nCount.bool) SpatialSer$log_nCount_Spatial <- log(SpatialSer$nCount_Spatial)
  
  
  return(SpatialSer)
}

#' Process single cell transcriptomics data with SCT transformation
#'
#' This function takes in a seurat object and applies various transformations and clustering algorithms. The end goal is to identify cell types and subpopulations within the data.
#'
#' @param serObj A seurat object containing single cell transcriptomics data
#' @param ncells The number of cells used to build NB regression, default is 3000
#' @param assay The assay used to run the SCT transformation, default is "Spatial"
#' @param dims The dimensions to use for running PCA, t-SNE, UMAP, default is 1:30
#' @param verbose Logical, whether to print progress messages, default is TRUE
#' @param clusterResolutions A numeric vector of resolutions to use for finding clusters, default is c(0.2, 0.4, 0.6, 0.8, 1.2)
#' @return A seurat object with added dimensionality reductions, clustering, and cell type annotations.
#' @export
process_SCTbase = function(serObj = NULL, ncells = 3000, assay = "Spatial",
                           dims=1:30, verbose = T, clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2)){
  print("SCT transformation start")
  serObj <- SCTransform(serObj, assay = assay, 
                        ncells = ncells, # of cells used to build NB regression def is 5000
                        verbose = verbose)
  print("Running PCA start")
  serObj <- RunPCA(serObj)
  serObj <- RunTSNE(serObj, dims = dims, perplexity = CellMembrane:::.InferPerplexityFromSeuratObj(serObj, perplexity = 30),  check_duplicates = FALSE)
  print("Running UMAP start")
  serObj <- RunUMAP(serObj, dims = dims)
  print("Running Clustering start")
  serObj <- FindNeighbors(serObj, dims = dims)
  for (resolution in clusterResolutions) {
    serObj <- FindClusters(object = serObj, resolution = resolution, verbose = verbose)
    serObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = serObj)
  }
  return(serObj)
}