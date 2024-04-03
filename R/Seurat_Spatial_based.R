

#' process Seurat Spatial object via SCT
#'
#' process Seurat Spatial object via SCT
#' 
#' @param SerObj a seurat obj 
#' @return A Processed Seurat Obj
#' @export
process_SCTbase = function(SerObj = NULL, ncells = 3000, assay = "Spatial", 
                           vst.flavor = "v2",
                           dims=1:30, verbose = T, 
                           do_tSNE = F,
                           clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2)){
  
  
  print("SCT transformation start")
  SerObj <- SCTransform(SerObj, assay = assay, 
                        ncells = ncells, # of cells used to build NB regression def is 5000
                        verbose = verbose,
                        vst.flavor = vst.flavor)
  print("Running PCA start")
  SerObj <- RunPCA(SerObj)
  if(do_tSNE) SerObj <- RunTSNE(SerObj, dims = dims, 
                                perplexity = CellMembrane:::.InferPerplexityFromSeuratObj(SerObj, perplexity = 30),  
                                check_duplicates = FALSE)
  print("Running UMAP start")
  SerObj <- RunUMAP(SerObj, dims = dims)
  print("Running Clustering start")
  SerObj <- FindNeighbors(SerObj, dims = dims)
  for (resolution in clusterResolutions) {
    SerObj <- FindClusters(object = SerObj, resolution = resolution, verbose = verbose)
    SerObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = SerObj)
  }
  return(SerObj)
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
#' @param SerObj A seurat object containing single cell transcriptomics data
#' @param ncells The number of cells used to build NB regression, default is 3000
#' @param assay The assay used to run the SCT transformation, default is "Spatial"
#' @param dims The dimensions to use for running PCA, t-SNE, UMAP, default is 1:30
#' @param verbose Logical, whether to print progress messages, default is TRUE
#' @param clusterResolutions A numeric vector of resolutions to use for finding clusters, default is c(0.2, 0.4, 0.6, 0.8, 1.2)
#' @return A seurat object with added dimensionality reductions, clustering, and cell type annotations.
#' @export
process_SCTbase = function(SerObj = NULL, ncells = 3000, assay = "Spatial",
                           dims=1:30, verbose = T, clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2)){
  print("SCT transformation start")
  SerObj <- SCTransform(SerObj, assay = assay, 
                        ncells = ncells, # of cells used to build NB regression def is 5000
                        verbose = verbose)
  print("Running PCA start")
  SerObj <- RunPCA(SerObj)
  SerObj <- RunTSNE(SerObj, dims = dims, perplexity = CellMembrane:::.InferPerplexityFromSeuratObj(SerObj, perplexity = 30),  check_duplicates = FALSE)
  print("Running UMAP start")
  SerObj <- RunUMAP(SerObj, dims = dims)
  print("Running Clustering start")
  SerObj <- FindNeighbors(SerObj, dims = dims)
  for (resolution in clusterResolutions) {
    SerObj <- FindClusters(object = SerObj, resolution = resolution, verbose = verbose)
    SerObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = SerObj)
  }
  return(SerObj)
}





#' load a h5ad SCE object and convert to seurat
#'
#' @param SCE_h5ad_file path to h5ad file
#' @param SER_rds_file path to save rds of ser obj
#' @param ReturnRawSer if T no normalization or filtering is done
#' @param verbose to print and plot or not T/F 
#'
#' @return a Seurat object for Spatial SlideSeq data
#'
#' @export
LoadSCE_ReturnSer <-function(SCE_h5ad_file,
                             ReturnRawSer =T, 
                             SER_rds_file,
                             verbose =T , 
                             nCount_Spatial_max=1000, nFeature_Spatial_min=20, 
                             nFeature_Spatial_max=700, percent.mt_max=10, 
                             LogNorm = F #if F SCT norm is done else logp1
){
  
  library(SingleCellExperiment)
  library(zellkonverter)
  
  
  sce <- readH5AD(SCE_h5ad_file, 
                  reader = "R",
                  verbose = T, 
                  version = "0.8.0")
  
  if(verbose) print(sce@metadata)
  
  
  if(verbose) print(rowData(sce)) # genes vs meta
  if(verbose) print(colData(sce))
  if(verbose) print(assayNames(sce))
  
  assayNames(sce) = "RNA"
  
  counts_data <- assays(sce)[["RNA"]]
  
  # if(verbose) print(dim(counts_data))
  # if(verbose) print(nrow(rowData(sce)))
  
  # dim(sce@assays@data$RNA)
  # spatialCoords(sce)
  # imgData(sce)
  # head(reducedDim(sce))
  
  # counts_data[1:10, 1:10]
  # sce@assays@data$RNA[1:10, 1:10]
  
  if(verbose) print(colnames(counts_data))
  if(verbose) print(rownames(counts_data))
  
  if( rowData(sce)$row_names[1] == "" ){
    print("gene names starts with blank position, odd, check processing")
    #this is odd
    counts_data = counts_data[-1,]
    rownames(counts_data) = rowData(sce)$row_names[-1]
  } else {
    rownames(counts_data) = rowData(sce)$row_names
  }
  
  
  # 
  # counts_data[-1,]
  
  # seurat_obj <- CreateSeuratObject(counts = counts_data, project = "STseq_Mouse01")
  # meta_data <- as.data.frame(colData(sce))
  # seurat_obj <- AddMetaData(seurat_obj, metadata = meta_data)
  # rownames(seurat_obj)
  # colnames(seurat_obj)
  
  
  SerObj = CreateSeuratObject(counts = counts_data, assay="Spatial")
  
  
  meta_data <- as.data.frame(colData(sce))
  SerObj <- AddMetaData(SerObj, metadata = meta_data)
  
  coord.df = data.frame(x=SerObj$x, y=SerObj$y, stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
  rownames(coord.df) = colnames(SerObj)
  
  # sfs <- scalefactors(spot = 138.656, fiducial = 223.9828, hires = 0.1139861, lowres = 0.03419583)
  SerObj@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    # scale.factors = sfs,
    coordinates = coord.df
  )
  rownames(SerObj)
  colnames(SerObj)
  
  if(verbose) print(SpatialDimPlot(SerObj, pt.size.factor = .2) + NoLegend())
  
  
  
  SerObj[["percent.mt"]] <- PercentageFeatureSet(SerObj, 
                                                 pattern = "^mt-")
  
  if(verbose) print( VlnPlot(
    SerObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
    pt.size = 0.1, ncol = 3) & 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()))
  
  
  # Jointly (rather than separately) consider the QC metrics when filtering
  if(verbose) print(FeatureScatter(
    SerObj, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend() + geom_hline(yintercept = c(10)) |
      FeatureScatter(
        SerObj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
      NoLegend() + geom_hline(yintercept = c(20, 700))  + geom_vline(xintercept = c(1000)))
  
  
  if(!ReturnRawSer){
    
    SerObj <- subset(SerObj, subset = nFeature_Spatial > nFeature_Spatial_min &
                       nFeature_Spatial < nFeature_Spatial_max & 
                       percent.mt < percent.mt_max )
    SerObj <- subset(SerObj, subset = nCount_Spatial < nCount_Spatial_max )
    
    if(verbose) print(VlnPlot(
      SerObj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
      pt.size = 0.1, ncol = 3) & 
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()))
    
    if(verbose) print(FeatureScatter(
      SerObj, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend() + geom_hline(yintercept = c(10)) |
        FeatureScatter(
          SerObj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
        NoLegend() + geom_hline(yintercept = c(20, 700))  + geom_vline(xintercept = c(1000)))
    
    
    if(verbose) print(SpatialFeaturePlot(
      SerObj, features = c("nFeature_Spatial"), pt.size.factor = .5) &
        theme(legend.position = "bottom") )
    
    
    if(verbose) print(SpatialFeaturePlot(
      SerObj, features = c("percent.mt"), pt.size.factor = .5) &
        theme(legend.position = "bottom")  )
    
    
    if(verbose) print(SpatialFeaturePlot(
      SerObj, features = c("nCount_Spatial"), pt.size.factor = .5) &
        theme(legend.position = "bottom")  )
    
    
    
    if (LogNorm){
      SerObj <- NormalizeData(SerObj, normalization.method = "LogNormalize", scale.factor = 10000)
      SerObj <- FindVariableFeatures(SerObj, selection.method = "vst", nfeatures = 2000)
      SerObj <- ScaleData(SerObj, features = VariableFeatures(SerObj))
      
    } else {
      SerObj <- SCTransform(SerObj, assay = "Spatial") #, clip.range = c(-10, 10)
      
    }
    
    top10 <- head(VariableFeatures(SerObj), 10)
    
    # plot variable features with and without labels
    
    if(verbose) print(LabelPoints(plot = VariableFeaturePlot(SerObj), points = top10, repel = TRUE))
    
    # all.genes <- rownames(SerObj)
    
    SerObj <- RunPCA(SerObj, npcs = 30, features = VariableFeatures(SerObj))
    if(verbose) print(ElbowPlot(SerObj, ndims = 30))
    
    SerObj <- RunUMAP(SerObj, dims = 1:30)
    SerObj <- FindNeighbors(SerObj, reduction = "pca", dims = 1:30)
    SerObj <- FindClusters(SerObj, resolution = 0.3)
    
    
    
    if(verbose) print(DimPlot(SerObj) + NoLegend())
    
    if(verbose) print(FeaturePlot(SerObj, "nCount_Spatial"))
    if(verbose) print(FeaturePlot(SerObj, "nFeature_Spatial"))
    
  }
  
  
  
  if(file.exists(SER_rds_file)){
    print(paste0(SER_rds_file, " file exists already"))
  } else {
    saveRDS(SerObj, SER_rds_file)
  }
  
  
  
  return(SerObj)
  
}



#' Filter Spatial data and plot before and after cell loss
#'
#' @param SerObj A SerObj
#' @param nCount_Spatial_max pnCount_Spatial_max
#' @param nFeature_Spatial_min pnCount_Spatial_max
#' @param nFeature_Spatial_max pnCount_Spatial_max
#' @param percent.mt_max pnCount_Spatial_max
#'
#' @return a Seurat object for Spatial SlideSeq data
#'
#' @export
FilterSpatial <- function(SerObj, nCount_Spatial_max=1000, nFeature_Spatial_min=20, 
                          nFeature_Spatial_max=700, percent.mt_max=10){
  
  gg1 = (FeatureScatter(
    SerObj, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend() + geom_hline(yintercept = c(percent.mt_max)) |
      FeatureScatter(
        SerObj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
      NoLegend() + geom_hline(yintercept = c(nFeature_Spatial_min, nFeature_Spatial_max))  + geom_vline(xintercept = c(nCount_Spatial_max)))
  
  print(gg1)
  
  print(ncol(SerObj))
  
  prevNcell = ncol(SerObj)
  
  SerObj <- subset(SerObj, subset = nFeature_Spatial > nFeature_Spatial_min &
                     nFeature_Spatial < nFeature_Spatial_max & 
                     percent.mt < percent.mt_max )
  
  SerObj <- subset(SerObj, subset = nCount_Spatial < nCount_Spatial_max )
  
  
  print(ncol(SerObj))
  
  afterNcell = ncol(SerObj)
  
  perLoss = round((prevNcell - afterNcell)/prevNcell * 100 , digits = 2)
  
  print(paste0("% loss: ", perLoss))
  
  gg1 = (FeatureScatter(
    SerObj, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend() + geom_hline(yintercept = c(percent.mt_max)) |
      FeatureScatter(
        SerObj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
      NoLegend() + geom_hline(yintercept = c(nFeature_Spatial_min, nFeature_Spatial_max))  + geom_vline(xintercept = c(nCount_Spatial_max)))
  
  print(gg1)
  
  return(SerObj)
  
}