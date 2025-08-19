

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


ConvertFOVToSpatialSeurat <- function(object, assay.name = "Spatial", project.name = "ConvertedSpatial") {
  
  # **1. Validate Input**
  if (is.null(object)) stop("Seurat object is NULL.")
  #TODO this limits to only 1 image which for now is all we have but some seurat objects could have more.. 
  if (!"FOV" %in% class(object@images[[1]])) stop("object@images[[1]] is not an FOV object.")
  
  # **2. Extract Expression Data from Old Seurat Object**
  expr <- GetAssayData(object, layer = "counts")  # Extract raw counts
  if (is.null(expr) || ncol(expr) == 0) stop("Expression data is missing in object.")
  
  # **3. Create a New Seurat Object**
  object_new <- CreateSeuratObject(
    counts = expr,
    project = project.name,
    assay = assay.name
  )
  
  # **4. Extract Spatial Coordinates from the FOV Object**
  positions <- data.frame(
    barcode = object@images[[1]]@boundaries$centroids@cells,  # Extract barcodes
    x = object@images[[1]]@boundaries$centroids@coords[, "x"],
    y = object@images[[1]]@boundaries$centroids@coords[, "y"]
  )
  
  # Ensure barcode alignment
  rownames(positions) <- positions$barcode
  positions <- positions[, c("x", "y")]  # Keep only coordinates
  
  # **5. Attach Spatial Image Data**
  object_new@images$spatial_image <- new(
    Class = "SlideSeq",
    assay = assay.name,
    key = "image_",
    coordinates = positions
  )
  
  # **6. Transfer Metadata from Old Seurat Object**
  object_new@meta.data <- object@meta.data
  
  # **7. Ensure Proper Rowname Alignment**
  rownames(object_new@meta.data) <- colnames(object_new)  # Align metadata rownames with cell names
  
  # **8. Return the New Spatial Seurat Object**
  return(object_new)
}



#' SpatialPlotGridAligned
#'
#' Arrange multiple samples from a Seurat object in a grid layout using their spatial coordinates,
#' preserving within‐sample morphology and annotating each panel with a metadata label.
#'
#' @param seurat_obj A \code{\link[Seurat]{Seurat}} object containing spatial metadata columns
#'   \code{x_slide_mm}, \code{y_slide_mm}, \code{SampleID}, and the grouping and labeling columns.
#' @param group.by Character. Name of the metadata column to use for coloring points (e.g.,
#'   \code{"seurat_clusters"}).
#' @param spacing Numeric. Vertical and horizontal spacing (in same units as slide coordinates)
#'   between panels in the grid.
#' @param ncol Integer. Number of columns in the output grid.
#' @param label.by Character. Name of the metadata column to use for panel titles
#'   (e.g., \code{"Tissue"} or \code{"SampleID"}).
#' @param col_vector Named character vector. Custom colors for each level of \code{group.by}.
#'   Names must match the factor levels; if \code{NULL}, a \code{viridis} discrete palette is used.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object showing all samples arranged in a grid,
#'   colored by the specified grouping and labeled by the specified metadata.
#'
#' @examples
#' \dontrun{
#'   # Default palette, 2 columns
#'   SpatialPlotGridAligned(SerObjIntg, group.by = "seurat_clusters",
#'                          spacing = 20, ncol = 2, label.by = "Tissue")
#'
#'   # Custom color vector
#'   my_colors <- c("0" = "#440154", "1" = "#414487", "2" = "#2A788E")
#'   SpatialPlotGridAligned(SerObjIntg, col_vector = my_colors)
#' }
#'
#' @import ggplot2 dplyr tibble
#' @export
SpatialPlotGridAligned <- function(seurat_obj,
                                   group.by = "seurat_clusters",
                                   spacing = 2,
                                   ncol = 2,
                                   label.by = "Tissue",
                                   col_vector = NULL) {
  meta <- seurat_obj@meta.data
  coords <- meta %>%
    tibble::rownames_to_column(var = "cell") %>%
    dplyr::select(cell, SampleID, x_slide_mm, y_slide_mm,
                  group = !!sym(group.by), label = !!sym(label.by))
  
  bounds <- coords %>%
    group_by(SampleID, label) %>%
    summarise(
      min_x = min(x_slide_mm),
      min_y = min(y_slide_mm),
      width = max(x_slide_mm) - min(x_slide_mm),
      height = max(y_slide_mm) - min(y_slide_mm),
      .groups = "drop"
    )
  
  bounds <- bounds %>%
    mutate(index = row_number() - 1,
           col = index %% ncol,
           row = index %/% ncol) %>%
    group_by(col) %>%
    mutate(offset_x = cumsum(lag(width + spacing, default = 0))) %>%
    ungroup() %>%
    group_by(row) %>%
    mutate(offset_y = cumsum(lag(height + spacing, default = 0))) %>%
    ungroup()
  
  coords_adj <- coords %>%
    left_join(bounds, by = c("SampleID", "label")) %>%
    mutate(
      x_adj = (x_slide_mm - min_x) + offset_x,
      y_adj = -(y_slide_mm - min_y) - offset_y
    )
  
  titles <- coords_adj %>%
    group_by(SampleID, label) %>%
    summarise(
      x_center = mean(range(x_adj)),
      y_top = min(y_adj) - spacing / 2,
      .groups = "drop"
    )
  
  p <- ggplot(coords_adj, aes(x = x_adj, y = y_adj, color = group)) +
    geom_point(size = 0.5, alpha = 0.6) +
    geom_text(data = titles, aes(x = x_center, y = y_top, label = label),
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    coord_fixed(expand = TRUE) +
    theme_void() +
    theme(legend.position = "right") +
    ggtitle(paste("Spatial Grid View:", group.by))
  
  # Apply color mapping
  if (!is.null(col_vector)) {
    p <- p + scale_color_manual(values = col_vector)
  } else {
    p <- p + scale_color_viridis_d(option = "C")
  }
  
  return(p)
}
