#' Convert GeoMxSet Object to Seurat
#'
#' @param x A GeoMxSet object.
#' @param normData Name of the normalized expression matrix (e.g., "exprs").
#' @param ident Optional identity column name in sample metadata to set Seurat identities.
#' @param coordinates Optional character vector of length 2 indicating x and y coordinate column names.
#' @param sequencingMetrics Optional vector of sequencing metric column names to retain in @misc.
#' @param QCMetrics Optional QC metrics column name to retain in @misc.
#'
#' @return A Seurat object.
convert_GeoMxSet_to_Seurat <- function(x, normData, ident = NULL, coordinates = NULL,
                                       sequencingMetrics = NULL, QCMetrics = "QCFlags") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed.")
  }
  
  meta <- as.data.frame(as.data.table(sData(x)[, !colnames(sData(x)) %in% c(sequencingMetrics, QCMetrics)]))
  
  projectName <- expinfo(experimentData(x))[["title"]]
  if (projectName == "") projectName <- "GeoMx"
  
  seuratConvert <- if (packageVersion("Seurat") < numeric_version("5.0.0")) {
    obj <- suppressWarnings(Seurat::CreateSeuratObject(
      counts = assayDataElement(x, normData),
      assay = "GeoMx",
      project = projectName
    ))
    obj <- suppressWarnings(Seurat::AddMetaData(obj, metadata = meta))
    obj@assays$GeoMx <- Seurat::AddMetaData(obj@assays$GeoMx, metadata = fData(x))
    obj
  } else {
    obj <- suppressWarnings(Seurat::CreateSeuratObject(
      counts = assayDataElement(x, "exprs"),
      assay = "GeoMx",
      project = projectName
    ))
    obj <- Seurat::SetAssayData(obj, layer = "data", new.data = assayDataElement(x, normData))
    obj <- suppressWarnings(Seurat::AddMetaData(obj, metadata = meta))
    obj@assays$GeoMx <- Seurat::AddMetaData(obj@assays$GeoMx, metadata = fData(x))
    obj
  }
  
  # Set Idents
  if (!is.null(ident)) {
    if (!ident %in% colnames(seuratConvert@meta.data)) {
      stop(paste0("ident \"", ident, "\" not found in GeoMxSet metadata"))
    }
    Seurat::Idents(seuratConvert) <- if (packageVersion("Seurat") < numeric_version("5.0.0")) {
      seuratConvert[[ident]]
    } else {
      as.factor(seuratConvert@meta.data[[ident]])
    }
  }
  
  # Attach misc metadata
  seuratConvert@misc <- otherInfo(experimentData(x))
  seuratConvert@misc[["sequencingMetrics"]] <- sData(x)[, colnames(sData(x)) %in% sequencingMetrics, drop = FALSE]
  seuratConvert@misc[["QCMetrics"]] <- sData(x)[, colnames(sData(x)) %in% QCMetrics, drop = FALSE]
  if (ncol(seuratConvert@misc[["QCMetrics"]]) == 0) {
    seuratConvert@misc[["QCMetrics"]] <- NULL
  }
  
  # Spatial coordinates
  if (!is.null(coordinates) && length(coordinates) == 2) {
    xcoord <- coordinates[1]
    ycoord <- coordinates[2]
    
    if (!all(c(xcoord, ycoord) %in% colnames(seuratConvert@meta.data))) {
      stop("One or both coordinate columns not found in metadata.")
    }
    
    coord.df <- data.frame(
      x = seuratConvert@meta.data[[xcoord]],
      y = seuratConvert@meta.data[[ycoord]]
    )
    colnames(coord.df) <- coordinates
    rownames(coord.df) <- rownames(seuratConvert@meta.data)
    
    # Remove coordinate columns from metadata
    seuratConvert@meta.data <- seuratConvert@meta.data[, !colnames(seuratConvert@meta.data) %in% coordinates]
    
    seuratConvert@images$image <- new(
      Class = "SlideSeq",
      assay = "GeoMx",
      key = "image_",
      coordinates = coord.df
    )
  }
  
  return(seuratConvert)
}