
#' This function to add imputed scores from gene loading
#'
#' @param serObj A seurat Object, with assay RNA and slot data
#' @param sda_loadings matrix of sda loadings
#' @param keepComps if NULL all else numeric or char vec
#' @param sdaObjID id name to apped to comp names
#' @param plot to plot or not 
#' @return seurat obj with new SDA score comps in metadata
#' @export
ImputeSDA2Ser <- function(serObj, sda_loadings, keepComps=NULL, sdaObjID="", plot=F, doAsinh = T){
  
  genes_overlap = intersect(rownames(serObj), colnames(sda_loadings))
  
  if(is.null(keepComps))keepComps = 1:nrow(sda_loadings)
  
  sda_scores = Matrix::as.matrix(Matrix::t(sda_loadings[keepComps, genes_overlap] %*% serObj@assays$RNA@data[genes_overlap, ]))
  colnames(sda_scores) = paste0("sda.", sdaObjID, ".V", keepComps)

  if(doAsinh) {
    serObj = AddMetaData(serObj, as.data.frame(asinh(sda_scores)))
  } else {
    serObj = AddMetaData(serObj, as.data.frame(sda_scores))
  }
  
  if(plot){
    print(Seurat::FeaturePlot(serObj, features = colnames(sda_scores), order = T) & 
            ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
  }
  
  
  return(serObj)
  
}

#' Combine loading matrices across multiple SDA runs
#'
#' This function takes a list of SDA runs and returns a matrix of loadings that are common across all runs.
#'
#' @param sdaRuns A list of SDA runs, each containing a loadings matrix.
#'
#' @return A matrix of loadings that are common across all SDA runs.
#' 
#' @examples
#' ## Create a matrix of common loadings across multiple SDA runs
#' my_loadings <- CommonLoadingSDAMat(sdaRuns = my_sda_runs)
#'
#' @export
CommonLoadingSDAMat = function(sdaRuns){
  
  CommonGenes = lapply(sdaRuns, function(x){
    colnames(x$loadings[[1]])
  })
  CommonGenes = Reduce(intersect, CommonGenes)
  
  CommonLoadingMat = lapply(sdaRuns, function(x){
    t(x$loadings[[1]][,CommonGenes])
  })
  CommonLoadingMat = do.call(cbind, CommonLoadingMat)
  
  return(CommonLoadingMat)
}


#' Create a new dimension reduction object using SDA scores
#'
#' This function creates a new dimension reduction object using SDA scores from a Seurat object, and adds it to the Seurat object as a new assay.
#'
#' @param seuratObj A Seurat object containing SDA scores to use for the reduction.
#' @param loadingMat A matrix of loadings to use for the reduction. If not specified, the loadings will be calculated from the SDA scores.
#' @param sdaComps A vector of column names or indices from seuratObj@meta.data that correspond to the SDA scores to use for the reduction.
#' @param reduction.key A prefix string to use for the key of the new dimension reduction object.
#' @param assayName The name of the assay in seuratObj to use for the reduction.
#' @param reduction.name The name to use for the new reduction in the Seurat object.
#'
#' @return A modified Seurat object with a new dimension reduction object added.
#'
#'
#' @examples
#' ## Create a new SDA dimension reduction object from a Seurat object
#' my_seurat <- SDAScoreMeta2Reduction(seuratObj = my_seurat, sdaComps = c("SDA_1", "SDA_2"))
#'
#' @export
SDAScoreMeta2Reduction <- function(seuratObj, 
                                   loadingMat = NULL,
                                   sdaComps = NULL, reduction.key = 'SDA_' , 
                                   assayName = "RNA", reduction.name = "SDA"){
  
  embeddings <- Matrix::as.matrix(seuratObj@meta.data[,sdaComps] )
  projected <- Matrix::as.matrix(loadingMat)
  
  colnames(projected) <- paste0(reduction.key, 1:ncol(projected))
  
  sda.reduction <- Seurat::CreateDimReducObject(embeddings = embeddings, 
                                                projected = projected, key = reduction.key, assay = assayName)
  
  seuratObj[[reduction.name]] <- sda.reduction
  return(seuratObj)
}



#' Plot a heatmap of tabulation of SDA score vs meta
#'
#' THis fx makes a heatmap tabulating SDA scores within a seurat obj metadata
#'
#' @param seuratObj A Seurat object containing SDA scores to use for the reduction.
#' @param componentNames a vector of meta features, the first is the tabulation meta the rest are sda scores... e.g. componentNames = c("ClusterNames_0.4", sdaComps[1:10])
#' @param direction direction to inspect
#' @param sdThr threshold to show sd on col labs default 0.2
#' @param doASINH boolean to do asinh transformation of residuals or not default F
#'
#' @return a complex heatmap
#' @export
Plot_SDAScoresPerFeature_ser <- function (seuratObj, 
                                          componentNames, 
                                          direction = "Both", doASINH = F,
                                          sdThr = 0.2) 
{
  # componentNames = c("ClusterNames_0.4", sdaComps[1:10])
  
  
  if (direction == "Both") {
    directions <- c("Pos", "Neg")
  }
  else {
    directions <- direction
  }
  for (direction in directions) {
    # direction = "Neg"
    SDAScores <- seuratObj[[componentNames]]
    MetaDF <- seuratObj@meta.data[, componentNames, drop = FALSE]
    SDAScores <- SDAScores[rownames(MetaDF), , drop = FALSE]
    lvl <- levels(factor(MetaDF[, 1, drop = TRUE]))
    
    if (length(lvl) < 2) {
      print(paste0("Only one value for field, skipping: ", 
                   componentNames))
      return()
    }
    
    CompsDF <- as.data.frame(lapply(lvl, function(CondX) {
      
      cellIDs = rownames(MetaDF)[which(MetaDF[, 1, 
                                              drop = TRUE] == CondX)]
      
      apply(SDAScores[cellIDs, 2:ncol(SDAScores), drop = FALSE], 2, 
            function(x) {
              if (direction == "Neg") {
                round(sum(x < 0)/nrow(SDAScores) * 100, 2)
              }
              else if (direction == "Pos") {
                round(sum(x > 0)/nrow(SDAScores) * 100, 2)
              }
              else {
                stop(paste0("Direction must be Pos or Neg: ", 
                            direction))
              }
            })
    }))
    
    colnames(CompsDF) <- levels(factor(MetaDF[, 1]))
    CompsDF <- as.data.frame(CompsDF[naturalsort::naturalsort(rownames(CompsDF)), 
    ])
    
    # CompsDF$SDA <- factor(rownames(CompsDF), levels = rownames(CompsDF))
    
    # ChiT <- stats::chisq.test(CompsDF[, 1:(ncol(CompsDF) - 
    #                                          1)])
    
    ChiT <- stats::chisq.test(CompsDF[])
    
    ChiTres <- ChiT$residuals
    ChiTres[which(is.na(ChiTres))] <- 0
    ChiResSD <- round(apply(ChiTres, 1, sd), 2)
    ChiResSD[which(is.na(ChiResSD))] <- 0
    ChiResSD[ChiResSD < sdThr ] <- ""
    
    if(doASINH) plotM = asinh(t(ChiTres)) else plotM = asinh(t(ChiTres))
    
    HM <- ComplexHeatmap::Heatmap(plotM, name = direction, 
                                  column_title = paste0("Direction: ", direction), 
                                  cluster_columns = TRUE, cluster_rows = TRUE, 
                                  col = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                                                                                                  name = "RdBu"))))(10), column_labels = paste0(rownames(CompsDF), 
                                                                                                                                                " sd_", ChiResSD))
    return(HM)
  }
}
