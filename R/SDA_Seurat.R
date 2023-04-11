
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