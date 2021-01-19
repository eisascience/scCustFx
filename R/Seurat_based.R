
#' @title ReduceSerObj_HTO
#'
#' @description Remove unwanted HTOS from Seurat objects
#' @param so A Seurat Object 
#' @param removeHTOs Vector of names of HTOs to remove, if NULL, removeHTOs = c("Discordant", "Doublet", "ND")
#' @return A subset Seurat object.
#' @export
ReduceSerObj_HTO <- function(so, removeHTOs = NULL){
  if(is.null(removeHTOs)){
    print("removeHTOs is null, therefore, deep cleaning")
    removeHTOs = c("Discordant", "Doublet", "ND")
  }
  if (sum(removeHTOs %in% names(table(so$HTO)))==0){
    stop("removeHTOs defined not in seurat object ")
  } else {
    keepHTOs <- setdiff(names(table(so$HTO)), removeHTOs)
  }
  
  print("Keept HTOs: ")
  print(keepHTOs)
  return(subset(so, HTO %in% keepHTOs))
}

#' @title DownSample_SerObj_perLevel
#'
#' @description downsample per feature levels
#' @param so A Seurat Object 
#' @param featureSplit a feature of seurat object to split. "BarcodePrefix" # this is a good way to sample per seq run if possible
#' @param totalCell  0 # is the default way to downsample by min common cells
#' @return A subseted downsampled Seurat object.
#' @export
DownSample_SerObj_perLevel <- function(so, featureSplit=NULL, totalCell = 300){

  availFeats = colnames(so@meta.data)
  
  if(!is.null(featureSplit)){
    
    if(length(featureSplit)>1) {
      warning("length of featureSplit > 1 only 1st is used")
      featureSplit = featureSplit[1]
      
      if(!(featureSplit %in% availFeats)) stop("featureSplit not found in meta.data")
      
    }
    
  } else stop("featureSplit is NULL")
  
  if(totalCell == 300) warning("totalCell = 300 is the default")
  
  CellCountMat = table(so@meta.data[,featureSplit]) %>% unlist() %>% as.matrix()
  
  if(totalCell == 0) {
    totalCell = min(CellCountMat[,1])
    print(paste0("total cells set by min to:", totalCell))
  }
  
  if(!("barcode" %in% availFeats)) {
    warning("barcode not found in meta.data")
    so$barcode = paste("cell_", 1:nrow(so@meta.data))
  }
  
  
  splitLevs = levels(factor(so@meta.data[,featureSplit]))
  
  barcodeLS = lapply(splitLevs, function(xSL){
    availBarcodes = so@meta.data[which(so@meta.data[,featureSplit] == xSL),]$barcode
    if(length(availBarcodes)<totalCell) {
      warning(paste0(xSL, " had less cells than totalCell"))
      availBarcodes
    }else {
      sample(availBarcodes, totalCell, replace = F)
    }
    
    
  })
  
  return(subset(so, barcode %in% unlist(barcodeLS)))
  
}
