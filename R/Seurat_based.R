#' @title SubsetSerByCoordinates
#'
#' @description Remove unwanted HTOS from Seurat objects
#' @param so A Seurat Object 
#' @param reduction_method an input of  = c("tSNE", "UMAP", "PCA") default tSNE
#' @param x_range an input of  two numbers defining the x range default c(-1, 1)
#' @param y_range aan input of  two numbers defining the y range default c(-1, 1)
#' @return A subset Seurat object.
#' @export
SubsetSerByCoordinates <- function(so, reduction_method = "tSNE", x_range = c(-1, 1), y_range = c(-1, 1)){
  if(reduction_method == "tSNE"){
    dim_1 <- "TSNE1"
    dim_2 <- "TSNE2"
  }else if(reduction_method == "UMAP"){
    dim_1 <- "UMAP1"
    dim_2 <- "UMAP2"
  }else if(reduction_method == "PCA"){
    dim_1 <- "PCA1"
    dim_2 <- "PCA2"
  }
  subsetted_obj <- SubsetData(so, subset = (so@reductions[[reduction_method]]@cell.embeddings[,dim_1] > x_range[1] & seurat_obj@reductions[[reduction_method]]@cell.embeddings[,dim_1] < x_range[2]) & (seurat_obj@reductions[[reduction_method]]@cell.embeddings[,dim_2] > y_range[1] & seurat_obj@reductions[[reduction_method]]@cell.embeddings[,dim_2] < y_range[2]))
  return(subsetted_obj)
}



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



#' @title compressSerHD5
#'
#' @description compressSerHD5 will save space on your drives.
#' @param so A Seurat Object 
#' @param load.path the path to a seurat .rds file
#' @param save.path the path to the desired hd5 seurat file. Default is NULL which then load.path is used to save the new object next to it.
#' @param overwrite  default F overwirtes the new hd5 if exists
#' @param updateSerObj  default F if T runs UpdateSeuratObject()
#' @param returnSerObj  default F if T returns the Seurat object to be used in a pipeline
#' @return Saves HD5 Seurat object to path given
#' @export
compressSerHD5 <- function(so = NULL, load.path = NULL, overwrite = F,
                           save.path = NULL, updateSerObj = F, returnSerObj = F){
  if(is.null(so) & is.null(load.path)) stop("give seurat object so or load.path to one")
  
  if((!is.null(so)) & (!is.null(load.path))) warning("both so and load.path given, load.path is ignored")
  
  
  if(is.null(so)){
    print("loading in RDS")
    so = readRDS(load.path)
    print("load of RDS from drives complete")
  }
  
  if(updateSerObj) so = SeuratDisk::UpdateSeuratObject(so)
  
  if(is.null(save.path)) save.path = gsub(".rds", "hd5Ser", load.path)
  print("saving file path:")
  print(save.path)
  SeuratDisk::SaveH5Seurat(so, filename=save.path,  overwrite = overwrite)
  print("save complete")
  if(returnSerObj) return(so)
  
}