
#' @import CellMembrane
#' @import Seurat

#' @title NormalizeAndScaleV2 a customized version from CellMembrane
#'
#' @description This is the entry point for processing raw scRNAseq data with Seurat & CellMembrane. 
#' @param SerObjLS, A named list of Seurat object.
#' @return Normalized and scaled seurat object
NormalizeAndScaleV2 <- function(SerObj, variableFeatureSelectionMethod = 'vst', 
                                nVariableFeatures = 2000, mean.cutoff = c(0.0125, 3), 
                                dispersion.cutoff = c(0.5, Inf), block.size = 1000,
                                normalization.method = "LogNormalize",
                                Regress_basic=F, Regress_pmito=T){
  
  SerObj <- NormalizeData(object = SerObj, normalization.method = normalization.method, verbose = F)
  
  print('Find variable features:')
  SerObj <- FindVariableFeatures(object = SerObj, mean.cutoff = mean.cutoff, 
                                    dispersion.cutoff = dispersion.cutoff , 
                                    verbose = F, selection.method = variableFeatureSelectionMethod, 
                                    nfeatures = nVariableFeatures)
  
  if ('p.mito' %in% names(SerObj@meta.data)) {
    totalPMito = length(unique(SerObj$p.mito))
  } else {
    totalPMito <- -1
  }
  
  toRegress <- c("nCount_RNA")
  if (totalPMito > 1) {
   if(Regress_pmito) toRegress <- c(toRegress, "p.mito")
  }
  
  if(!Regress_basic) toRegress = NULL
  
  print('Scale data:')
  SerObj <- ScaleData(object = SerObj, features = rownames(x = SerObj), vars.to.regress = toRegress, block.size = block.size, verbose = F)
  
  return(SerObj)
}




#' @title autoProcessRawSerObjs.basic: An automated and basic pipeline to process raw Seurat objects CellMembrane package.
#'
#' @description This is the entry point for processing raw scRNAseq data with Seurat & CellMembrane. 
#' @param SerObjLS, A named list of Seurat object.
#' @param verbose, set T/F for print to screen msgs
#' @param doDoubletFinder, removed doublets from the seurat objects
#' @param doRawCountFilter, filter cells based on defined parameters
#' @param doMergeObjects, filter cells based on defined parameters
#' @param doNomalizeAndScale, filter cells based on defined parameters
#' @param doRegressBasic, filter cells based on defined parameters
#' @param doRegresspmito, filter cells based on defined parameters
#' @param doRemoveCellCycle, filter cells based on defined parameters
#' @param doDimReduxTrio, filter cells based on defined parameters
#' @param doSingleR, filter cells based on defined parameters
#' @param nCount_RNA.high If doCellFilter=T, cells with nCount_RNA above this value will be filtered
#' @param nCount_RNA.low If doCellFilter=T, cells with nCount_RNA below this value will be filtered
#' @param nFeature.high If doCellFilter=T, cells with nFeature above this value will be filtered
#' @param nFeature.low If doCellFilter=T, cells with nFeature below this value will be filtered
#' @param pMito.high If doCellFilter=T, cells with percent mito above this value will be filtered
#' @param pMito.low If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableFeatureSelectionMethod If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param nVariableFeatures If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param mean.cutoff If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param dispersion.cutoff If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param block.size If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableGenesInclusion If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableGenesExclusion If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param npcs If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param minDimsToUse If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @return A processed seruat object
#' @export
autoProcessRawSerObjs.basic <- function(SerObjLS, verbose = T, 
                                        doDoubletFinder = T, dropDoublets = F, 
                                        featureData = NULL,
                                        doRawCountFilter = T,
                                        doMergeObjects = F, mergeObjectName = "MergedObjs",
                                        doNomalizeAndScale = T, doRegressBasic=T, doRegresspmito=T,
                                        doRemoveCellCycle = F,
                                        doDimReduxTrio = T,
                                        doSingleR = T,
                                        nCount_RNA.low = 0, nFeature.low = 200, pMito.low = 0, 
                                        nCount_RNA.high = 20000, nFeature.high = 5000, pMito.high = 0.15,
                                        variableFeatureSelectionMethod = 'vst', 
                                        nVariableFeatures = 2000, 
                                        mean.cutoff = c(0.0125, 3), 
                                        dispersion.cutoff = c(0.5, Inf), 
                                        block.size = 1000, 
                                        variableGenesInclusion = NULL, variableGenesExclusion = NULL, 
                                        npcs = 50, minDimsToUse = 15){
  
  if(!is.list(SerObjLS)) stop("SerObjLS needs to be a named list of seurat objects")
  
  #Doublet finder
  if(doDoubletFinder){
    if(verbose) print("starting doublet finder")
    newSerObjLS = list()
    for (datasetId in names(SerObjLS)) {
      if(verbose) print(datasetId)
      
      SerObj <- SerObjLS[[datasetId]]
      SerObj <- CellMembrane::FindDoublets(SerObj, dropDoublets = dropDoublets)
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  

  
  #Raw count filter
  if(doRawCountFilter){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::FilterRawCounts(SerObj, 
                                                 nCount_RNA.high = nCount_RNA.high, 
                                                 nCount_RNA.low = nCount_RNA.low, 
                                                 nFeature.high = nFeature.high, 
                                                 nFeature.low = nFeature.low, 
                                                 pMito.high = pMito.high, pMito.low = pMito.low) 
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #merge of ser objects to a single ser obj
  if(doMergeObjects){
    newSerObjLS = list()
    #the original SerObjLS is not needed as the next steps will be don on the alamgam object
    newSerObjLS[[mergeObjectName]] <- CellMembrane::MergeSerObjs(SerObjLS, projectName = mergeObjectName)
    SerObjLS = newSerObjLS
  }
  
  
  #Normalization and scaling steps 
  if(doNomalizeAndScale){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      gc()
      
      SerObj <- NormalizeAndScaleV2(SerObj,
                                       variableFeatureSelectionMethod = variableFeatureSelectionMethod, 
                                       nVariableFeatures = nVariableFeatures, 
                                       mean.cutoff = mean.cutoff, 
                                       dispersion.cutoff = dispersion.cutoff, 
                                       block.size = block.size,
                                       Regress_basic=doRegressBasic, Regress_pmito=doRegresspmito)
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #cell cycle removal
  if(doRemoveCellCycle){
    
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RemoveCellCycle(SerObj, block.size = block.size, min.genes = 10)
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #PCA, tSNE, UMAP
  if(doDimReduxTrio){
    
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RunPcaSteps(SerObj, npcs = npcs)
      SerObj <- CellMembrane::FindClustersAndDimRedux(SerObj, minDimsToUse = minDimsToUse)
      
      newSerObjLS[[datasetId]] <- SerObj
      
    }
    SerObjLS = newSerObjLS
  }
  
  
  if(doSingleR){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RunSingleR(SerObj)
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  
  return(SerObjLS)
  
}





#' @title autoProcessRawSerObjs.prime: An automated and basic pipeline to process raw Seurat objects obtained form prime-seq with CellMembrane package.
#'
#' @description This is the entry point for processing raw scRNAseq data with Seurat & CellMembrane. 
#' @param SerObjLS, A named list of Seurat object.
#' @param verbose, set T/F for print to screen msgs
#' @param doDoubletFinder, removed doublets from the seurat objects
#' @param doAppendCellHashing, filter cells based on defined parameters
#' @param doAppendCITEseq, filter cells based on defined parameters
#' @param doRawCountFilter, filter cells based on defined parameters
#' @param doMergeObjects, filter cells based on defined parameters
#' @param doNomalizeAndScale, filter cells based on defined parameters
#' @param doRegressBasic, filter cells based on defined parameters
#' @param doRegresspmito, filter cells based on defined parameters
#' @param doRemoveCellCycle, filter cells based on defined parameters
#' @param doDimReduxTrio, filter cells based on defined parameters
#' @param doSingleR, filter cells based on defined parameters
#' @param nCount_RNA.high If doCellFilter=T, cells with nCount_RNA above this value will be filtered
#' @param nCount_RNA.low If doCellFilter=T, cells with nCount_RNA below this value will be filtered
#' @param nFeature.high If doCellFilter=T, cells with nFeature above this value will be filtered
#' @param nFeature.low If doCellFilter=T, cells with nFeature below this value will be filtered
#' @param pMito.high If doCellFilter=T, cells with percent mito above this value will be filtered
#' @param pMito.low If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableFeatureSelectionMethod If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param nVariableFeatures If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param mean.cutoff If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param dispersion.cutoff If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param block.size If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableGenesInclusion If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param variableGenesExclusion If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param npcs If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param minDimsToUse If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @return A processed seruat object
#' @export
autoProcessRawSerObjs.prime <- function(SerObjLS, verbose = T, 
                                        doDoubletFinder = T, dropDoublets = F,
                                        doAppendCellHashing = T, featureData = NULL,
                                        doAppendCITEseq = T,
                                        doRawCountFilter = T,
                                        doMergeObjects = F, mergeObjectName = "MergedObjs",
                                        doNomalizeAndScale = T, doRegressBasic=T, doRegresspmito=T,
                                        doRemoveCellCycle = F,
                                        doDimReduxTrio = T,
                                        doSingleR = T,
                                        nCount_RNA.low = 0, nFeature.low = 200, pMito.low = 0, 
                                        nCount_RNA.high = 20000, nFeature.high = 5000, pMito.high = 0.15,
                                        variableFeatureSelectionMethod = 'vst', 
                                        nVariableFeatures = 2000, 
                                        mean.cutoff = c(0.0125, 3), 
                                        dispersion.cutoff = c(0.5, Inf), 
                                        block.size = 1000, 
                                        variableGenesInclusion = NULL, variableGenesExclusion = NULL, 
                                        npcs = 50, minDimsToUse = 15){
  
  if(!is.list(SerObjLS)) stop("SerObjLS needs to be a named list of seurat objects")
  
  #Doublet finder
  if(doDoubletFinder){
    if(verbose) print("starting doublet finder")
    newSerObjLS = list()
    for (datasetId in names(SerObjLS)) {
      if(verbose) print(datasetId)
      
      SerObj <- SerObjLS[[datasetId]]
      SerObj <- CellMembrane::FindDoublets(SerObj, dropDoublets = dropDoublets)
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #Cell hashing (unique to .prime)
  if(is.null(featureData)) {
    if(doAppendCellHashing) print("featureData is NULL skipping AppendCellHashing")
    doAppendCellHashing = F
  }
  
  if(doAppendCellHashing){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      if (!(datasetId %in% names(featureData))) {
        stop(paste0('No hashing information found for datasetId: ', datasetId))
      }
      
      callFile <- featureData[[datasetId]]
      if (!is.null(callFile)) {
        SerObj <- cellhashR::AppendCellHashing(SerObj, barcodeCallFile = callFile, barcodePrefix = datasetId)
      } else {
        # Add empty columns to keep objects consistent
        SerObj$HTO <- c(NA)
        SerObj$consensuscall.global <- c(NA)
      }
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #doAppendCITEseq
  
  #doAppendTCR
  
  #Raw count filter
  if(doRawCountFilter){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::FilterRawCounts(SerObj, 
                                                 nCount_RNA.high = nCount_RNA.high, 
                                                 nCount_RNA.low = nCount_RNA.low, 
                                                 nFeature.high = nFeature.high, 
                                                 nFeature.low = nFeature.low, 
                                                 pMito.high = pMito.high, pMito.low = pMito.low) 
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #merge of ser objects to a single ser obj
  if(doMergeObjects){
    newSerObjLS = list()
    #the original SerObjLS is not needed as the next steps will be don on the alamgam object
    newSerObjLS[[mergeObjectName]] <- CellMembrane::MergeSerObjs(SerObjLS, projectName = mergeObjectName)
    SerObjLS = newSerObjLS
  }
  

  #Normalization and scaling steps 
  if(doNomalizeAndScale){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      gc()
      
      SerObj <- NormalizeAndScaleV2(SerObj,
                                                   variableFeatureSelectionMethod = variableFeatureSelectionMethod, 
                                                   nVariableFeatures = nVariableFeatures, 
                                                   mean.cutoff = mean.cutoff, 
                                                   dispersion.cutoff = dispersion.cutoff, 
                                                   block.size = block.size,
                                                   Regress_basic=doRegressBasic, Regress_pmito=doRegresspmito)
        
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #cell cycle removal
  if(doRemoveCellCycle){
    
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RemoveCellCycle(SerObj, block.size = block.size, min.genes = 10)
      
      newSerObjLS[[datasetId]] <- SerObj
  
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  #PCA, tSNE, UMAP
  if(doDimReduxTrio){

    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RunPcaSteps(SerObj, npcs = npcs)
      SerObj <- CellMembrane::FindClustersAndDimRedux(SerObj, minDimsToUse = minDimsToUse)
    
      newSerObjLS[[datasetId]] <- SerObj
  
    }
    SerObjLS = newSerObjLS
  }
  
  
  if(doSingleR){
    for (datasetId in names(SerObjLS)) {
      SerObj <- SerObjLS[[datasetId]]
      
      SerObj <- CellMembrane::RunSingleR(SerObj)
      
      newSerObjLS[[datasetId]] <- SerObj
      
      gc()
    }
    SerObjLS = newSerObjLS
  }
  
  
 return(SerObjLS)
  
}


