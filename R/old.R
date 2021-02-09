
#' @title Run a variation (1b) of the primary seurat processing steps found in the CellMembrane package.
#'
#' @description This is the 1b entry point for processing (v1b) scRNAseq data with Seurat
#' @param seuratObj, A Seurat object.
#' @param saveFile If provided, the seuratObj will be saved here as it is processed, providing some ability to resume if there is a failure
#' @param doCellCycle If true, CellCycle genes will be regressed
#' @param doCellFilter If true, basic filtering will be performed using nCount_RNA, nFeature_RNA, and pMito
#' @param nCount_RNA.high If doCellFilter=T, cells with nCount_RNA above this value will be filtered
#' @param nCount_RNA.low If doCellFilter=T, cells with nCount_RNA below this value will be filtered
#' @param nFeature.high If doCellFilter=T, cells with nFeature above this value will be filtered
#' @param nFeature.low If doCellFilter=T, cells with nFeature below this value will be filtered
#' @param pMito.high If doCellFilter=T, cells with percent mito above this value will be filtered
#' @param pMito.low If doCellFilter=T, cells with percent mito  below this value will be filtered
#' @param forceReCalc If true, all steps will be repeated even if already marked as complete
#' @param variableGeneTable If provided, a table of variable genes will be written to this file
#' @param variableFeatureSelectionMethod The selection method to be passed to FindVariableFeatures()
#' @param useSCTransform If true, SCTransform will be used in place of the standard Seurat workflow (NormalizeData, ScaleData, FindVariableFeatures)
#' @param nVariableFeatures The number of variable features to find
#' @param dispersion.cutoff Passed directly to FindVariableFeatures
#' @param mean.cutoff Passed directly to FindVariableFeatures
#' @param spikeGenes If provided these will be appended to the set of VariableFeatures
#' @param excludedGenes If provided these will be removed from the set of VariableFeatures
#' @param printDefaultPlots If true, the default set of QC plots will be printed
#' @param npcs Number of PCs to use for RunPCA()
#' @param ccPcaResultFile If provided, the PCA results from cell cycle regression will be written here
#' @param Regress_basic if true p.mito and nCount_RNA are used to do a basic regression to remove the associated 'noise'. By default, this is F because by nature of the method, we are biasing the data.
#' @return A modified Seurat object.
# ProcessSeurat1b <- function(seuratObj, saveFile = NULL, doCellCycle = T, doCellFilter = F,
#                             nCount_RNA.high = 20000, nFeature.high = 3000, pMito.high = 0.15,
#                             nCount_RNA.low = 0.99, nFeature.low = 200, pMito.low = -Inf, forceReCalc = F,
#                             variableGeneTable = NULL, variableFeatureSelectionMethod = 'vst', 
#                             nVariableFeatures = 2000, printDefaultPlots = T,
#                             npcs = 50, ccPcaResultFile = NULL, useSCTransform = F, 
#                             mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), 
#                             spikeGenes = NULL, excludedGenes = NULL, verbose = FALSE, 
#                             Regress_basic = F){
#   
#   if (!forceReCalc && HasStepRun(seuratObj, 'ProcessSeurat1', forceReCalc = forceReCalc)) {
#     if (printDefaultPlots){
#       .PrintSeuratPlots(seuratObj, doCellCycle)
#     }
#     
#     return(seuratObj)
#   }
#   
#   if (doCellFilter & (forceReCalc | !HasStepRun(seuratObj, 'FilterCells', forceReCalc = forceReCalc))) {
#     seuratObj <- .DoCellFilter(seuratObj = seuratObj,
#                                nCount_RNA.high = nCount_RNA.high,
#                                nFeature.high = nFeature.high,
#                                pMito.high = pMito.high,
#                                nCount_RNA.low = nCount_RNA.low,
#                                nFeature.low = nFeature.low,
#                                pMito.low = pMito.low
#     )
#     
#     seuratObj <- MarkStepRun(seuratObj, 'FilterCells', saveFile)
#   }
#   
#   totalPMito = length(unique(seuratObj[['p.mito']]))
#   toRegress <- c("nCount_RNA")
#   if (totalPMito > 1) {
#     if(doCellCycle) toRegress <- c(toRegress, "p.mito")
#   }
#   
#   if(!Regress_basic) toRegress = NULL
#   
#   if (!useSCTransform) {
#     if (forceReCalc | !HasStepRun(seuratObj, 'NormalizeData', forceReCalc = forceReCalc)) {
#       print('Normalizing data:')
#       seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", verbose = F)
#       seuratObj <- MarkStepRun(seuratObj, 'NormalizeData', saveFile)
#     }
#     
#     if (forceReCalc | !HasStepRun(seuratObj, 'FindVariableFeatures', forceReCalc = forceReCalc)) {
#       print('Find variable features:')
#       seuratObj <- FindVariableFeatures(object = seuratObj, mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff , verbose = F, selection.method = variableFeatureSelectionMethod, nfeatures = nVariableFeatures)
#       seuratObj <- MarkStepRun(seuratObj, 'FindVariableFeatures', saveFile)
#     }
#     
#     if (forceReCalc | !HasStepRun(seuratObj, 'ScaleData', forceReCalc = forceReCalc)) {
#       print('Scale data:')
#       seuratObj <- ScaleData(object = seuratObj, features = rownames(x = seuratObj), vars.to.regress = toRegress, verbose = F)
#       seuratObj <- MarkStepRun(seuratObj, 'ScaleData', saveFile)
#     }
#     
#     if (doCellCycle & (forceReCalc | !HasStepRun(seuratObj, 'CellCycle', forceReCalc = forceReCalc))) {
#       seuratObj <- RemoveCellCycle(seuratObj, pcaResultFile = ccPcaResultFile, useSCTransform = F)
#       seuratObj <- MarkStepRun(seuratObj, 'CellCycle', saveFile)
#     }
#   } else {
#     print('Using SCTransform')
#     seuratObj <- SCTransform(seuratObj, vars.to.regress = toRegress, verbose = FALSE, return.only.var.genes = F)
#     
#     if (doCellCycle & (forceReCalc | !HasStepRun(seuratObj, 'CellCycle', forceReCalc = forceReCalc))) {
#       seuratObj <- RemoveCellCycle(seuratObj, pcaResultFile = ccPcaResultFile, useSCTransform = T)
#       seuratObj <- MarkStepRun(seuratObj, 'CellCycle', saveFile)
#     }
#   }
#   
#   if (!all(is.null(spikeGenes))){
#     missing <- spikeGenes[!(spikeGenes %in% rownames(seuratObj))]
#     if (length(missing) > 0) {
#       print(paste0("The following ", length(missing)," spikeGenes were found not in the seurat object: ", paste0(missing, collapse = ", ")))
#       spikeGenes <- spikeGenes[(spikeGenes %in% rownames(seuratObj))]
#     }
#     
#     VariableFeatures(seuratObj) <- unique(c(VariableFeatures(seuratObj), spikeGenes))
#   }
#   
#   if (!all(is.null(excludedGenes))) {
#     missing <- excludedGenes[!(excludedGenes %in% rownames(seuratObj))]
#     if (length(missing) > 0) {
#       print(paste0("The following ", length(missing)," excludedGenes were not found in the seurat object: ", paste0(missing, collapse = ", ")))
#       excludedGenes <- excludedGenes[(excludedGenes %in% rownames(seuratObj))]
#     }
#     
#     VariableFeatures(seuratObj) <- setdiff(VariableFeatures(seuratObj), excludedGenes)
#   }
#   
#   vg <- VariableFeatures(object = seuratObj)
#   
#   if (verbose) {
#     print("Final Variable Genes:")
#     print(paste(gtools::mixedsort(vg), collapse = ", "))
#   }
#   
#   if (forceReCalc | !HasStepRun(seuratObj, 'RunPCA', forceReCalc = forceReCalc)) {
#     seuratObj <- RunPCA(object = seuratObj, features = vg, verbose = F, npcs = npcs)
#     seuratObj <- MarkStepRun(seuratObj, 'RunPCA', saveFile)
#   }
#   
#   if (forceReCalc | !HasStepRun(seuratObj, 'ProjectDim', forceReCalc = forceReCalc)) {
#     seuratObj <- ProjectDim(object = seuratObj)
#     seuratObj <- MarkStepRun(seuratObj, 'ProjectDim', saveFile)
#   }
#   
#   #Verify data exists.  This appears to get reset, possibly by DimRedux steps
#   runJackStraw <- forceReCalc | !HasStepRun(seuratObj, 'JackStraw', forceReCalc = forceReCalc)
#   if (!useSCTransform) {
#     if (!runJackStraw && length(seuratObj@reductions$pca@jackstraw$empirical.p.values) == 0) {
#       warning('JackStraw marked as complete, but seurat object lacks data')
#       runJackStraw <- TRUE
#     }
#   }
#   
#   if (!useSCTransform && runJackStraw) {
#     seuratObj <- JackStraw(object = seuratObj, num.replicate = 100, verbose = F)
#     seuratObj <- ScoreJackStraw(object = seuratObj, dims = 1:20)
#     seuratObj <- MarkStepRun(seuratObj, 'JackStraw', saveFile)
#   }
#   
#   print(paste0('Total variable genes: ', length(vg)))
#   if (!is.null(variableGeneTable)){
#     write.table(sort(vg), file = variableGeneTable, sep = '\t', row.names = F, quote = F, col.names = F)
#   }
#   
#   if (printDefaultPlots){
#     .PrintSeuratPlots(seuratObj, doCellCycle)
#   }
#   
#   seuratObj <- MarkStepRun(seuratObj, 'ProcessSeurat1', saveFile = saveFile)
#   
#   print(seuratObj)
#   
#   return(seuratObj)
# }

