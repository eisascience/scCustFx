
#' Preprocess Seurat Object for Single-cell Data Analysis
#'
#' This function prepares a Seurat object for single-cell data analysis (SDA) by filtering features (genes)
#' and cells based on various criteria such as minimum feature count, minimum number of cells expressing a feature,
#' expression thresholds, inclusion and exclusion lists, and library size. It also checks and warns about 
#' the maximum number of cells the SDA can handle. The function utilizes the Seurat package for accessing
#' assay data and performs initial quality control checks and data normalization.
#'
#' @param seuratObj A Seurat object containing single-cell RNA sequencing data.
#' @param assayName The name of the assay to use, default is "RNA".
#' @param Layer The data layer to use, default is "counts".
#' @param minFeatureCount Minimum count for a feature to be included, default is 1.
#' @param minCellsExpressingFeature Minimum number of cells that must express a feature, default is 0.
#' @param perCellExpressionThreshold Threshold for expression level per cell, default is 0.
#' @param featureInclusionList A list of features to forcibly include, default is NULL.
#' @param featureExclusionList A list of features to exclude, default is NULL.
#' @param maxFeaturesDiscarded The maximum proportion (0-1) or absolute number of features to discard, default is 0.75.
#' @param minLibrarySize The minimum library size required for a cell to be included, default is 0.
#'
#' @return A list containing the filtered assay data (`SerObj.DGE`), the list of features used (`featuresToUse`),
#'         and the minimum library size considered (`minLibrarySize`).
#'
#' @export
SDA_ser_preproc <- function(seuratObj, 
                            assayName = "RNA", 
                            Layer = "counts",
                            minFeatureCount = 1,
                            minCellsExpressingFeature = 0,
                            perCellExpressionThreshold = 0,
                            featureInclusionList = NULL,
                            featureExclusionList = NULL,
                            maxFeaturesDiscarded = 0.75,
                            minLibrarySize = 0){
  

  
  SerObj.DGE <- Seurat::GetAssayData(seuratObj, 
                                     assay = assayName, 
                                     layer = Layer)
  
  n_cells <- ncol(SerObj.DGE)
  
  if (n_cells > 200000) {
    stop("SDA has shown to max handle ~200K cells ")
  } else if (n_cells > 150000) {
    warning("SDA has shown to max handle ~150K cells ")
  }
  
  print(paste0("Initial features: ", nrow(SerObj.DGE)))
  print(paste0("Initial cells: ", ncol(SerObj.DGE)))
  
  featuresToUse <- rownames(SerObj.DGE)
  
  
  
  
  if (!is.na(minFeatureCount) && minFeatureCount > 0) {
    
    
    numFeatures <- length(featuresToUse)
    
    hist(asinh(Matrix::rowSums(SerObj.DGE[featuresToUse, ])), breaks = 100, main = "lib size pre filter", xlab = "asinh(lib size)")
    
    
    featuresToUse <- featuresToUse[Matrix::rowSums(SerObj.DGE[featuresToUse, 
    ]) >= minFeatureCount]
    
    
    hist(asinh(Matrix::rowSums(SerObj.DGE[featuresToUse, ])), breaks = 100, main = "lib size post filter", xlab = "asinh(lib size)")
    
    
    print(paste0("After filtering to features with total counts of at least ", 
                 minFeatureCount, ": ", length(featuresToUse), " features remain (", 
                 scales::percent(length(featuresToUse)/numFeatures), 
                 " of input)"))
    
    
  }
  
  
  if (!is.na(minCellsExpressingFeature) && minCellsExpressingFeature > 
      0) {
    
    
    if (is.na(perCellExpressionThreshold)) {
      stop("Must provide perCellExpressionThreshold when minCellsExpressingFeature is above zero")
    }
    
    print("Filtering on minCellsExpressingFeature")
    
    if (minCellsExpressingFeature < 1) {
      minCellsExpressingFeatureRaw <- minCellsExpressingFeature
      minCellsExpressingFeature <- floor(minCellsExpressingFeatureRaw * 
                                           ncol(seuratObj))
      print(paste0("Interpreting minCellsExpressingFeature as a percentage of total cells (", 
                   ncol(seuratObj), "), converted from ", minCellsExpressingFeatureRaw, 
                   " to ", minCellsExpressingFeature))
    }
    
    numFeatures <- length(featuresToUse)
    numNonZeroCells <- Matrix::rowSums(SerObj.DGE >= perCellExpressionThreshold)
    
    featuresToUse <- names(numNonZeroCells[which(numNonZeroCells >= 
                                                   minCellsExpressingFeature)])
    print(paste0("After limiting to features with expression GTE ", 
                 perCellExpressionThreshold, " in at least ", minCellsExpressingFeature, 
                 " cells: ", length(featuresToUse), " features remain (", 
                 scales::percent(length(featuresToUse)/numFeatures), 
                 " of input)"))
    rm(numFeatures)
  }
  
  if (!all(is.null(featureInclusionList))) {
    featureInclusionList <- RIRA::ExpandGeneList(featureInclusionList)
    preExisting <- intersect(featuresToUse, featureInclusionList)
    print(paste0("Adding ", length(featureInclusionList), 
                 " features, of which ", length(preExisting), " are already present"))
    featuresToUse <- unique(c(featuresToUse, featureInclusionList))
    print(paste0("Total after: ", length(featuresToUse)))
  }
  
  if (!all(is.null(featureExclusionList))) {
    featureExclusionList <- RIRA::ExpandGeneList(featureExclusionList)
    preExisting <- intersect(featuresToUse, featureExclusionList)
    print(paste0("Excluding ", length(featureExclusionList), 
                 " features(s), of which ", length(preExisting), " are present"))
    featuresToUse <- unique(featuresToUse[!(featuresToUse %in% 
                                              featureExclusionList)])
    print(paste0("Total after: ", length(featuresToUse)))
  }
  
  if (length(featuresToUse) == 0) {
    stop("No features remain after filtering")
  }
  
  if (!is.null(maxFeaturesDiscarded)) {
    if (maxFeaturesDiscarded < 1) {
      maxFeaturesDiscarded <- maxFeaturesDiscarded * nrow(SerObj.DGE)
    }
    featsDiscarded <- nrow(SerObj.DGE) - length(featuresToUse)
    if (featsDiscarded > maxFeaturesDiscarded) {
      stop(paste0("The total number of features discarded, ", 
                  featsDiscarded, " exceeds the threshold of ", 
                  maxFeaturesDiscarded))
    }
  }
  
  df <- data.frame(x = Matrix::colSums(SerObj.DGE[featuresToUse, 
  ]))
  dens <- stats::density(df$x)
  mx <- dens$x[which.max(dens$y)]
  
  P1 <- ggplot(df, aes(x = x)) + scale_x_sqrt() + geom_density() + 
    ggtitle("Total features per cell") + labs(x = "Features/Cell", 
                                              y = "Density") + geom_vline(xintercept = mx, color = "red") + 
    ggtitle(paste0("Library Size: Peak = ", mx))
  
  if (!is.null(minLibrarySize)) {
    P1 <- P1 + geom_vline(xintercept = minLibrarySize, color = "red")
  }
  print(P1)
  print("starting dropsim normaliseDGE")
  
  return(list(SerObj.DGE = SerObj.DGE, featuresToUse = featuresToUse, minLibrarySize = minLibrarySize))
 
}



#' Prepare and Save Normalized Gene Expression Data for SDA
#'
#' This function prepares the environment for single-cell data analysis (SDA) by saving
#' normalized gene expression data to a specified output folder. It ensures that the
#' output directory exists, appending a trailing slash if necessary, and utilizes
#' SDAtools for data export. The function also prepares a results directory within
#' the output folder for further analysis.
#'
#' @param normedDGE A matrix or data frame containing normalized gene expression data,
#'   ready for analysis. This data is expected to be preprocessed and normalized appropriately.
#' @param outputFolder A character string specifying the path to the output folder where
#'   the normalized data and results will be saved. The function will create the directory
#'   if it does not exist, and will ensure it ends with a slash for consistency in path handling.
#'
#' @return A list containing two elements: `resultsDir` indicating the path to the results
#'   directory, and `rawDataFile` specifying the path to the saved normalized gene expression data file.
#'
#' @export
prepNormedDGE_SDA <- function(normedDGE, outputFolder){
  
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  if (!endsWith(outputFolder, "/")) {
    outputFolder <- paste0(outputFolder, "/")
  }
  
  
  if(class(normedDGE)[1]!="matrix") {
    print("normedDGE needs to be a non-sparse matrix")
    normedDGE= as.matrix(normedDGE)
  }
  print(paste0("Saving raw data to: ", outputFolder))
  
  
  SDAtools::export_data(normedDGE, path = outputFolder, name = "rawData")
  
  rawDataFile <- paste0(outputFolder, "rawData")
  resultsDir <- paste0(outputFolder, "results/")
  print(paste0("Saving results to: ", resultsDir))
  
  if (dir.exists(resultsDir)) {
    message(paste0("existing result folder: ", resultsDir))
    # unlink(resultsDir, recursive = TRUE)
  }
  
  # print("Running SDA")
  # if (!file.exists(path.sda)) {
  #   x <- unname(Sys.which(path.sda))
  #   if (x != "") {
  #     print(paste0("Found SDA under PATH: ", x))
  #     path.sda <- x
  #   }
  # }
  
  print(resultsDir)
  return(list(resultsDir = resultsDir, rawDataFile = rawDataFile))
  
}


#' Run SDA V1.1 in 2D i.e., matrix form
#'
#' Executes the SDA algorithm on preprocessed and normalized single-cell gene expression data,
#' allowing for dimensionality reduction and analysis. This function is designed to work with
#' the specified version of SDA, taking in several parameters to customize the analysis process,
#' including the number of components, iteration settings, and parallel processing options.
#'
#' @param path.sda A character string specifying the path to the SDA software executable.
#' @param resultsDir A character string specifying the directory where the SDA results will be saved.
#' @param rawDataFile A character string specifying the path to the preprocessed and normalized gene
#'   expression data file to be analyzed.
#' @param numComps Integer, the number of components to be used in the analysis, default is 30.
#' @param max_iter Integer, the maximum number of iterations for the SDA algorithm, default is 100.
#' @param save_freq Integer, frequency of saving interim results, default is 20 iterations.
#' @param set_seed Integer, the seed for random number generation to ensure reproducibility, default is 1234.
#' @param n_cells Integer, the number of cells in the dataset to be analyzed.
#' @param nThreads Integer, the number of threads for parallel processing, default is 6.
#'   If greater than 1, eigen decomposition will be done in parallel.
#' @param num_blocks Integer, the number of blocks to divide the dataset into for parallel processing, default is 6.
#'
#' @return Invisible NULL. The function is called for its side effects, including the execution
#'   of the SDA analysis and saving of results to the specified directory.
#' @export
runSDAv1_2D <- function(path.sda, resultsDir, rawDataFile, numComps = 30, max_iter = 100, save_freq = 20, 
                        set_seed = 1234, n_cells, nThreads = 6, num_blocks = 6){
  
  
  SDAtools::run_SDA(sda_location = path.sda, out = resultsDir, 
                    data = rawDataFile, num_comps = numComps, max_iter = max_iter, 
                    save_freq = save_freq, set_seed = set_seed, N = n_cells, 
                    eigen_parallel = (nThreads > 1), ignore_missing = FALSE, 
                    num_blocks = num_blocks, num_openmp_threads = nThreads)
  
  
  
}

#' Post-process Single-cell Data Analysis (SDA) Results
#'
#' This function loads the results from an SDA analysis, enriches them with additional metadata such
#' as cell barcodes and feature (gene) names, computes component statistics, and optionally performs
#' Gene Ontology (GO) enrichment analysis on the components. It is designed to work with results
#' generated by the SDAtools and provides a comprehensive overview of the analysis for further interpretation.
#'
#' @param resultsDir A character string specifying the directory where the SDA analysis results are stored.
#' @param outputFolder A character string specifying the path to the output folder where additional
#'   processed results and analyses will be saved.
#' @param normedDGE A matrix or data frame containing the normalized gene expression data used in the
#'   SDA analysis. Row names should correspond to features (genes), and column names should correspond to cell barcodes.
#' @param storeGoEnrichment A logical value indicating whether to perform and store Gene Ontology (GO)
#'   enrichment analysis on the analysis components. Default is `TRUE`.
#'
#' @return A list containing the enriched SDA results, including cell barcodes, feature names, component
#'   statistics, and optionally GO enrichment analysis results.
#' @export
SDApostprocessing <- function(resultsDir, 
                              outputFolder, 
                              normedDGE, 
                              orgDb = "org.Hs.eg.db", 
                              # mouse org.Mm.eg.db
                              # rhesus macaque org.Mmu.eg.db
                              storeGoEnrichment = T){
  
  results <- SDAtools::load_results(results_folder = resultsDir, 
                                    data_path = outputFolder)
  
  results$CellBarcodes <- colnames(normedDGE)
  results$Features <- rownames(normedDGE)
  results <-  SDA.AddCompStats(results)
  
  if (storeGoEnrichment) {
    results$goEnrichment <- SDA.GO_Enrichment(results, components = 1:nrow(results$loadings[[1]]), orgDb=orgDb)
  }
  
 
  return(results)
}




#' Diagnostic Checks on SDA Results
#'
#' Executes a series of diagnostic checks on Single-cell Data Analysis (SDA) results to
#' assess convergence, distribution of loadings and scores, and Posterior Inclusion Probability (PIP)
#' distributions. The function uses `SDAtools` for generating diagnostic outputs and handles
#' errors gracefully for each diagnostic step, ensuring that all possible checks are attempted.
#'
#' @details This function is designed to be run with SDA results available in the global environment
#'   or within its scope. It attempts to run convergence checks, loading and score distributions, and
#'   PIP distributions. If an error occurs during any of the checks, it catches the error, prints
#'   an error message specific to the failed check, and continues with the next checks.
#'
#' @return Invisible NULL. The function is called for its side effects, which include printing
#'   diagnostic information and error messages directly to the console.
#'
#' @export
SDAresChecks <- function(results){
  tryCatch({
    print(SDAtools::check_convergence(results))
    print(SDAtools::loading_distribution(results))
    print(SDAtools::scores_distribution(results))
  }, error = function(e) {
    print(paste0("Error generating SDA first plots"))
    print(conditionMessage(e))
  })
  tryCatch({
    print(SDAtools::PIP_distribution(results))
  }, error = function(e) {
    print(paste0("Error generating SDA PIP_distribution"))
    print(conditionMessage(e))
  })
  tryCatch({
    print(SDAtools::PIP_threshold_distribution(results))
  }, error = function(e) {
    print(paste0("Error generating SDA PIP_threshold_distribution"))
    print(conditionMessage(e))
  })
}


#' This function to add imputed scores from gene loading, V3 is one component at a time selecting topN genes
#'
#' @param SerObj A seurat Object, with assay RNA and slot data
#' @param sda_loadings matrix of sda loadings
#' @param keepComps if NULL all else numeric or char vec
#' @param sdaObjID id name to apped to comp names
#' @param plot to plot or not 
#' @return seurat obj with new SDA score comps in metadata
#' @export
ImputeSDA2SerV3 <- function(SerObj, sda_loadings, keepComps=NULL, sdaObjID="", plot=F, doAsinh = T,
                            TopNpos = 20, TopNneg=20, MetaExtraName=""){
  
  genes_overlap = intersect(rownames(SerObj), colnames(sda_loadings))
  
  
  if(is.null(keepComps)) keepComps = 1:nrow(sda_loadings)
  
  for (KC in keepComps){
    # KC  = keepComps[1]
    
    print(paste0("sda.", sdaObjID, ".V", KC, MetaExtraName))
    
    

    
    
    if(! (is.null(TopNpos) & is.null(TopNneg)) ) {
      print(paste0("Top N neg: ", TopNneg)); print(paste0("Top N pos: ", TopNpos))
      
      TopPos = names(sort(sda_loadings[KC, genes_overlap], decreasing = T))[1:TopNpos]
      TopNeg = names(sort(sda_loadings[KC, genes_overlap], decreasing = F))[1:TopNneg]
      
      genes_overlap = c(TopPos, TopNeg)
      print(genes_overlap)
    }
    
    sda_scores = Matrix::as.matrix(Matrix::t(sda_loadings[KC, genes_overlap] %*% SerObj@assays$RNA$data[genes_overlap, ]))
    colnames(sda_scores) = paste0("sda.", sdaObjID, ".V", KC, MetaExtraName)
    

    
    
    if(doAsinh) {
      SerObj = AddMetaData(SerObj, as.data.frame(asinh(sda_scores)))
    } else {
      SerObj = AddMetaData(SerObj, as.data.frame(sda_scores))
    }
    
    
    if(plot){
      print(Seurat::FeaturePlot(SerObj, features = colnames(sda_scores), order = T) & 
              ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
    }
    
  }
  return(SerObj)

}



#' This function to add imputed scores from gene loading
#'
#' @param SerObj A seurat Object, with assay RNA and slot data
#' @param sda_loadings matrix of sda loadings
#' @param keepComps if NULL all else numeric or char vec
#' @param sdaObjID id name to apped to comp names
#' @param plot to plot or not 
#' @param assay Default "RNA" can be "SCT", etc 
#' @return seurat obj with new SDA score comps in metadata
#' @export
ImputeSDA2SerV2 <- function(SerObj, sda_loadings, keepComps=NULL, sdaObjID="", 
                            plot=F, doAsinh = T, MakeReduc=T, assay = "RNA"){
  
  genes_overlap = intersect(rownames(SerObj), colnames(sda_loadings))
  
  if(is.null(keepComps)) keepComps = 1:nrow(sda_loadings)
  
  if(!assay %in% names(SerObj@assays)){
    print(paste0("Available Assays: ", paste0(names(SerObj@assays), collapse = ", ") )  )
    stop("Assay not in object ")
  }
  
  sda_scores = Matrix::as.matrix(Matrix::t(sda_loadings[keepComps, genes_overlap] %*% SerObj@assays[[assay]]$data[genes_overlap, ]))
  colnames(sda_scores) = paste0("sda.", sdaObjID, ".V", keepComps)
  
  if(doAsinh) {
    SerObj = AddMetaData(SerObj, as.data.frame(asinh(sda_scores)))
  } else {
    SerObj = AddMetaData(SerObj, as.data.frame(sda_scores))
  }
  
  if(plot){
    print(Seurat::FeaturePlot(SerObj, features = colnames(sda_scores), order = T) & 
            ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
  }
  
  
  
  
  
  CommonLoadingMat = t(as.matrix(sda_loadings)[keepComps, genes_overlap])
  colnames(CommonLoadingMat) = paste0("sda.", sdaObjID, ".V", keepComps)
  
  
  

  if(MakeReduc){
    SerObj = SDAScoreMeta2Reduction(SerObj = SerObj,
                                    sdaComps = colnames(CommonLoadingMat), 
                                    loadingMat = as.matrix(CommonLoadingMat),
                                    reduction.key =  paste0("sda", sdaObjID, "_"), 
                                    includeLoading = T,
                                    assayName = "RNA", reduction.name = paste0("sda_", sdaObjID) )
    
  }
 
  
  
  
  
  return(SerObj)
  
}



#' This function to add imputed scores from gene loading
#'
#' @param SerObj A seurat Object, with assay RNA and slot data
#' @param sda_loadings matrix of sda loadings
#' @param keepComps if NULL all else numeric or char vec
#' @param sdaObjID id name to apped to comp names
#' @param plot to plot or not 
#' @return seurat obj with new SDA score comps in metadata
#' @export
ImputeSDA2Ser <- function(SerObj, sda_loadings, keepComps=NULL, sdaObjID="", plot=F, doAsinh = T){
  
  genes_overlap = intersect(rownames(SerObj), colnames(sda_loadings))
  
  if(is.null(keepComps))keepComps = 1:nrow(sda_loadings)
  
  sda_scores = Matrix::as.matrix(Matrix::t(sda_loadings[keepComps, genes_overlap] %*% SerObj@assays$RNA$data[genes_overlap, ]))
  colnames(sda_scores) = paste0("sda.", sdaObjID, ".V", keepComps)

  if(doAsinh) {
    SerObj = AddMetaData(SerObj, as.data.frame(asinh(sda_scores)))
  } else {
    SerObj = AddMetaData(SerObj, as.data.frame(sda_scores))
  }
  
  if(plot){
    print(Seurat::FeaturePlot(SerObj, features = colnames(sda_scores), order = T) & 
            ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
  }
  
  
  return(SerObj)
  
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
#' @param SerObj A Seurat object containing SDA scores to use for the reduction.
#' @param loadingMat A matrix of loadings to use for the reduction. If not specified, the loadings will be calculated from the SDA scores.
#' @param sdaComps A vector of column names or indices from SerObj@meta.data that correspond to the SDA scores to use for the reduction.
#' @param reduction.key A prefix string to use for the key of the new dimension reduction object.
#' @param assayName The name of the assay in SerObj to use for the reduction.
#' @param reduction.name The name to use for the new reduction in the Seurat object.
#'
#' @return A modified Seurat object with a new dimension reduction object added.
#'
#'
#' @examples
#' ## Create a new SDA dimension reduction object from a Seurat object
#' my_seurat <- SDAScoreMeta2Reduction(SerObj = my_seurat, sdaComps = c("SDA_1", "SDA_2"))
#'
#' @export
SDAScoreMeta2Reduction <- function(SerObj, 
                                   loadingMat = NULL,
                                   sdaComps = NULL, reduction.key = 'SDA_' , 
                                   assayName = "RNA", reduction.name = "SDA", includeLoading = F){
  
  embeddings <- Matrix::as.matrix(ComboSerObj@meta.data[,sdaComps] )
  projected <- Matrix::as.matrix(loadingMat[,sdaComps])
  
  colnames(projected) <- paste0(reduction.key, 1:ncol(projected))
  colnames(embeddings) <-  colnames(projected)
  
  
  if(includeLoading) {
    sda.reduction <- Seurat::CreateDimReducObject(embeddings = embeddings, 
                                                  loadings = as.matrix(loadingMat),
                                                  projected = projected, key = reduction.key, assay = assayName)
  } else {
    sda.reduction <- Seurat::CreateDimReducObject(embeddings = embeddings, 
                                                  projected = projected, key = reduction.key, assay = assayName)
  }
  
  
  SerObj[[reduction.name]] <- sda.reduction
  return(SerObj)
}



#' Plot a heatmap of tabulation of SDA score vs meta
#'
#' THis fx makes a heatmap tabulating SDA scores within a seurat obj metadata
#'
#' @param SerObj A Seurat object containing SDA scores to use for the reduction.
#' @param componentNames a vector of meta features, the first is the tabulation meta the rest are sda scores... e.g. componentNames = c("ClusterNames_0.4", sdaComps[1:10])
#' @param direction direction to inspect
#' @param sdThr threshold to show sd on col labs default 0.2
#' @param doASINH boolean to do asinh transformation of residuals or not default F
#'
#' @return a complex heatmap
#' @export
Plot_SDAScoresPerFeature_ser <- function (SerObj, 
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
    SDAScores <- SerObj[[componentNames]]
    MetaDF <- SerObj@meta.data[, componentNames, drop = FALSE]
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
