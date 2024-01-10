#' Plot Gene Expression in Seurat Object
#'
#' This function generates a violin plot or pie chart for the expression of a specified gene
#' in a Seurat object, based on a given threshold.
#'
#' @param SerObj A Seurat object containing single-cell RNA-seq data.
#' @param gene The name of the gene for which expression is to be plotted.
#' @param threshold The expression threshold for selecting cells (default is 0.5).
#' @param groupBy The variable by which the data is grouped (default is "RNA_snn_res.0.4").
#' @param doPie Logical, if TRUE, a pie chart is generated instead of a violin plot (default is FALSE).
#'
#'
#'
#' @export
plotGeneExpression <- function(serObj, gene, threshold = 0.5, groupBy = "RNA_snn_res.0.4", doPie = FALSE, col_vector = col_vector) {
  
  # Get the expression values for the specified gene
  geneExpression <- GetAssayData(object = serObj, assay = "RNA", layer = "counts")[gene, ]
  
  # Get cell ids based on the threshold and gene expression
  cellIds <- rownames(geneExpression)[geneExpression > threshold]
  
  # Subset the Seurat object
  subsetObj <- subset(serObj, cells = cellIds)
  
  # Create a table and print it
  cat("Table of", groupBy, "distribution:\n")
  phenoTable <- table(subsetObj[[groupBy]])
  print(phenoTable)
  
  # Create a pie chart if doPie is TRUE
  if (doPie) {
    gg_pie_table(phenoTable)
  } else {
    # Create a violin plot
    VlnPlot(subsetObj, 
            group.by = groupBy, 
            features = gene, 
            cols = col_vector, 
            raster = FALSE)
  }
}



#' Most_Common_Clones_Per
#' Function to compute the most common clones for a TCR feature in a Seurat object, grouped by a second feature
#'
#' @param SerObj A Seurat object containing the data.
#' @param tcr_feature A metadata feature e.g. TRA_J
#' @param grouping_feature A metadata feature to group by e.g. Tissue
#' @param top_n a numerical value, default 10
#' 
#' @export
Most_Common_Clones_Per <- function(SerObj, tcr_feature, grouping_feature, top_n = 10) {
  # Extract TCR feature and grouping feature from the Seurat object
  tcr_data <- SerObj[[tcr_feature]]
  tcr_data$grouping_data <- as.character(SerObj[[grouping_feature]][,1])
  
  # Initialize an empty list to store top clones and their counts for each group
  top_clones_list <- list()
  
  # Get unique levels of the grouping feature
  grouping_levels <- unique(tcr_data$grouping_data)
  
  # Iterate through each level of the grouping feature
  for (level in grouping_levels) {
    # Subset data for the current level of the grouping feature
    subset_data <- tcr_data[tcr_data$grouping_data == level, 1]
    
    # Split clones separated by comma and remove NA values
    clones <- unlist(strsplit(as.character(subset_data), ","))[!is.na(subset_data)]
    
    # Count the frequency of each clone
    clone_counts <- table(clones)
    
    # Sort clones by frequency in descending order
    sorted_clones <- sort(clone_counts, decreasing = TRUE)
    
    # Get the top N most common clones for the current level
    top_clones <- names(sorted_clones)[1:min(length(sorted_clones), top_n)]
    
    # Remove any leading or trailing spaces in clone names and remove NA
    top_clones <- trimws(gsub("\"| NA", "", top_clones))
    top_clones <- top_clones[top_clones != ""]
    
    # Store top clones and their counts in a data frame
    top_clones_df <- data.frame(Clone = top_clones, Count = sorted_clones[top_clones])
    
    # Add top clones and their counts for the current level to the list
    top_clones_list[[as.character(level)]] <- top_clones_df
  }
  
  return(top_clones_list)
}


#' Most_Common_Clones_2_Matrix
#' A function that takes in a list, from Most_Common_Clones_Per()
#'
#' @param top_clones_list A list of data frames with Clone and Count.Freq
#' 
#' @export
Most_Common_Clones_2_Matrix <- function(top_clones_list) {
  # Extract clone names and their counts for each group
  clone_names <- unique(unlist(lapply(top_clones_list, function(df) df$Clone)))
  grouping_levels <- names(top_clones_list)
  
  # Initialize an empty matrix with rows as clones and columns as grouping levels
  freq_matrix <- matrix(0, nrow = length(clone_names), ncol = length(grouping_levels),
                        dimnames = list(clone_names, grouping_levels))
  
  # Fill the matrix with clone frequencies
  for (i in 1:length(grouping_levels)) {
    group <- as.character(grouping_levels[i])
    clone_df <- top_clones_list[[group]]
    clone_indices <- match(clone_df$Clone, clone_names)
    freq_matrix[clone_indices, i] <- clone_df$Count.Freq
  }
  
  return(freq_matrix)
}


#' Most_Common_Clones
#' Function to compute the most common clones for a TCR feature in a Seurat object
#'
#' @param SerObj A Seurat object containing the data.
#' @param tcr_feature A metadata feature e.g. TRA_J
#' @param top_n a numerical value, default 10
#'
#' @export
Most_Common_Clones <- function(SerObj, tcr_feature, top_n = 10) {
  # Extract TCR feature from the Seurat object
  tcr_data <- SerObj[[tcr_feature]]
  
  # Split clones separated by comma and remove NA values
  clones <- unlist(strsplit(as.character(tcr_data), ","))[!is.na(tcr_data)]
  
  # Count the frequency of each clone
  clone_counts <- table(clones)
  
  # Sort clones by frequency in descending order
  sorted_clones <- sort(clone_counts, decreasing = TRUE)
  
  # Get the top N most common clones
  top_clones <- names(sorted_clones)[1:min(length(sorted_clones), top_n)]
  
  # Remove any leading or trailing spaces in clone names and remove NA
  top_clones <- trimws(gsub("\"| NA", "", top_clones))
  top_clones = top_clones[top_clones!=""]
  return(top_clones)
}

#' Find_FOXP3
#' Identify FOXP3+ cells
#'
#' @param SerObj A Seurat object containing the data.
#' @return A Seurat object with the IsFOXP3 metadata column added
#'
#' @export
Find_FOXP3 <- function(SerObj, plot = FALSE) {
  
  SerObj$IsFOXP3 = "FOXP3-"
  SerObj$IsFOXP3[WhichCells(SerObj, expression = FOXP3 > 0)] = 'FOXP3+'
  

    if (plot) {
      print(DimPlot(SerObj, group.by = 'IsFOXP3', raster = TRUE, raster.dpi = c(800, 800), cols = c("gray", "maroon")))
      print(pheatmap::pheatmap(asinh(chisq.test(table(SerObj$IsFOXP3, SerObj$Population))$res)))
    }
    
    SerObj$IsMAITclassic <- "Not TRAV1-2"
    SerObj$IsMAITclassic[SerObj$Has_TRAV1_2 == "TRAV1-2"] <- "TRAV1-2"
    SerObj$IsMAITclassic[SerObj$Has_TRAV1_2 == "TRAV1-2" & SerObj$TRA_J == "TRAJ33"] <- "TRAV1-2 TRAJ33"
    

  
  return(SerObj)
}


#' Find_TCR_MAITS
#' Identify TRAV1-2 cells and Classical MAITS as TRAV1-2 TRAJ33
#'
#' @param SerObj A Seurat object containing the data.
#' @return A Seurat object with the IsMAITclassic and Has_TRAV1_2 metadata columns added
#'
#' @export
Find_TCR_MAITS <- function(SerObj, plot = FALSE) {
  
  TRAV1_2cells <- grepl("TRAV1-2", SerObj$TRA_V)
  if(!is.null(TRAV1_2cells)){
    SerObj$Has_TRAV1_2 <- "Not TRAV1-2"
    SerObj$Has_TRAV1_2[TRAV1_2cells] <- "TRAV1-2"
    
    if (plot) {
      print(DimPlot(SerObj, group.by = 'Has_TRAV1_2', raster = TRUE, raster.dpi = c(800, 800), cols = c("gray", "maroon")))
      print(pheatmap::pheatmap(asinh(chisq.test(table(SerObj$Has_TRAV1_2, SerObj$Population))$res)))
    }
    
    SerObj$IsMAITclassic <- "Not TRAV1-2"
    SerObj$IsMAITclassic[SerObj$Has_TRAV1_2 == "TRAV1-2"] <- "TRAV1-2"
    SerObj$IsMAITclassic[SerObj$Has_TRAV1_2 == "TRAV1-2" & SerObj$TRA_J == "TRAJ33"] <- "TRAV1-2 TRAJ33"
    
    
  } else {
    print("TRAV1-2 not found in SerObj$TRA_V")
  }
  
  return(SerObj)
}



#' Bool.CellsGT
#'
#' This function checks whether the values in a specific gene's assay data of a Seurat object are greater than a given threshold.
#'
#' @param SerObj A Seurat object.
#' @param thresh A numeric threshold value.
#' @param layer The layer name of the Seurat object to retrieve the data from.
#' @param assay The assay name to retrieve the data from.
#' @param gene The name of the gene to subset the data.
#'
#' @export
Bool.CellsGT <- function(SerObj, thresh, layer = "data", assay = "RNA", gene=NULL){
  tempTab <- GetAssayData(SerObj, layer = layer, assay = assay)[gene, ]
  print(summary(tempTab))
  tempTab>thresh
}


#' CalculateModuleScoreAvg
#'
#' This function calculates the module scores for each cell in a Seurat object based on the average expression of a given gene list.
#'
#' @param SerObj A Seurat object.
#' @param genes.list A character vector specifying the genes to calculate the module scores.
#' @param genes.pool A character vector specifying the genes to use for calculating the average expression.
#' @param n.bin An integer specifying the number of bins for partitioning the average expression values.
#' @param ctrl.size An integer specifying the size of the control group for calculating module scores.
#' @param assay The assay name to retrieve the data from. If NULL, the default assay of the Seurat object is used.
#' @param seed.use An integer specifying the seed value for random number generation.
#'
#' @export
CalculateModuleScoreAvg <- function (SerObj, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
                                     ctrl.size = 100, assay = NULL, seed.use = 1) 
{
  set.seed(seed = seed.use)
  genes.old <- genes.list
  if (is.null(x = genes.list)) {
    stop("Missing input gene list")
  }
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(SerObj)
  }
  genes.list <- intersect(x = genes.list, y = rownames(SerObj))
  if (length(genes.list) == 0) {
    print(paste0("No matching genes found in the seurat object from the gene list, attempting to match case."))
    genes.list <- Seurat::CaseMatch(genes.old, match = rownames(SerObj))
  }
  if (length(genes.list) == 0) {
    print(paste0("No matching genes found in the seurat object from the gene list, aborting"))
    return(NA)
  }
  if (is.null(x = genes.pool)) {
    genes.pool = rownames(SerObj)
  }
  data.avg <- Matrix::rowMeans(Seurat::GetAssayData(SerObj, 
                                                    assay = assay)[genes.pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
  names(x = data.cut) <- names(x = data.avg)
  genes.use <- genes.list
  ctrl.use <- character()
  for (j in 1:length(x = genes.use)) {
    ctrl.use <- c(ctrl.use, names(x = sample(x = data.cut[which(x = data.cut == 
                                                                  data.cut[genes.use[j]])], size = ctrl.size, replace = FALSE)))
  }
  ctrl.use <- unique(ctrl.use)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = Seurat::GetAssayData(SerObj, assay = assay)))
  data <- Seurat::GetAssayData(SerObj, assay = assay)
  ctrl.scores <- Matrix::colMeans(x = data[ctrl.use, , drop = F])
  genes.scores <- Matrix::colMeans(x = data[genes.list, , drop = FALSE])
  
  
  return(data.frame(CellBarcode = colnames(x = Seurat::GetAssayData(SerObj, 
                                                                    assay = assay)), Score = (genes.scores - ctrl.scores)))
}




#' Rank a feature and plot a swarm grouped by select categorical metadata
#' 
#' @param SerObj A Seurat object to perform gating on.
#' @param RankFeat A feature to rank, like SDA comp
#' @param group.by A categorical meta feature
#' @param color.by A categorical meta feature default NULL to color the point
#' @param title A title
#' 
#' @return ggplot viz of ranking per group
#' 
#' @export
Rank_Swarm <- function(SerObj, RankFeat, group.by, title = "Cell SerObjrted", color.by = NULL){
  
  #TODO if PC or UMAP RANKFeat given to pull that seperately 
  
  library(ggbeeswarm)
  
  SerObj[[paste0(RankFeat, "_rank")]] <- rank(SerObj[[RankFeat]][,1]) 
  
  if(is.null(color.by)){
    MetaDF = SerObj@meta.data[, c(paste0(RankFeat, "_rank"), group.by )]
    
    colnames(MetaDF) = c("pseudoT", "Grp")
    
   gg1 = ggplot(MetaDF, aes(x = pseudoT, y = Grp, 
                       colour = Grp)) 
      
  } else {
    
    MetaDF = SerObj@meta.data[, c(paste0(RankFeat, "_rank"), c(group.by, color.by) )]
    
    colnames(MetaDF) = c("pseudoT", "Grp", "Col")
    
    gg1 = ggplot(MetaDF, aes(x = pseudoT, y = Grp, 
                             colour = Col)) 
    
  }
  
  gg1 + geom_quasirandom(groupOnX = FALSE) +
    ggthemes::scale_color_tableau() + theme_classic() +
    xlab(RankFeat) + ylab(group.by) +
    ggtitle(title)
  
}


#' Perform dual-feature gating on a Seurat object
#' 
#' @param SerObj A Seurat object to perform gating on.
#' @param feat1 A character string representing the name of the first feature to gate on. Defaults to NULL.
#' @param feat2 A character string representing the name of the second feature to gate on. Defaults to NULL.
#' @param thr1 A numeric threshold value for the first feature. Defaults to 0.
#' @param thr2 A numeric threshold value for the second feature. Defaults to 0.
#' @param dir1 A character string representing the direction of the threshold for the first feature, either "pos" or "neg". Defaults to "pos".
#' @param dir2 A character string representing the direction of the threshold for the second feature, either "pos" or "neg". Defaults to "pos".
#' @param plotDimplot A logical value indicating whether to plot the result using Seurat's DimPlot function. Defaults to TRUE.
#' @param returnSerObj A logical value indicating whether to return the modified Seurat object or just the gated cell names. Defaults to TRUE.
#' 
#' @return A Seurat object with the compGate metadata column added if \code{returnSerObj} is TRUE, otherwise a character vector of gated cell names.
#' 
#' @export
DualFeatGate <- function(SerObj, 
                         feat1 = NULL, 
                         feat2 = NULL, 
                         thr1 = 0, thr2 = 0, 
                         dir1 = "pos", dir2 = "pos",
                         cols = c("maroon", "gray", "blue", "forestgreen"),
                         plotDimplot = T, returnSerObj=T, reduction = "umap",
                         label = F, label.size = 6, repel = T, 
                         raster = F) {
  
  score1 <- SerObj[[feat1]]
  score2 <- SerObj[[feat2]]
  
  if(dir1=="pos") {
    cells1 <- rownames(score1)[score1[,1] >= thr1]
  } else {
    cells1 <- rownames(score1)[score1[,1] <= thr1]
  }
  
  if(dir2=="pos") {
    cells2 <- rownames(score2)[score2[,1] >= thr2]
  } else {
    cells2 <- rownames(score2)[score2[,1] <= thr2]
  }
  
  SerObj[["compGate"]] <- "neither"
  SerObj[["compGate"]][cells1,] <- feat1
  SerObj[["compGate"]][cells2,] <- feat2
  SerObj[["compGate"]][intersect(cells1, cells2),] <- "both"
  
  if(plotDimplot){
    p <- DimPlot(SerObj, 
                 group.by = "compGate",
                 cols = cols, 
                 reduction = reduction,
                 label = label, label.size = label.size, repel = repel, 
                 raster = raster) + 
      NoLegend() +
      facet_wrap(~compGate)
    
    print(p)
  }
  
  if(returnSerObj) {
    return(SerObj)
  } else {
    return(  p  )
  }
  
}



#' Perform dual-feature gating on a Seurat object
#' 
#' @param SerObj A Seurat object to perform gating on.
#' @param feat1 A character string representing the name of the first feature to gate on. Defaults to NULL.
#' @param feat2 A character string representing the name of the second feature to gate on. Defaults to NULL.
#' @param thr1 A numeric threshold value for the first feature. Defaults to 0.
#' @param thr2 A numeric threshold value for the second feature. Defaults to 0.
#' @param dir1 A character string representing the direction of the threshold for the first feature, either "pos" or "neg". Defaults to "pos".
#' @param dir2 A character string representing the direction of the threshold for the second feature, either "pos" or "neg". Defaults to "pos".
#' @param plotDimplot A logical value indicating whether to plot the result using Seurat's DimPlot function. Defaults to TRUE.
#' @param returnSerObj A logical value indicating whether to return the modified Seurat object or just the gated cell names. Defaults to TRUE.
#' 
#' @return A Seurat object with the compGate metadata column added if \code{returnSerObj} is TRUE, otherwise a character vector of gated cell names.
#' 
#' @export
SingleFeatGate <- function(SerObj, 
                         feat1 = NULL, 
                         thr1 = 0, 
                         dir1 = "pos",
                         cols = c("maroon", "gray"),
                         plotDimplot = T, returnSerObj=T) {
  
  score1 <- SerObj[[feat1]]

  if(dir1=="pos") {
    cells1 <- rownames(score1)[score1[,1] >= thr1]
  } else {
    cells1 <- rownames(score1)[score1[,1] <= thr1]
  }
  
  
  SerObj[["compGate"]] <- "not"
  SerObj[["compGate"]][cells1,] <- feat1

  if(plotDimplot){
    p <- DimPlot(SerObj, 
                 group.by = "compGate",
                 cols = cols, 
                 label = F, label.size = 6, repel = T, 
                 raster = F) + 
      NoLegend() +
      facet_wrap(~compGate)
    
    print(p)
  }
  
  if(returnSerObj) {
    return(SerObj)
  } else {
    return(  SerObj[["compGate"]]  )
  }
  
}




#' Get feature windows
#'
#' This function generates feature windows for a given feature and groups in a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param feat A character vector representing the feature name.
#' @param group.by A character vector representing the group to use for grouping the data.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{winDF}{A table with the counts of cells falling into each window, for each group.}
#'   \item{vlines}{A numeric vector with the vertical lines representing the boundaries of each window.}
#' }
#'
#'
#' @export
getFeatWindows <- function(SerObj, feat, group.by){
  
  
  
  v1 <- Seurat:::VlnPlot(SerObj, features = feat, group.by = group.by, flip = T) #timepoint
  
  v1 <- v1$data
  colnames(v1) = c("Feat", "Grp")
  
  # v1$Grp = factor(as.character(v1$Grp), levels = c("PBMC_Baseline", "PBMC_Peak", "LN_Baseline", "LN_Peak"))
  
  tmpTbl = table(cut(v1$Feat, c(min(v1$Feat),
                                # round(mean(v1$Feat) - 2*sd(v1$Feat), 2),
                                round(mean(v1$Feat) - sd(v1$Feat), 2),
                                round(mean(v1$Feat), 2),
                                round(mean(v1$Feat) + sd(v1$Feat), 2),
                                # round(mean(v1$Feat) + 2*sd(v1$Feat), 2),
                                max(v1$Feat))), v1$Grp)
  
  if(any(rowSums(tmpTbl)==0)){
    
    if(abs(round(min(v1$Feat), 2)) - abs(round(mean(v1$Feat) - sd(v1$Feat), 2)) < 0.5){
      tmpTbl = table(cut(v1$Feat, c(min(v1$Feat),
                                    # round(mean(v1$Feat) - 2*sd(v1$Feat), 2),
                                    # round(mean(v1$Feat) - sd(v1$Feat), 2),
                                    round(mean(v1$Feat), 2),
                                    # ThrX,
                                    round(mean(v1$Feat) + sd(v1$Feat), 2),
                                    round(mean(v1$Feat) + 2*sd(v1$Feat), 2),
                                    max(v1$Feat))), v1$Grp)
    } else {
      tmpTbl = table(cut(v1$Feat, c(min(v1$Feat),
                                    # round(mean(v1$Feat) - 2*sd(v1$Feat), 2),
                                    round(mean(v1$Feat) - sd(v1$Feat), 2),
                                    round(mean(v1$Feat), 2),
                                    # ThrX,
                                    round(mean(v1$Feat) + sd(v1$Feat), 2),
                                    # round(mean(v1$Feat) + 2*sd(v1$Feat), 2),
                                    max(v1$Feat))), v1$Grp)
    }
    
    
    
  }
  
  if(any(rowSums(tmpTbl)==0)){
    
    tmpTbl = table(cut(v1$Feat, c(min(v1$Feat),
                                  round(mean(v1$Feat) - 2*sd(v1$Feat), 2),
                                  round(mean(v1$Feat) - sd(v1$Feat), 2),
                                  round(mean(v1$Feat), 2),
                                  # ThrX,
                                  # round(mean(v1$Feat) + sd(v1$Feat), 2),
                                  # round(mean(v1$Feat) + 2*sd(v1$Feat), 2),
                                  max(v1$Feat))), v1$Grp)
  }
  
  vlines = lapply( strsplit(rownames(tmpTbl), ","), function(x){
    c(as.numeric(gsub("\\(", "", x[1])),
      as.numeric(gsub("\\]", "", x[2])))
  }) %>% unlist() %>% unique()
  
  return(list(winDF = tmpTbl, vlines = vlines))
  
}



#' Generate a FeaturePlot_cust which is a scatter plot with gene expression as color scale
#'
#' This function generates a scatter plot of the specified dimensions (features) of a Seurat object,
#' colored by expression values of the specified feature.
#'
#' @param SerObj A Seurat object containing the data to plot
#' @param dim1 Name of the first feature to plot on the x-axis
#' @param dim2 Name of the second feature to plot on the y-axis
#' @param pt.size Point size for the plot (default: 1)
#' @param group.by Name of the grouping variable for the plot (optional)
#' @param Feat Name of the feature to use for coloring the plot
#' @param colors Vector of colors for the color gradient (default: c('white', 'gray', 'gold', 'red', 'maroon'))
#' @param base_size Base font size for the plot (default: 14)
#'
#' @return A ggplot object
#' 
#' @export
FeaturePlot_cust <- function(SerObj, dim1, dim2, pt.size = 1, 
                             group.by, Feat, colors=c('white', 'gray', 'gold', 'red', 'maroon'), 
                             base_size = 14) {
  
  gg1 <- suppressMessages(Seurat::FeatureScatter(SerObj, feature1 = dim1, 
                                                 feature2 = dim2, 
                                                 group.by = group.by))
  
  
  gg1$data$GeneExpr <- FetchData(SerObj, Feat, layer = "data")[rownames(gg1$data), 1]
  
  gg1$data <- gg1$data[order(gg1$data$GeneExpr, decreasing = F), ]
  
  ggplot(gg1$data, aes(x = .data[[dim1]], y = .data[[dim2]], color = GeneExpr)) +
    geom_point(size = pt.size) + 
    theme_classic(base_size = base_size) +
    scale_colour_gradientn(colours = colors) + 
    ggtitle(Feat)
  
}

#' Perform unsupervised clustering and UMAP reduction 
#'
#' This function performs unsupervised clustering on single-cell RNA-seq data and then it does the Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction method. It returns and updated seurat object.
#'
#' @param SerObj A Seurat obj
#' @param Feats A vector of character strings specifying which features to use for clustering. If NULL, all features are used.
#' @param reduction The type of dimensionality reduction to apply. If NULL, no reduction is applied. Options are "PCA", "TSNE", "UMAP", or "none".
#' @param dims The number of dimensions to use for the reduction. If NULL, the default number of dimensions for the chosen method is used.
#' @param reSerObjlution The UMAP reSerObjlution parameter.
#' @param verbose A boolean indicating whether to print progress updates during the clustering process.
#' @param random.seed The random seed to use for reproducibility.
#' @param doUMAP Logical if F unsupervised clustering is only done default (T)
#' @param doClust Logical if F unsupervised clustering is only done default (T)
#' @param n.neighbors The number of neighbors to use for UMAP construction.
#' @param n.components The number of umap components default 2, or 3 or try 1
#' @param n.epochs The number of epochs to use for UMAP construction.
#' @param min.dist The minimum distance between UMAP points.
#' @param spread The spread parameter for UMAP.
#' @param reduction.key A character string specifying the prefix for the output object names.
#'
#' @return seurat obj
#'
#' @export
UnsupClust_UMAP_proc <- function(SerObj, 
                                 doClust = T,
                                 Feats = NULL, #character names
                                 reduction = NULL,
                                 dims = NULL, #numeric 
                                 reSerObjlution = 0.4, 
                                 verbose = F,
                                 random.seed = 1234,
                                 doUMAP = T, 
                                 n.components = 2,
                                 n.neighbors = 60, n.epochs = 300,
                                 min.dist = 0.2, spread = 1,
                                 reduction.key = 'UMAPCust_'){
  if(!doClust & !doUMAP ) stop("choose SerObjmething to do cluster or umap")
  
  if(doClust){
    print("Finding Neighbors")
    SerObj <- FindNeighbors(SerObj, 
                               reduction = reduction, 
                               dims = dims, 
                               verbose = verbose)
    
    print("Finding Clusters")
    SerObj <- FindClusters(object = SerObj,
                              reSerObjlution = reSerObjlution,
                              verbose = verbose,
                              random.seed = random.seed)
    
    print("Updating Ser Obj with:")
    print(paste0("ClusterNames", reduction, "_", reSerObjlution))
    
    SerObj[[paste0("ClusterNames", reduction, "_", reSerObjlution)]] <- Idents(object = SerObj)
    
  }
  
  
  
  if(doUMAP){
    UMAPDF = scCustFx:::RunUMAP.Matrix(DGEmat = SerObj@meta.data[, Feats],
                                       n_threads = 8,
                                       assay = NULL,
                                       n.neighbors = n.neighbors, #40
                                       n.components = n.components,
                                       metric = "cosine",
                                       n.epochs = n.epochs,
                                       learning.rate = 1.0,
                                       min.dist = min.dist,
                                       spread = spread,
                                       set.op.mix.ratio = 1.0,
                                       local.connectivity = 1L,
                                       repulsion.strength = 1,
                                       negative.sample.rate = 5,
                                       a = NULL,
                                       b = NULL,
                                       seed.use = random.seed,
                                       metric.kwds = NULL,
                                       angular.rp.forest = FALSE,
                                       reduction.key = reduction.key,
                                       verbose = TRUE)
    
    
    obID = paste0("umap_", reduction, "keep")
    
    
    
    SerObj[[obID]] = Seurat::CreateDimReducObject(
      embeddings = UMAPDF, 
      key = paste0(tolower(gsub("_", "", obID)), "_"), 
      assay = NULL, global = TRUE)
  }
  
  
  return(SerObj)
  
}

#' BiSplit_DE
#'
#' A function to perform differential expression analysis using BiSplit algorithm.
#'
#' @param SerObj An object of class Seurat
#' @param gate_feat A character vector of feature names to use as gatekeepers for the BiSplit algorithm.
#' @param gate_thr The threshold value for the gatekeepers. Default is 0.
#' @param plot_DimPlot Logical value indicating whether to generate a DimPlot of the results. Default is FALSE.
#' @param doDE set doDE to F from T (default) when just visualizing the gating dimplot
#' @param cols A character vector of colors to use for plotting the results. Default is c("gray", "red").
#' @param raster Logical value indicating whether to rasterize the plot. Default is FALSE.
#' @param assay The assay type to use for the analysis. Default is "RNA".
#' @param layer The layer to use for the analysis. Default is "data".
#' @param only.pos Logical value indicating whether to only consider positive differential expression. Default is FALSE.
#' @param min.pct The minimum percentage of cells expressing a feature to consider it. Default is 0.65.
#' @param min.diff.pct The minimum percentage difference in expression between groups to consider a feature differentially expressed. Default is 0.2.
#' @param logfc.threshold The log-fold change threshold to use for identifying differentially expressed features. Default is 0.25.
#'
#' @return A DF containing the results of the differential expression analysis.
#'
#' @export
BiSplit_DE <- function(SerObj, gate_feat = NULL, gate_thr = 0, 
                       plot_DimPlot = F, doDE = T,
                       cols = c("gray", "red"), raster = F, 
                                assay = "RNA", layer = "data",
                                only.pos = F, min.pct = 0.65, min.diff.pct = 0.2,
                                logfc.threshold = 0.25){
  
  if(is.null(gate_feat)) stop("gate_feat is null")
  
  SerObj$DE_bool = ifelse(SerObj[[gate_feat]] >= 0, "R", "L") #right or left of
  
  
  if(plot_DimPlot) {
    print(DimPlot(SerObj, group.by = "DE_bool",
                  cols = cols, 
                  label = F, label.size = 6, repel = T, raster = raster) + NoLegend())
    
  } else {
    if(!doDE) {
      print("plot_DimPlot is F doDE = F :: setting plot_DimPlot = T")
      print(DimPlot(SerObj, group.by = "DE_bool",
                    cols = cols, 
                    label = F, label.size = 6, repel = T, raster = raster) + NoLegend())
    }
   
    
  }
  
  if(doDE) {
    DefaultAssay(SerObj)  = assay
    Idents(SerObj) = "DE_bool"
    
    
    DEgenes_unsupclusts = FindMarkers(SerObj,
                                      logfc.threshold = logfc.threshold,
                                      ident.1 = "R",
                                      ident.2 = "L",
                                      test.use = "wilcox", 
                                      layer = layer,
                                      min.pct = min.pct, min.diff.pct = min.diff.pct,
                                      only.pos = only.pos)
    
    DEgenes_unsupclusts$PctDiff <- abs(DEgenes_unsupclusts$pct.1 - DEgenes_unsupclusts$pct.2)
    
    return(DEgenes_unsupclusts)
  }
  
  
}



#' Plot Violin Plots of Gene Expression by Group combo version
#' 
#' @param SerObj A Seurat object containing the data to plot.
#' @param Feats_pos The name of the feature to plot from pos side
#' @param Feats_neg The name of the feature to plot from neg side
#' @param group.by The name of the metadata field to group cells by.
#' @param NumFeatName The name of the variable containing the number of features for each cell. If provided, cells with a number of features greater than the median will be excluded from the plot.
#' @param cutGT  If TRUE (default), cells with a number of features greater than the median will be excluded from the plot. If FALSE, cells with a number of features less than or equal to the median will be excluded from the plot.
#' @param ThrCut  The threshold for feature filtering. Cells with a value in the sdaCompName column greater than this threshold will be excluded from the plot. Default is 0, which will skip this step.
#' @param GrpFacLevels  The levels of the factor variable to use for grouping. If NULL (default), all levels will be included.
#' @param compariSerObjns  A list of compariSerObjns to perform using \code{stat_compare_means}. Each compariSerObjn should be a vector of two group names.
#' @param xlab (Optional) The x-axis label.
#' @param ylab (Optional) The y-axis label.
#' @param palette  The color palette to use for the plot.
#' @param addJitter  If TRUE, jitter points will be added to the plot. Default is FALSE.
#' @param getWindows  Bool default T to get bins to split 
#' @param decreasing  prder of pseudotime .
#' 
#' @return A ggplot object.
#' 
#' 
#' @export
Plot_Pseudotime_V_Gene_combo <- function(SerObj, SerObjrtByName, Feats_pos = NULL, Feats_neg = NULL, base_size = 20, col_vector = col_vector,
                                         showScatter = TRUE, downsampleScatter = TRUE, scatterAlpha = 0.5, 
                                         getWindows = T, decreasing=F,  group.by = "Population") {
  
  if (is.null(Feats_pos) && is.null(Feats_neg))
    stop("Enter values for Feats_pos or Feats_neg")
  
  tempDF <- Seurat::FetchData(SerObj, SerObjrtByName)
  tempDF$orig.ord <- 1:nrow(tempDF)
  tempDF <- tempDF[order(tempDF[, 1], decreasing = decreasing), , drop = FALSE]
  tempDF$ord <- 1:nrow(tempDF)
  
  if(getWindows){
    
    vlines = scCustFx::getFeatWindows( SerObj, 
                                       feat = SerObjrtByName, 
                                       group.by = group.by)$vlines
    tempDF$bin = "low"
    
    tempDF$bin[tempDF[,1]>vlines[3]] = "mid-high"
    tempDF$bin[tempDF[,1]>vlines[4]] = "high"
    
    tempDF$bin[tempDF[,1]<=vlines[3]] = "mid-low"
    tempDF$bin[tempDF[,1]<=vlines[2]] = "low"
    tempDF$bin = factor(tempDF$bin, levels=c("low", "mid-low", "mid-high", "high"))
    
    
  }
  
  
  
  
  
  baseColnames <- colnames(tempDF)
  
  if (!is.null(Feats_pos)) {
    if (length(Feats_pos) == 1) {
      tempDF <- cbind(tempDF, Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats_pos, rownames(tempDF)])
    } else if (length(Feats_pos) > 1) {
      
      tempDF <- cbind(tempDF, Matrix::as.matrix(Matrix::t(Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats_pos, rownames(tempDF)])))
      
      # colnames(tempDF) <- c(baseColnames, Feats_pos)
    }
  }
  
  if (!is.null(Feats_neg)) {
    if (length(Feats_neg) == 1) {
      tempDF <- cbind(tempDF, Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats_neg, rownames(tempDF)])
    } else if (length(Feats_neg) > 1) {
      
      tempDF <- cbind(tempDF, Matrix::as.matrix(Matrix::t(Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats_neg, rownames(tempDF)])))
      # colnames(tempDF) <- c(baseColnames, Feats_neg)
    }
  }
  
  tempDF_long <- reshape2::melt(tempDF, id.vars = baseColnames, measure.vars = c(Feats_pos, Feats_neg))
  
  tempDF_long$feat_type <- ifelse(tempDF_long$variable %in% Feats_pos, "Positive", "Negative")
  
  
  gg1 <- ggplot(tempDF_long, aes(x = sqrt(ord), y = value, color = variable, linetype = feat_type))
  
  if (showScatter) {
    if (downsampleScatter) {
      gg1 <- gg1 + geom_point(data = tempDF_long %>% sample_frac(0.2), alpha = scatterAlpha)
    } else {
      gg1 <- gg1 + geom_point(alpha = scatterAlpha)
    }
  }
  
  
  
  gg1 = gg1 + geom_smooth(aes(fill = feat_type), method = "auto", se = TRUE, level = 0.95) +
    scale_color_manual(values = c(col_vector[1:length(Feats_pos)], col_vector[1:length(Feats_neg)]), name = "Top Genes") +
    scale_linetype_manual(values = c("dashed", "SerObjlid"), name = "Feature Type", guide = "none") +
    scale_fill_manual(values = c("skyblue2", "pink2"), name = "Feature Type", guide = "none") +
    theme_classic(base_size = base_size) +
    ggtitle(SerObjrtByName) +
    # scale_y_sqrt()+
    xlab("Sqrt-ordered score high -> low") +
    ylab("Gene Expression") 
  
  if(getWindows){
    gg1 = gg1 +  geom_vline(xintercept = c(
      # sqrt(min(subset(tempDF, bin=="low")$ord)),
      sqrt(min(subset(tempDF, bin=="mid-low")$ord)),
      sqrt(min(subset(tempDF, bin=="mid-high")$ord)),
      sqrt(min(subset(tempDF, bin=="high")$ord))
    ))
    
  }
  
  return(gg1)
}



#' Plot Violin Plots of Gene Expression by Group
#' 
#' @param SerObj A Seurat object containing the data to plot.
#' @param Feature The name of the feature to plot.
#' @param group.by The name of the metadata field to group cells by.
#' @param NumFeatName The name of the variable containing the number of features for each cell. If provided, cells with a number of features greater than the median will be excluded from the plot.
#' @param cutGT  If TRUE (default), cells with a number of features greater than the median will be excluded from the plot. If FALSE, cells with a number of features less than or equal to the median will be excluded from the plot.
#' @param ThrCut  The threshold for feature filtering. Cells with a value in the sdaCompName column greater than this threshold will be excluded from the plot. Default is 0, which will skip this step.
#' @param GrpFacLevels  The levels of the factor variable to use for grouping. If NULL (default), all levels will be included.
#' @param compariSerObjns  A list of compariSerObjns to perform using \code{stat_compare_means}. Each compariSerObjn should be a vector of two group names.
#' @param xlab (Optional) The x-axis label.
#' @param ylab (Optional) The y-axis label.
#' @param palette  The color palette to use for the plot.
#' @param addJitter  If TRUE, jitter points will be added to the plot. Default is FALSE.
#' 
#' @return A ggplot object.
#' 
#' 
#' @export
plot_violin_wFeatFilter <- function(SerObj, Feature, group.by, 
                                    NumFeatName = NULL, 
                                    cutGT = T,
                                    ThrCut = 0, 
                                    GrpFacLevels=NULL,
                                    compariSerObjns = NULL, xlab="", ylab="", 
                                    palette = col_vector, addJitter = F) {
  
  
  v1 <- Seurat:::VlnPlot(SerObj, 
                         features = Feature, 
                         group.by = group.by, flip = T) 
  
  scoreDF <- SerObj[[NumFeatName]]
  
  if(cutGT) {
    KeepCells <- scoreDF[scoreDF[,1] > ThrCut, , drop=F] %>% rownames()
  } else {
    KeepCells <- scoreDF[scoreDF[,1] < ThrCut, , drop=F] %>% rownames()
    
  }
  
  v1 <- v1$data[KeepCells,]
  colnames(v1) <- c("Feat", "Grp")
  
  if (!is.null(GrpFacLevels)) {
    v1$Grp <- factor(as.character(v1$Grp), levels = GrpFacLevels)
  }
  # nrow(v1)
  
  
  if(addJitter) {
    ggv <- ggpubr::ggviolin(v1, x = "Grp", y = "Feat",
                            fill = "Grp", palette = palette,
                            add = "jitter",
                            add.params = list(size = .005, alpha = 0.05),
                            title = paste0(Feature, "\nWilcox p."))
  } else {
    ggv <- ggpubr::ggviolin(v1, x = "Grp", y = "Feat",
                            fill = "Grp", palette = palette,
                            # add = "jitter",
                            add.params = list(size = .005, alpha = 0.05),
                            title = paste0(Feature, "\nWilcox p."))
  }
  
  
  if(!is.null(compariSerObjns)){
    ggv = ggv + ggpubr::stat_compare_means(label = "p.signif",
                                           compariSerObjns =
                                             compariSerObjns )
  }
  
  ggv = ggv + xlab(xlab) + ylab(ylab) + 
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1.05, hjust=1))
  
  return(ggv)
}







#' Plot a scatter plot with smoothed lines for selected features
#'
#' This function creates a scatter plot of a Seurat object, with smoothed lines for each selected feature. The x-axis is defined by the selected SerObjrtByName parameter, and the y-axis is defined by the values of the selected features added to the data frame. The function uses ggplot2 for visualization.
#'
#' @param SerObj A Seurat object to plot.
#' @param SerObjrtByName The name of the variable to use for SerObjrting cells on the x-axis.
#' @param Feats A vector of feature names to plot on the y-axis. Default is NULL, which will plot all features in the object.
#' @param base_size Base size of points and lines. Default is 20.
#' @param col_vector color vector
#' @param showScatter boolean default T to show scatter points
#' @param dowsampleScatter boolean default T to downsample the scatter points
#' @param scatterAlpha alpha value default 0.5 to show scatter points
#'
#' @return A ggplot object.
#'
#' @examples
#' Plot_Pseudotime_V_Gene(SerObj, "nCount_RNA", c("TNFRSF4", "CD69"))
#'
#' @export
Plot_Pseudotime_V_Gene <- function(SerObj, SerObjrtByName, Feats=NULL, base_size = 20, col_vector=col_vector,
                                   showScatter = T, dowsampleScatter = T, scatterAlpha = 0.5) {
  
  if(is.null(Feats)) stop("Enter values for Feats")
  
  tempDF = Seurat::FetchData(SerObj, SerObjrtByName)
  tempDF$orig.ord = 1:nrow(tempDF)
  tempDF = tempDF[order(tempDF[,1], decreasing = T), , drop=F]
  tempDF$ord = 1:nrow(tempDF)
  
  baseColnames = colnames(tempDF)
  # head(tempDF)
  
  if(length(Feats) == 1){
    
    tempDF = cbind(tempDF, Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats, rownames(tempDF)])
    
  } else if(length(Feats) > 1){
    
    tempDF = cbind(tempDF, Matrix::as.matrix( Matrix::t(Seurat::GetAssayData(SerObj, layer = "data", assay = "RNA")[Feats, rownames(tempDF)]) ) )
    
    colnames(tempDF) = c(baseColnames, Feats)
    
  }
  
  # head(tempDF)
  
  # Create smoothed lines for each item in Feat
  tempDF_long <- reshape2::melt(tempDF, 
                                id.vars = baseColnames, 
                                measure.vars = Feats)
  head(tempDF_long)
  
  gg1 = ggplot(tempDF_long, aes(x = sqrt(ord), y = value, color = variable))
  
  if(showScatter) {
    if(dowsampleScatter){
      gg1 = gg1 + geom_point(data = tempDF_long %>% sample_frac(0.2), alpha = scatterAlpha)
    } else {
      gg1 = gg1 + geom_point(alpha = scatterAlpha)
    }
  }
  
  gg1 + geom_smooth(method = "auto", se = TRUE, level = 0.95) +
    scale_color_manual(values = col_vector, name = "Feature")+
    theme_bw(base_size = base_size) + 
    ggtitle(SerObjrtByName) + 
    xlab("Sqrt-ordered score high -> low") + 
    ylab ("Gene Expression")
  
  
  
  # # Create a factor variable with 3 bins
  # tempDF_long2 <- tempDF_long #%>% mutate(ord_bin = cut(sqrt(ord), breaks = 3, labels = c("A", "B", "C")))
  # 
  # # Create a factor variable with 3 bins based on the order variable
  # tempDF_long2$bin <- cut(tempDF_long2$ord, breaks = 3, labels = c("A", "B", "C"), include.lowest = TRUE)
  # 
  # tempDF_long2 = subset(tempDF_long2, value > 0)
  # 
  # # tempDF_long2$value = exp(tempDF_long2$value )^2
  # 
  # library(ggpubr)
  # # Add pairwise compariSerObjns p-value
  # my_compariSerObjns <- list(c("A", "B"), c("B", "C"), c("A", "C"))
  # 
  # 
  # # Create a box plot with colored segments
  # ggboxplot(tempDF_long2, x = "bin", y = "value", color = "variable", 
  #                        palette = col_vector,#  c("#E69F00", "#56B4E9", "#009E73"),
  #                        ggtheme = theme_classic()) +
  #   ggtitle("Gene Expression by Score Bins") +
  #   xlab("Score Bins") +
  #   ylab("Gene Expression") #+ 
  #  # stat_compare_means(compariSerObjns = my_compariSerObjns, method = "wilcox.test",
  #  #                              label = "p.signif"
  #  #                              #label.y = max(tempDF_long$value) * 1.05
  #  #                              )
  
  
}



#' BetterViolins
#'
#' This function generates violin plots with customizable features.
#'
#' @param SerObj a data frame or a Seurat object containing the expression data
#' @param feature2plot a string indicating the feature to plot on the y-axis (default: "nCount_RNA")
#' @param featuredummy a string indicating a dummy feature to plot on the x-axis (default: "nFeature_RNA")
#' @param group.by a string indicating the grouping variable (default: "ExpID")
#' @param downsampleN an integer indicating the number of cells to downsample to (default: 500)
#' @param title a string indicating the title of the plot (default: "")
#' @param ComplexPlot a boolean indicating whether to generate a complex plot or not (default: FALSE)
#' @param log10Y a boolean indicating whether to use log10 transformation on the y-axis (default: TRUE)
#'
#' @return a ggplot object
#' @export
BetterViolins = function(SerObj=NULL, 
                         feature2plot = "nCount_RNA", 
                         featuredummy = "nFeature_RNA", 
                         group.by = "ExpID", 
                         downsampleN = 500, title="", ComplexPlot=F, log10Y=T){
  
  if(is.null(SerObj)) stop()
  
  if(!ComplexPlot) downsampleN = "none" # no need to downsample since ggstatplot is slow and needs downsampling 
  if(ComplexPlot) require(ggstatsplot)
  
  
  
  ggScat1 = FeatureScatter(SerObj, feature1 =feature2plot, feature2 = featuredummy, group.by = group.by)
  
  plotLS = lapply(levels(ggScat1$data$colors), function(xL){
    #xL = "T1"
    
    tempDF = subset(ggScat1$data[,c(feature2plot, "colors")], colors == xL)
    
    
    
    labs = reshape2::melt(tempDF) %>%
      mutate(text = forcats::fct_reorder(colors, value)) %>%
      group_by(variable, colors) %>%
      summarize(mean = mean(value, na.rm = TRUE),
                median = median(value, probs=5))
    
    
    if(!ComplexPlot) {
      
      
      
      gg <- ggplot(reshape2::melt(tempDF),
                   aes(x=colors, y=value)) +
        # geom_jitter(color="black", aes(size=log10(value)), alpha=0.2) +
        # geom_violin( size=0.2) +
        geom_boxplot(color = "black", size = .2, alpha = .6) +
        # scale_y_log10()+
        ggrepel::geom_label_repel(data = labs, aes(colors, mean,  label=paste0("mean = ", round(mean))),
                                  nudge_x = .15,
                                  box.padding = 0.5,
                                  nudge_y = 1,
                                  segment.curvature = -0.1,
                                  segment.ncp = 3,
                                  segment.angle = 20) + 
        scale_fill_viridis(discrete=TRUE) +
        scale_color_viridis(discrete=TRUE) +
        theme_bw() +
        theme(
          legend.position="none"
        ) +
        ggrepel::geom_label_repel(data = labs, aes(colors, median,  label=paste0("median = ", round(median))),
                                  nudge_x = -.15,
                                  box.padding = 0.5,
                                  nudge_y = 1,
                                  segment.curvature = -0.1,
                                  segment.ncp = 3,
                                  segment.angle = 20) + 
        # scale_fill_viridis(discrete=TRUE) +
        # scale_color_viridis(discrete=TRUE) +
        theme_bw() +
        theme(
          legend.position="none"
        ) +
        coord_flip() +
        xlab("") +
        ylab("")
      
      if(log10Y) gg = gg + scale_y_log10()
      
    }
    
    if(ComplexPlot) {
      if(!is.null(downsampleN)){
        if(nrow(tempDF)>downsampleN){
          tempDF = tempDF[sample(1:nrow(tempDF), downsampleN),]
        }}
      
      gg <- ggbetweenstats(
        data  = reshape2::melt(tempDF),
        x     = colors,
        y     = value,
        title = "",
        xlab = "",
        ylab = "")
      
    }
    
    gg
    
    
  })
  
  patchwork::wrap_plots(plotLS) + 
    plot_annotation(
      title = title,
      subtitle = feature2plot,
      caption = paste0("downsampled N=", downsampleN, ifelse(log10Y, " | log10 scale", ""))
    )
  
}




#' scatter-bubble style of selected PCs and genes
#'
#'scatter-bubble style of selected PCs and genes
#' 
#' @param SerObj a seurat obj 
#' @param features features to show
#' @param pcs a numeric vector of 2 pc nnumbers e.g. c(1, 2) 
#' @param base_size base_size of plot
#' @param quaT color pink via quantiles
#' @return ggplot scatter-bubble plot style 
#' @export
plotPCgeneLoadings = function(SerObj, features, pcs = c(1, 2), quaT = 0.8, base_size =10){
  tempDF = as.data.frame(SerObj@reductions$pca@feature.loadings[features, pcs])
  colnames(tempDF) = c("Dim1", "Dim2")
  tempDF$gene = rownames(tempDF)
  
  tempDF$col = "gray"
  
  # tempDF[abs(tempDF$Dim1) > quantile(abs(tempDF$Dim1), quaT),]$col = "pink"
  # tempDF[abs(tempDF$Dim2) > quantile(abs(tempDF$Dim2), quaT),]$col = "pink"
  
  tempDF[(tempDF$Dim1) > quantile((tempDF$Dim1), quaT),]$col = "pink"
  tempDF[(tempDF$Dim1) < quantile((tempDF$Dim1), 1-quaT),]$col = "pink"
  
  tempDF[(tempDF$Dim2) > quantile((tempDF$Dim2), quaT),]$col = "pink"
  tempDF[(tempDF$Dim2) < quantile((tempDF$Dim2), 1-quaT),]$col = "pink"
  
  
  tempDF$lab = ifelse(tempDF$col == "gray", "", tempDF$gene)
  
  tempDF$col <- factor(tempDF$col, levels = c("pink", "gray"))
  
  tempDF$meanW = (as.numeric(tempDF$Dim1) + as.numeric( tempDF$Dim2))/2
  
  ggplot(tempDF, aes(Dim1, Dim2, col=col))+
    geom_point(aes(size = meanW)) + 
    theme_bw(base_size = base_size) +
    # scale_color_distiller(palette = "Spectral") +
    theme(legend.position = "none") +
    scale_color_manual(values = c("pink", "gray")) +
    ggrepel::geom_text_repel(aes(label=lab), 
                             size=4, max.overlaps =10,
                             colour = "black",
                             # segment.color="grey50", 
                             nudge_y = 0.01) +
    ggtitle("Top loaded genes")+
    xlab(paste0("PC", pcs[1])) + ylab(paste0("PC", pcs[2]))
  
  
  
}







#' Calculate the percent expressed cells for each group in a Seurat object
#'
#' @param SerObj A Seurat object
#' @param group.by A meta feature used to group the cells
#' @param features A vector of features (genes) to calculate the percent expressed cells for
#' @param plot_perHM A logical indicating whether to plot a heatmap of the results
#' @return A data frame with the percent expressed cells for each group and feature
#' @export
PercentExpressed <- function(SerObj, group.by, features, plot_perHM = F){
  
  print("removing:")
  print(features[!features %in% rownames(SerObj)])
  
  
  features = features[features %in% rownames(SerObj)]
  MyDot = DotPlot(SerObj, group.by = group.by, features = (features) )
  
  PctExprDF = lapply(levels(MyDot$data$id), function(xL){
    subset(MyDot$data, id == xL)$pct.exp
  }) %>% as.data.frame()
  
  rownames(PctExprDF) = features
  colnames(PctExprDF) = levels(MyDot$data$id)
  rowSums(PctExprDF)
  
  if(plot_perHM) pheatmap::pheatmap(PctExprDF)
  
  return(PctExprDF)
}


#' Make an RNA heatmap
#' 
#' @param SerObj A Seurat object.
#' @param markerVec A vector of marker genes.
#' @param labelColumn Column name to use for labeling rows.
#' @param rowsplit Column name to split rows by.
#' @param columnsplit Column name to split columns by.
#' @param size Size of plot.
#' @param coldendside Dendrogram on the column side.
#' @param rowdendside Dendrogram on the row side.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param asGG returns GGplot not list if true
#' @return A heatmap of RNA expression.
#' @export
make_RNA_heatmap = function(SerObj, markerVec, labelColumn, rowsplit=NULL, columnsplit=NULL,
                            size=6,
                            clus_cols = TRUE, show_column_dend = T,
                            clus_rows=TRUE, show_row_dend = T,
                            fontsize=10, titlefontsize=10, pairedList2=NULL, asGG = F,
                            column_names_side = "bottom",
                            column_dend_side = "bottom",
                            row_names_side = "left",
                            row_dend_side = "left",
                            assays = 'RNA', useCDnames = T){
  avgSeurat <- Seurat::AverageExpression(SerObj, group.by = labelColumn,
                                         features = markerVec,
                                         layer = 'counts', return.seurat = T,
                                         assays = assays)
  avgSeurat <- Seurat::NormalizeData(avgSeurat)
  mat <- t(as.matrix(Seurat::GetAssayData(avgSeurat, layer = 'data')))
  mat <- mat %>% pheatmap:::scale_mat(scale = 'column')
  
  if(useCDnames) colnames(mat) <- CellMembrane::RenameUsingCD(colnames(mat))
  
  if(!is.null(pairedList2)){
    for (nam in names(pairedList2)){
      foundnam_pos <- grep(nam, colnames(mat))
      rep <- pairedList2[[nam]]
      colnames(mat)[foundnam_pos] <- rep
    }
  }
  
  col_RNA = circlize::colorRamp2(c(min(mat), 0, max(mat)), c(Seurat::BlueAndRed(20)[1], "gray85", Seurat::BlueAndRed(20)[20]), space = "sRGB")
  # col_RNA = c(Seurat::BlueAndRed(20)[c(1,3,5,7)], Seurat::BlueAndRed(20)[c(14,16,18,20)])
  #Above emulates Seurat's BlueAndRed color scheme.
  P1 <- ComplexHeatmap::Heatmap(mat,
                            width = ncol(mat)*unit(size, "mm"), 
                            height = nrow(mat)*unit(size, "mm"),
                            row_names_side = row_names_side,
                            row_dend_side = row_dend_side,,
                            column_names_rot = 45,
                            col = col_RNA,
                            column_names_side = column_names_side,
                            column_dend_side = column_dend_side,
                            row_split = rowsplit,
                            column_split = columnsplit,
                            row_title=NULL,
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_names_gp = grid::gpar(fontsize = fontsize),
                            column_title = "RNA Markers\n", 
                            column_title_gp = grid::gpar(fontsize = titlefontsize), 
                            name = "Scaled Avg. Exp.",
                            cluster_columns = clus_cols,
                            cluster_rows = clus_rows,
                            show_row_dend = show_row_dend, 
                            show_column_dend = show_column_dend, 
                            show_heatmap_legend = FALSE
    )
  
  if(asGG) {
    P1 <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap:::draw(P1))) 
    return(P1)
  } else {
    return(list(col_RNA = col_RNA, plot = P1))
    }
  
 
}


#' Plot RNA heatmap v2, no re-normalization of the data
#' 
#' @param SerObj A Seurat object.
#' @param labelColumn Column for cell labels.
#' @param markerVec Marker genes to use.
#' @param assay Assay type.
#' @param clus_cols Cluster columns.
#' @param show_column_dend Show column dendrogram.
#' @param clus_rows Cluster rows.
#' @param show_row_dend Show row dendrogram.
#' @param rowsplit Split rows.
#' @param columnsplit Split columns.
#' @param show_heatmap_legend Show heatmap legend.
#' @param size Plot size.
#' @param column_names_rot Rotation for column names.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param coldendside Position of column dendrogram.
#' @param asGG Output as ggplot object.
#' @return A heatmap of RNA expression.
#' @export
make_RNA_heatmap2 = function(SerObj, labelColumn = 'ClusterNames_0.2', 
                            markerVec = NULL, assay = "RNA",
                            clus_cols = TRUE, show_column_dend = T,
                            clus_rows=TRUE, show_row_dend = T,
                            rowsplit = NULL, columnsplit=NULL, 
                            show_heatmap_legend = T,
                            size = 4, column_names_rot = 45,
                            fontsize=8, titlefontsize = 12,
                            coldendside = "bottom",
                            column_names_side = "bottom",
                            row_names_side = "left",
                            row_dend_side = "right",
                            asGG = T){
  
  
  mat <- asinh(scale(t(AverageExpression(SerObj, group.by = labelColumn, features = markerVec)[[assay]])))
  col_RNA = circlize::colorRamp2(c(min(mat), 0, max(mat)), c(Seurat::BlueAndRed(20)[1], "gray85", Seurat::BlueAndRed(20)[20]), space = "sRGB")
  
  P1 <-
    ComplexHeatmap::Heatmap(mat,
                            row_names_side = row_names_side,
                            row_dend_side = row_dend_side,
                            col = col_RNA,
                            column_names_side = column_names_side,
                            column_names_side = column_names_side,
                            column_dend_side = column_dend_side,
                            row_split = rowsplit,
                            column_split = columnsplit,
                            width = ncol(mat)*unit(size, "mm"),
                            height = nrow(mat)*unit(size, "mm"),
                            column_names_rot = column_names_rot,
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_title_gp = grid::gpar(fontsize = titlefontsize),
                            # name = "Scaled Avg. Exp.",
                            show_row_dend = show_row_dend,
                            show_column_dend = show_column_dend,
                            show_heatmap_legend = show_heatmap_legend,
                            row_title=NULL,
                            cluster_columns = clus_cols,
                            cluster_rows = clus_rows,
                            column_names_gp = grid::gpar(fontsize = 10),
                            column_title = NULL
    )
  

  
  if(asGG) {
    P1 <- ggplotify::as.ggplot(grid::grid.grabExpr(ComplexHeatmap:::draw(P1))) 
    return(P1)
  } else {
    return(list(col_RNA = col_RNA, plot = P1))
  }
  
}

#' Produce Combo Heatmap From Matrix
#'
#' This function produces a combination heatmap from a matrix using Seurat object.
#'
#' @param SerObj A Seurat object containing the single-cell RNA-seq data.
#' @param mat A matrix object containing the data for generating the heatmap.
#' @param markerVec A vector of marker genes to be used for generating the heatmap.
#' @param pairedList A list of paired conditions for the first set of markers.
#' @param pairedList2 A list of paired conditions for the second set of markers.
#' @param labelColumn A column name or index in the metadata of SerObj to be used as labels for the heatmap.
#' @param prefix A prefix for the output file name.
#' @param adtoutput Output type for ADT features. Default is "unpaired".
#' @param rowsplit A vector of row splits for the heatmap.
#' @param columnsplit A vector of column splits for the heatmap.
#' @param size Size of the heatmap. Default is NULL.
#' @param coldend A logical value indicating whether to add dendrogram on the columns. Default is TRUE.
#' @param rowdend A logical value indicating whether to add dendrogram on the rows. Default is TRUE.
#' @param coldendside A character value specifying the side to place the column dendrogram. Default is "bottom".
#' @param rowdendside A character value specifying the side to place the row dendrogram. Default is "left".
#' @param fontsize Font size for the heatmap labels. Default is 12.
#' @param titlefontsize Font size for the title of the heatmap. Default is 14.
#' @param gap Gap between the heatmaps. Default is 0.
#'
#' @return A combo heatmap plot.
#' @export
ProduceComboHeatmapFromMatrix <- function(SerObj, mat, markerVec, pairedList, pairedList2, labelColumn, prefix, adtoutput = "unpaired", rowsplit = NULL, columnsplit = NULL, size, coldend = TRUE, rowdend = TRUE, coldendside = "bottom", rowdendside = "left", fontsize = 12, titlefontsize = 14, gap = 0){
  
  library(circlize)
  
  mat <- mat %>% pheatmap:::scale_mat(scale = 'row') %>% as.data.frame() %>% filter(rownames(mat) %in% markerVec)
  colnames(mat) <- CellMembrane::RenameUsingCD(colnames(mat))
  mat <- t(as.matrix(mat)) %>% as.data.frame()
  for (nam in names(pairedList2)){
    foundnam_pos <- grep(nam, colnames(mat))
    rep <- pairedList2[[nam]]
    colnames(mat)[foundnam_pos] <- rep
  }
  
  col_RNA = circlize::colorRamp2(c(min(mat), 0, max(mat)), c(Seurat::BlueAndRed(20)[1], "gray85", Seurat::BlueAndRed(20)[20]), space = "sRGB")
  # col_RNA = c(Seurat::BlueAndRed(20)[c(1,3,5,7)], Seurat::BlueAndRed(20)[c(14,16,18,20)])
  #Above emulates Seurat's BlueAndRed color scheme.
  fullorder <- rownames(mat)
  P1 <-
    ComplexHeatmap::Heatmap(mat,
                            width = ncol(mat)*unit(size, "mm"), 
                            height = nrow(mat)*unit(size, "mm"),
                            row_names_side = "left",
                            row_dend_side = rowdendside,
                            column_names_rot = 45,
                            col = col_RNA,
                            column_names_side = "bottom",
                            column_dend_side = coldendside,
                            row_split = rowsplit,
                            column_split = columnsplit,
                            row_title=NULL,
                            cluster_columns = TRUE,
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_names_gp = grid::gpar(fontsize = fontsize),
                            column_title = "RNA Markers\n", column_title_gp = grid::gpar(fontsize = titlefontsize, fontface = "bold"), name = "Scaled Avg. Exp.", show_row_dend = rowdend, show_column_dend = coldend, show_heatmap_legend = FALSE
    )
  
  
  
  # P2 <- CreatePairedHeatmap(SerObj, pairedList, labelColumn)
  if (adtoutput == "unpaired"){
    res <- CreateProteinHeatmap(SerObj, proteinList = unname(unlist(pairedList)), labelColumn=labelColumn, prefix = prefix, size, fontsize = fontsize, titlefontsize = titlefontsize, fullorder = fullorder)
    P2 <- res[[1]]
    col_ADT <- res[[2]]
  } else {
    res <- CreatePairedHeatmap(SerObj, pairedList, labelColumn=labelColumn, prefix = prefix)
    
  }
  
  
  res2 <- plotEnrichment(SerObj, field1 = labelColumn, field2 = 'Tissue', size, fontsize = fontsize, titlefontsize = titlefontsize, fullorder = fullorder)
  P3 <- res2[[1]]
  col_tiss <- res2[[2]]
  plotlist <- P1 + P2 + P3
  # draw(plotlist, ht_gap = unit(0, "mm"))
  lgd1 <- Legend(title = "Scaled Avg. Exp.", col_fun = col_RNA, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  lgd2 <- Legend(title = "Scaled Avg. ADT", col_fun = col_ADT, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  lgd3 <- Legend(title = "Tissue Enrichment", col_fun = col_tiss, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  pd = packLegend(lgd1, lgd2, lgd3, direction = "vertical")
  draw(P1)
  draw(P2)
  draw(P3)
  draw(plotlist, heatmap_legend_list = pd, ht_gap = unit(gap,"mm"), heatmap_legend_side = "right")
  
}


#' Produce a combination heatmap
#' 
#' @param SerObj A seurat object
#' @param markerVec A vector of markers
#' @param pairedList A list of paired data
#' @param pairedList2 A second list of paired data
#' @param labelColumn Column name for sample labels
#' @param prefix Prefix for output file names
#' @param adtoutput Output type (default: "unpaired")
#' @param rowsplit Split data by row (default: NULL)
#' @param columnsplit Split data by column (default: NULL)
#' @param size Plot size
#' @param coldend Show dendrogram for columns (default: TRUE)
#' @param rowdend Show dendrogram for rows (default: TRUE)
#' @param coldendside Side to place column dendrogram (default: "bottom")
#' @param rowdendside Side to place row dendrogram (default: "left")
#' @param fontsize Font size (default: 12)
#' @param titlefontsize Title font size (default: 20)
#' @param gap Gap between panels (default: 0)
#' @return A list of differentially expressed genes
#' @export
ProduceComboHeatmap <- function(SerObj, markerVec, pairedList, pairedList2, labelColumn, prefix, adtoutput = "unpaired", rowsplit = NULL, columnsplit = NULL, size, coldend = TRUE, rowdend = TRUE, coldendside = "bottom", rowdendside = "left", fontsize = 12, titlefontsize = 20, gap = 0){
  avgSeurat <- Seurat::AverageExpression(SerObj, group.by = labelColumn,
                                         features = markerVec,
                                         layer = 'data', return.seurat = T,
                                         assays = 'RNA')
  # avgSeurat <- Seurat::NormalizeData(avgSeurat)
  mat <- t(as.matrix(Seurat::GetAssayData(avgSeurat, layer = 'data')))
  mat <- mat %>% pheatmap:::scale_mat(scale = 'column') %>% as.data.frame()
  colnames(mat) <- CellMembrane::RenameUsingCD(colnames(mat))
  for (nam in names(pairedList2)){
    foundnam_pos <- grep(nam, colnames(mat))
    rep <- pairedList2[[nam]]
    colnames(mat)[foundnam_pos] <- rep
  }
  col_RNA = circlize::colorRamp2(c(min(mat), 0, max(mat)), c(Seurat::BlueAndRed(20)[1], "gray85", Seurat::BlueAndRed(20)[20]), space = "sRGB")
  # col_RNA = c(Seurat::BlueAndRed(20)[c(1,3,5,7)], Seurat::BlueAndRed(20)[c(14,16,18,20)])
  #Above emulates Seurat's BlueAndRed color scheme.
  fullorder <- rownames(mat)
  P1 <-
    ComplexHeatmap::Heatmap(mat,
                            width = ncol(mat)*unit(size, "mm"), 
                            height = nrow(mat)*unit(size, "mm"),
                            row_names_side = "left",
                            row_dend_side = rowdendside,
                            column_names_rot = 45,
                            col = col_RNA,
                            column_names_side = "bottom",
                            column_dend_side = coldendside,
                            row_split = rowsplit,
                            column_split = columnsplit,
                            row_title=NULL,
                            cluster_columns = TRUE,
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            column_names_gp = grid::gpar(fontsize = fontsize),
                            column_title = "RNA Markers\n", column_title_gp = grid::gpar(fontsize = titlefontsize, fontface = "bold"), name = "Scaled Avg. Exp.", show_row_dend = rowdend, show_column_dend = coldend, show_heatmap_legend = FALSE
    )
  
  # P2 <- CreatePairedHeatmap(SerObj, pairedList, labelColumn)
  if (adtoutput == "unpaired"){
    res <- CreateProteinHeatmap(SerObj, proteinList = unname(unlist(pairedList)), labelColumn=labelColumn, prefix = prefix, size, fontsize = fontsize, titlefontsize = titlefontsize, fullorder = fullorder)
    P2 <- res[[1]]
    col_ADT <- res[[2]]
  } else {
    res <- CreatePairedHeatmap(SerObj, pairedList, labelColumn=labelColumn, prefix = prefix)
    
  }
  
  
  res2 <- plotEnrichment(SerObj, field1 = labelColumn, field2 = 'Tissue', size, fontsize = fontsize, titlefontsize = titlefontsize, fullorder = fullorder)
  P3 <- res2[[1]]
  col_tiss <- res2[[2]]
  plotlist <- P1 + P2 + P3
  # draw(plotlist, ht_gap = unit(0, "mm"))
  lgd1 <- Legend(title = "Scaled Avg. Exp.", col_fun = col_RNA, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  lgd2 <- Legend(title = "Scaled Avg. ADT", col_fun = col_ADT, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  lgd3 <- Legend(title = "Tissue Enrichment", col_fun = col_tiss, direction = "horizontal", title_gp = gpar(fontsize = 10),
                 labels_gp = gpar(fontsize = 8))
  pd = packLegend(lgd1, lgd2, lgd3, direction = "vertical")
  draw(P1)
  draw(P2)
  draw(P3)
  draw(plotlist, heatmap_legend_list = pd, ht_gap = unit(gap,"mm"), heatmap_legend_side = "right")
}

#' Create a protein heatmap
#' 
#' @param SerObj A Seurat object.
#' @param proteinList A list of proteins to plot.
#' @param labelColumn Column name to use for labeling rows.
#' @param prefix Prefix for plot title.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param fullorder 
#' @return A heatmap of protein expression.
#' @export
CreateProteinHeatmap <- function(SerObj, proteinList, labelColumn, prefix, size, fontsize, titlefontsize, fullorder) {
  print(proteinList)
  adt_columns <- proteinList # unname(unlist(pairedList))
  mat_ADT <- cbind("Label" = SerObj@meta.data[[labelColumn]], SerObj@meta.data[,adt_columns])
  colnames(mat_ADT) <- gsub(paste0(prefix, "."), "", colnames(mat_ADT))
  colnames(mat_ADT) <- gsub("_UCell_pos", "", colnames(mat_ADT)) %>% as.factor()
  colnames(mat_ADT) <- gsub("\\.", "-", colnames(mat_ADT)) %>% as.factor()
  
  
  colordering <- colnames(mat_ADT)
  mat_ADT <- mat_ADT %>% pivot_longer(cols = 2:length(mat_ADT), names_to = "ADT") %>% group_by(Label, ADT) %>% summarize(avg = mean(value)) %>% pivot_wider(id_cols = Label, names_from = "ADT", values_from = "avg") %>% as.data.frame()
  
  colordering_mat <- colordering[-1]
  
  
  
  mat_ADT <- mat_ADT[,colordering]
  rownames(mat_ADT) <- mat_ADT$Label
  mat_ADT <- mat_ADT %>% select(-Label)
  mat_ADT <- mat_ADT[fullorder,]
  
  rowordering <- rownames(mat_ADT)
  avgSeuratADT <- Seurat::AverageExpression(SerObj, group.by = labelColumn,
                                            layer = 'data', return.seurat = T,
                                            assays = prefix)
  # avgSeuratADT <- Seurat::NormalizeData(avgSeuratADT)
  mat <- t(as.matrix(avgSeuratADT@assays[[prefix]]@data)) %>% pheatmap:::scale_mat(scale = 'column') %>% as.data.frame()
  
  mat <- mat[,colnames(mat) %in% colnames(mat_ADT)]
  mat <- mat[fullorder,]
  mat <- mat[,colordering_mat]
  
  col_ADT = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("white", "gray85", "blue"))
  
  print(mat)
  print(mat_ADT)
  
  P2 <- Heatmap(as.matrix(mat_ADT), name = "Scaled\nAvg. ADT", col = col_ADT, rect_gp = gpar(type="none"), border_gp = gpar(col = "black", lty = 1), show_row_dend = FALSE, show_column_dend = FALSE, width = ncol(mat)*unit(size, "mm"), 
                height = nrow(mat)*unit(size, "mm"), column_names_gp = grid::gpar(fontsize = fontsize), row_names_gp = grid::gpar(fontsize = fontsize),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.circle(x = x, y = y, r = (mat_ADT[i,j])* min(unit.c(width)/2), #(mat_fin[i, j])/3 * min(unit.c(width, height)), 
                              gp = gpar(fill = col_ADT(mat[i, j]), col = NA))
                }, cluster_rows = TRUE, cluster_columns = FALSE,
                column_names_side = "bottom",
                column_dend_side = NULL,
                column_names_rot = 45,
                row_dend_side = NULL,
                column_title_gp = grid::gpar(fontsize = titlefontsize, fontface = "bold"), column_title = "Surface\nProtein",
                show_row_names = FALSE, show_column_names = TRUE, show_heatmap_legend = FALSE, row_split = NULL)
  print(P2)
  return(list(P2, col_ADT))
}


#' Plot Enrichment 
#' 
#' @param SerObj A Seurat object.
#' @param field1 Field 1 to compare.
#' @param field2 Field 2 to compare.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param fullorder 
#' @return A plot of enrichment analysis between two fields.
#' @export
plotEnrichment <- function(SerObj, field1, field2, size, fontsize, titlefontsize, fullorder) {
  mat <- asinh(chisq.test(table(SerObj[[field1]][[field1]], SerObj[[field2]][[field2]]))$res) %>% as.matrix()  
  mat <- mat[fullorder,]
  
  rowLabels <- apply(mat, 1, function(x){
    # ToDo:
    lab = ifelse(x >= 0.9 * max(x), "***", ifelse(x >= 0.75 * max(x), "**", ifelse(x >= max(x) * .5, "*", "")))
    x <- lab
    x
  }) %>% t()
  col_tiss <- colorRamp2(c(min(mat), 0, max(mat)), c("purple", "white", "darkorange"))
  P1 <- ComplexHeatmap::Heatmap(
    matrix = mat, border_gp = gpar(col = "black", lty = 1),
    col = colorRamp2(c(min(mat), 0, max(mat)), c("white", "white", "white")),
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      grid.circle(x = x, y = y, r = 0.1, #* min(unit.c(width)), #(mat_fin[i, j])/3 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_tiss(mat[i, j]), col = "white"))
      grid.text(rowLabels[i, j], x, y)
    },
    row_names_gp = grid::gpar(fontsize = fontsize),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = fontsize),
    width = ncol(mat)*unit(size, "mm"), 
    height = nrow(mat)*unit(size, "mm"),
    column_names_side = "bottom",
    column_dend_side = "bottom",
    column_title = "Tissue\nEnrichment", column_title_gp = grid::gpar(fontsize = titlefontsize, fontface = "bold"), show_row_dend = FALSE, show_column_dend = FALSE,
    show_row_names = FALSE, show_heatmap_legend = FALSE,
    name = "Tissue\nEnrichment"
  )
  
  
  P1
  return(list(P1, col_tiss))  
}




#' Produce a combination heatmap old ver
#' 
#' @param SerObj A seurat object
#' @param markerVec A vector of markers
#' @param pairedList A list of paired data
#' @param pairedList2 A second list of paired data
#' @param labelColumn Column name for sample labels
#' @param prefix Prefix for output file names
#' @param adtoutput Output type (default: "unpaired")
#' @param rowsplit Split data by row (default: NULL)
#' @param columnsplit Split data by column (default: NULL)
#' @param size Plot size
#' @param coldend Show dendrogram for columns (default: TRUE)
#' @param rowdend Show dendrogram for rows (default: TRUE)
#' @param coldendside Side to place column dendrogram (default: "bottom")
#' @param rowdendside Side to place row dendrogram (default: "left")
#' @param fontsize Font size (default: 12)
#' @param titlefontsize Title font size (default: 20)
#' @param gap Gap between panels (default: 0)
#' @return A list of differentially expressed genes
#' @export
ProduceComboHeatmap.old <- function(SerObj, markerVec, pairedList, pairedList2=NULL, 
                                labelColumn, prefix, adtoutput = "unpaired", 
                                rowsplit = NULL, columnsplit = NULL, 
                                size, coldend = TRUE, rowdend = TRUE, 
                                coldendside = "bottom", rowdendside = "left",
                                row_names_side = "left",  column_names_side = "bottom",
                                fontsize = 12, titlefontsize = 20, gap = 0,
                                show_column_dend = T, show_row_dend = T){
  

  
  
  P1 = scCustFx::make_RNA_heatmap(SerObj = SerObj, markerVec = markerVec, labelColumn = labelColumn, 
                                  rowsplit = rowsplit, columnsplit=columnsplit, size=size,
                                  clus_cols = rowdend, show_column_dend = show_column_dend,
                                  clus_rows=coldend, show_row_dend = show_row_dend,
                                  # coldendside=coldendside, rowdendside=rowdendside, 
                                  fontsize = fontsize, titlefontsize = titlefontsize,  pairedList2 = pairedList2,
                                  column_names_side = column_names_side,
                                  column_dend_side = coldendside,
                                  row_names_side =row_names_side,
                                  row_dend_side = rowdendside)
  col_RNA = P1$col_RNA
  P1 = P1$plot
  
  if (adtoutput == "unpaired"){
    P2 <- scCustFx::CreateProteinHeatmap.old(SerObj, proteinList = unname(unlist(pairedList)), labelColumn=labelColumn, prefix = prefix, size, fontsize = fontsize, titlefontsize = titlefontsize)
    col_ADT <- P2[[2]]
    P2 <- P2[[1]]
  } else {
    P2 <- scCustFx::CreatePairedHeatmap.old(SerObj, pairedList, labelColumn=labelColumn, prefix = prefix)
    
  }
  
  
  P3 <- scCustFx::plotEnrichment.old(SerObj, field1 = labelColumn, field2 = 'Tissue', 
                                 size, fontsize = fontsize, titlefontsize = titlefontsize,
                                 column_title = "Tissue\nEnrichment")
  col_tiss <- P3[[2]]
  P3 <- P3[[1]]
  
  plotlist <- P1 + P2 + P3
  # draw(plotlist, ht_gap = unit(0, "mm"))
  lgd1 <- Legend(title = "Scaled Avg. Exp.", col_fun = col_RNA, direction = "horizontal", title_gp = gpar(fontsize = 14),
                 labels_gp = gpar(fontsize = 12))
  lgd2 <- Legend(title = "Scaled Avg. ADT", col_fun = col_ADT, direction = "horizontal", title_gp = gpar(fontsize = 14),
                 labels_gp = gpar(fontsize = 12))
  lgd3 <- Legend(title = "Tissue Enrichment", col_fun = col_tiss, direction = "horizontal", title_gp = gpar(fontsize = 14),
                 labels_gp = gpar(fontsize = 12))
  pd = packLegend(lgd1, lgd2, lgd3, direction = "vertical")
  draw(plotlist, heatmap_legend_list = pd, ht_gap = unit(gap,"mm"), heatmap_legend_side = "right")
}

#' Create a protein heatmap old ver
#' 
#' @param SerObj A Seurat object.
#' @param proteinList A list of proteins to plot.
#' @param labelColumn Column name to use for labeling rows.
#' @param prefix Prefix for plot title.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @return A heatmap of protein expression.
#' @export
CreateProteinHeatmap.old <- function(SerObj, proteinList, labelColumn, 
                                 prefix, size, fontsize, titlefontsize) {
  print(proteinList)
  adt_columns <- proteinList # unname(unlist(pairedList))
  mat_ADT <- cbind("Label" = SerObj@meta.data[[labelColumn]], SerObj@meta.data[,adt_columns])
  colnames(mat_ADT) <- gsub(paste0(prefix, "."), "", colnames(mat_ADT))
  colnames(mat_ADT) <- gsub("_UCell_pos", "", colnames(mat_ADT)) %>% as.factor()
  colnames(mat_ADT) <- gsub("\\.", "-", colnames(mat_ADT)) %>% as.factor()
  ordering <- names(mat_ADT)
  mat_ADT <- mat_ADT %>% pivot_longer(cols = 2:length(mat_ADT), names_to = "ADT") %>% group_by(Label, ADT) %>% summarize(avg = mean(value)) %>% pivot_wider(id_cols = Label, names_from = "ADT", values_from = "avg") %>% as.data.frame()
  mat_ADT <- mat_ADT[, ordering]
  rownames(mat_ADT) <- mat_ADT$Label
  mat_ADT <- mat_ADT %>% select(-Label)
  ordering <- names(mat_ADT)
  
  avgSeuratADT <- Seurat::AverageExpression(SerObj, group.by = labelColumn,
                                            layer = 'data', return.seurat = T,
                                            assays = prefix)
  avgSeuratADT <- Seurat::NormalizeData(avgSeuratADT)
  mat <- t(as.matrix(avgSeuratADT@assays[[prefix]]@data)) %>% pheatmap:::scale_mat(scale = 'column') %>% as.data.frame()
  mat <- mat[,colnames(mat) %in% colnames(mat_ADT)]
  mat <- mat[,ordering]
  
  col_ADT = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("white", "gray85", "blue"))
  
  P2 <- Heatmap(mat_ADT, name = "Scaled\nAvg. ADT", col = col_ADT, rect_gp = gpar(type="none"), border_gp = gpar(col = "black", lty = 1), show_row_dend = FALSE, show_column_dend = FALSE, width = ncol(mat)*unit(size, "mm"), 
                height = nrow(mat)*unit(size, "mm"), column_names_gp = grid::gpar(fontsize = fontsize), row_names_gp = grid::gpar(fontsize = fontsize),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.circle(x = x, y = y, r = mat_ADT[i,j]/3* min(unit.c(width)), #(mat_fin[i, j])/3 * min(unit.c(width, height)), 
                              gp = gpar(fill = col_ADT(mat[i, j]), col = NA))
                }, cluster_rows = TRUE, cluster_columns = FALSE,
                column_names_side = "bottom",
                column_dend_side = NULL,
                column_names_rot = 45,
                row_dend_side = NULL,
                column_title_gp = grid::gpar(fontsize = titlefontsize), column_title = "Surface\nProtein",
                show_row_names = FALSE, show_column_names = TRUE, show_heatmap_legend = FALSE, row_split = NULL)
  return(list(P2, col_ADT))
}


#' Plot Enrichment 
#' 
#' @param SerObj A Seurat object.
#' @param field1 Field 1 to compare.
#' @param field2 Field 2 to compare.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param column_title Column title for plot.
#' @return A plot of enrichment analysis between two fields.
#' @export
plotEnrichment.old <- function(SerObj, field1, field2, size, fontsize, titlefontsize, column_title) {
  mat <- asinh(chisq.test(table(SerObj[[field1]][[field1]], SerObj[[field2]][[field2]]))$res) %>% as.matrix()  
  
  rowLabels <- apply(mat, 1, function(x){
    # ToDo:
    lab = ifelse(x >= 0.9 * max(x), "***", ifelse(x >= 0.75 * max(x), "**", ifelse(x >= max(x) * .5, "*", "")))
    x <- lab
    x
  }) %>% t()
  col_tiss <- circlize::colorRamp2(c(min(mat), 0, max(mat)), c("purple", "white", "darkorange"))
  P1 <- ComplexHeatmap::Heatmap(
    matrix = mat, border_gp = gpar(col = "black", lty = 1),
    col = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("white", "white", "white")),
    cell_fun = function(j, i, x, y, width, height, fill) {
      
      grid.circle(x = x, y = y, r = 0.5* min(unit.c(width)), #(mat_fin[i, j])/3 * min(unit.c(width, height)), 
                  gp = gpar(fill = col_tiss(mat[i, j]), col = "white"))
      grid.text(rowLabels[i, j], x, y)
    },
    row_names_gp = grid::gpar(fontsize = fontsize),
    column_names_rot = 45,
    column_names_gp = grid::gpar(fontsize = fontsize),
    width = ncol(mat)*unit(size, "mm"), 
    height = nrow(mat)*unit(size, "mm"),
    column_names_side = "bottom",
    column_dend_side = "bottom",
    column_title = column_title, 
    column_title_gp = grid::gpar(fontsize = titlefontsize), 
    show_row_dend = FALSE, show_column_dend = FALSE,
    show_row_names = FALSE, show_heatmap_legend = FALSE,
    name = column_title
  )
  
  
  P1
  return(list(P1, col_tiss))  
}




#' Subset a Seurat object by coordinates and find DE genes
#'
#' This function subsets a Seurat object by the coordinates of a dimensionality reduction, such as tSNE, UMAP, or PCA, and finds differentially expressed genes in the subsetted object.
#' 
#' @param SerObj A Seurat object to subset
#' @param reduction_method Character input for the dimensionality reduction method to subset by. Accepts "tSNE", "UMAP", or "PCA". Default is "tSNE"
#' @param x_range A numeric vector of length 2 for the x-coordinate range to subset by
#' @param y_range A numeric vector of length 2 for the y-coordinate range to subset by
#' @param logfc.threshold A numeric threshold for the log fold change of genes to be considered differentially expressed. Default is 0.25
#' @param test.use A character input for the test to use for finding differentially expressed genes. Default is "wilcox"
#' @param layer Character input for which layer to use for finding differentially expressed genes. Default is "data"
#' @param min.pct A numeric threshold for the minimum percentage of cells expressing a gene to be considered. Default is 0.65
#' @param min.diff.pct A numeric threshold for the minimum difference in percentage of cells expressing a gene between groups to be considered differentially expressed. Default is 0.2
#' @param only.pos Logical input indicating whether to only consider positively expressed genes. Default is TRUE
#' @param savePathName A character input for the path to save the differentially expressed genes as an RDS file. Default is NULL
#' @return A list of differentially expressed genes
#' @export
#' @examples
#' DEgenes_unsupclusts <- SubsetSerAndFindDEGenes(SerObj = pbmc, reduction_method = "tSNE", x_range = c(-10, 10), y_range = c(-20, 20), savePathName = "DEgenes_unsupclusts.rds")
FindAllMarkers_SubsetCoord <- function(SerObj, reduction_method = "tSNE", x_range = c(-1, 1), y_range = c(-1, 1),
                                       logfc.threshold = 0.25, test.use = "wilcox", layer = "data",
                                       min.pct = 0.65, min.diff.pct = 0.2, only.pos = T,
                                       savePathName = NULL){
  SerObj_sub = SubsetSerByCoordinates(SerObj=SerObj, reduction_method = reduction_method, x_range = x_range, y_range = y_range)
  
  DEgenes_unsupclusts = FindAllMarkers(SerObj,
                            logfc.threshold = logfc.threshold,
                            test.use = test.use, layer = layer,
                            min.pct = min.pct, min.diff.pct = min.diff.pct,
                            only.pos = only.pos)

  if(!is.null(savePathName)){
    saveRDS(DEgenes_unsupclusts, savePathName)
  }
  
  return(DEgenes_unsupclusts)
}


#' @title SubsetSerByCoordinates
#'
#' @description Remove unwanted HTOS from Seurat objects
#' @param SerObj A Seurat Object 
#' @param reduction_method an input of  = c("tSNE", "UMAP", "PCA") default tSNE
#' @param x_range an input of  two numbers defining the x range default c(-1, 1)
#' @param y_range aan input of  two numbers defining the y range default c(-1, 1)
#' @return A subset Seurat object.
#' @export
SubsetSerByCoordinates <- function(SerObj, reduction_method = "tSNE", assay = "RNA",
                                   x_range = c(-1, 1), y_range = c(-1, 1)){

  #TODO: deal with assay
  if(reduction_method == "tSNE"){
    dim_1 <- "rnaTSNE_1"
    dim_2 <- "rnaTSNE_2"
  }else if(reduction_method == "UMAP"){
    dim_1 <- "rnaUMAP_1"
    dim_2 <- "rnaUMAP_1"
  }else if(reduction_method == "PCA"){
    #TODO deal with additional PCA dims
    dim_1 <- "PC_1"
    dim_2 <- "PC_2"
  }
  
  cells = (SerObj@reductions[[reduction_method]]@cell.embeddings[,dim_1] > x_range[1] & 
             SerObj@reductions[[reduction_method]]@cell.embeddings[,dim_1] < x_range[2]) & 
    (SerObj@reductions[[reduction_method]]@cell.embeddings[,dim_2] > y_range[1] & 
       SerObj@reductions[[reduction_method]]@cell.embeddings[,dim_2] < y_range[2])
  
  print(paste0("Selection included: ", length(cells), " cells"))
  subsetted_obj <- SerObj[, cells]
  
  return(subsetted_obj)
}



#' @title ReduceSerObj_HTO
#'
#' @description Remove unwanted HTOS from Seurat objects
#' @param SerObj A Seurat Object 
#' @param removeHTOs Vector of names of HTOs to remove, if NULL, removeHTOs = c("Discordant", "Doublet", "ND")
#' @return A subset Seurat object.
#' @export
ReduceSerObj_HTO <- function(SerObj, removeHTOs = NULL){
  if(is.null(removeHTOs)){
    print("removeHTOs is null, therefore, deep cleaning")
    removeHTOs = c("Discordant", "Doublet", "ND")
  }
  if (sum(removeHTOs %in% names(table(SerObj$HTO)))==0){
    stop("removeHTOs defined not in seurat object ")
  } else {
    keepHTOs <- setdiff(names(table(SerObj$HTO)), removeHTOs)
  }
  
  print("Keept HTOs: ")
  print(keepHTOs)
  return(subset(SerObj, HTO %in% keepHTOs))
}

#' @title DownSample_SerObj_perLevel
#'
#' @description downsample per feature levels
#' @param SerObj A Seurat Object 
#' @param featureSplit a feature of seurat object to split. "BarcodePrefix" # this is a good way to sample per seq run if possible
#' @param totalCell  0 # is the default way to downsample by min common cells
#' @return A subseted downsampled Seurat object.
#' @export
DownSample_SerObj_perLevel <- function(SerObj, featureSplit=NULL, totalCell = 300){

  availFeats = colnames(SerObj@meta.data)
  
  if(!is.null(featureSplit)){
    
    if(length(featureSplit)>1) {
      warning("length of featureSplit > 1 only 1st is used")
      featureSplit = featureSplit[1]
      
      if(!(featureSplit %in% availFeats)) stop("featureSplit not found in meta.data")
      
    }
    
  } else stop("featureSplit is NULL")
  
  if(totalCell == 300) warning("totalCell = 300 is the default")
  
  CellCountMat = table(SerObj@meta.data[,featureSplit]) %>% unlist() %>% as.matrix()
  
  if(totalCell == 0) {
    totalCell = min(CellCountMat[,1])
    print(paste0("total cells set by min to:", totalCell))
  }
  
  if(!("barcode" %in% availFeats)) {
    warning("barcode not found in meta.data")
    SerObj$barcode = paste("cell_", 1:nrow(SerObj@meta.data))
  }
  
  
  splitLevs = levels(factor(SerObj@meta.data[,featureSplit]))
  
  barcodeLS = lapply(splitLevs, function(xSL){
    availBarcodes = SerObj@meta.data[which(SerObj@meta.data[,featureSplit] == xSL),]$barcode
    if(length(availBarcodes)<totalCell) {
      warning(paste0(xSL, " had less cells than totalCell"))
      availBarcodes
    }else {
      sample(availBarcodes, totalCell, replace = F)
    }
    
    
  })
  
  return(subset(SerObj, barcode %in% unlist(barcodeLS)))
  
}



#' @title compressSerHD5
#'
#' @description compressSerHD5 will save space on your drives.
#' @param SerObj A Seurat Object 
#' @param load.path the path to a seurat .rds file
#' @param save.path the path to the desired hd5 seurat file. Default is NULL which then load.path is used to save the new object next to it.
#' @param overwrite  default F overwirtes the new hd5 if exists
#' @param updateSerObj  default F if T runs UpdateSerObj()
#' @param returnSerObj  default F if T returns the Seurat object to be used in a pipeline
#' @return Saves HD5 Seurat object to path given
#' @export
compressSerHD5 <- function(SerObj = NULL, load.path = NULL, overwrite = F,
                           save.path = NULL, updateSerObj = F, returnSerObj = F){
  if(is.null(SerObj) & is.null(load.path)) stop("give seurat object SerObj or load.path to one")
  
  if((!is.null(SerObj)) & (!is.null(load.path))) warning("both SerObj and load.path given, load.path is ignored")
  
  
  if(is.null(SerObj)){
    print("loading in RDS")
    SerObj = readRDS(load.path)
    print("load of RDS from drives complete")
  }
  
  if(updateSerObj) SerObj = SeuratDisk::UpdateSerObj(SerObj)
  
  if(is.null(save.path)) save.path = gsub(".rds", "hd5Ser", load.path)
  print("saving file path:")
  print(save.path)
  SeuratDisk::SaveH5Seurat(SerObj, filename=save.path,  overwrite = overwrite)
  print("save complete")
  if(returnSerObj) return(SerObj)
  
}

#' @title scatter_plot_seurat
#'
#' @description create scatter plot
#' @param seurat_obj A Seurat Object 
#' @return ggplot
#' @export
scatter_plot_seurat <- function(seurat_obj, meta_feature, gene_name, x_split_value, y_split_value, 
                                assay = "RNA", layer="counts", base_size = 11) {
  require(broom)
  DefaultAssay(seurat_obj) = assay
  
  x_cell_ids <- which(seurat_obj@meta.data[,meta_feature] == x_split_value)
  y_cell_ids <- which(seurat_obj@meta.data[,meta_feature] == y_split_value)
  
  x_expr <- GetAssayData(seurat_obj, layer = layer)[gene_name, x_cell_ids]
  y_expr <- GetAssayData(seurat_obj, layer = layer)[gene_name, y_cell_ids]
  
  downSampMin = min(c(length(x_expr), length(y_expr)))
  
  x_expr = SerObjrt(sample(x_expr, downSampMin, replace = F), decreasing = F)
  y_expr = SerObjrt(sample(y_expr, downSampMin, replace = F), decreasing = F)
  
  
  r_squared <- cor(x=as.numeric(x_expr), 
                   y=as.numeric(y_expr))^2
  if(is.na(r_squared)) r_squared = 0 
  
  # linear_model <- lm(as.numeric(y_expr) ~ as.numeric(x_expr))
  # 
  # if(nrow(summary(linear_model)$coefficients)>1){
  #   p_value <- summary(linear_model)$coefficients[2,4]
  # } else { 
  #   p_value = 1
  # }
  # 
  
  
  
  ggplot(data.frame(x = x_expr, y = y_expr), aes(x = x, 
                                                 y = y )) +
    geom_jitter() +
    geom_point(aes(col="red")) +
    # stat_smooth(method = "lm", formula = y ~ x,geom = "smooth") +
    geom_smooth(method = "lm", se = FALSE, col="dodgerblue", lty=3) +
    
    xlab(paste0(meta_feature, " = " , x_split_value))+ 
    ylab(paste0(meta_feature, " = " , y_split_value)) + 
    ggtitle(paste(gene_name, #ifelse(p_value <= 0.01, " : p<=0.01", " : NS"),
                  "\nR^2 =", round(r_squared,2)#, " p=", round(p_value,5)
    )) +
    theme_bw(base_size = base_size) +
    theme(legend.position="none",
          legend.title = element_blank())
  
}