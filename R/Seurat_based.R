#' scatter-bubble style of selected PCs and genes
#'
#'scatter-bubble style of selected PCs and genes
#' 
#' @param serObj a seurat obj 
#' @return A Processed Seurat Obj
#' @export
process_SCTbase = function(serObj = NULL, ncells = 3000, assay = "Spatial",
                           dims=1:30, verbose = T, clusterResolutions = c(0.2, 0.4, 0.6, 0.8, 1.2)){
  print("SCT transformation start")
  serObj <- SCTransform(serObj, assay = assay, 
                        ncells = ncells, # of cells used to build NB regression def is 5000
                        verbose = verbose)
  print("Running PCA start")
  serObj <- RunPCA(serObj)
  serObj <- RunTSNE(serObj, dims = dims, perplexity = CellMembrane:::.InferPerplexityFromSeuratObj(serObj, perplexity = 30),  check_duplicates = FALSE)
  print("Running UMAP start")
  serObj <- RunUMAP(serObj, dims = dims)
  print("Running Clustering start")
  serObj <- FindNeighbors(serObj, dims = dims)
  for (resolution in clusterResolutions) {
    serObj <- FindClusters(object = serObj, resolution = resolution, verbose = verbose)
    serObj[[paste0("ClusterNames_", resolution)]] <- Idents(object = serObj)
  }
  return(serObj)
}



#' scatter-bubble style of selected PCs and genes
#'
#'scatter-bubble style of selected PCs and genes
#' 
#' @param serObj a seurat obj 
#' @param features features to show
#' @param pcs a numeric vector of 2 pc nnumbers e.g. c(1, 2) 
#' @param base_size base_size of plot
#' @param quaT color pink via quantiles
#' @return ggplot scatter-bubble plot style 
#' @export
plotPCgeneLoadings = function(serObj, features, pcs = c(1, 2), quaT = 0.8, base_size =10){
  tempDF = as.data.frame(serObj@reductions$pca@feature.loadings[features, pcs])
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
#' @param serobj A Seurat object
#' @param group.by A meta feature used to group the cells
#' @param features A vector of features (genes) to calculate the percent expressed cells for
#' @param plot A logical indicating whether to plot 
#' @param topN select top N genes  
#' @param cols color vector
#' @return A data frame with the percent expressed cells for each group and feature
#' @export
VlnPlot_subset <- function(serobj, group.by, features, plot=F, topN = 10, xlab = "",
                           cols=c("#8DD3C7", "#33A02C", "#F4CAE4", "#A6CEE3", "#FDCDAC",
                                  "#B2DF8A","#E78AC3", "#1F78B4", "#FB9A99", "#E31A1C", 
                                  "#66C2A5", "#FF7F00"), addDimplot = T, returnLs = T){
  
  
  #fix factor
  serobj@meta.data[,group.by] = naturalsort::naturalfactor(as.character(serobj@meta.data[,group.by]))
  
  features = unique(features)

  
  PctExprDF =  scCustFx:::PercentExpressed(serobj, 
                                           group.by=group.by, 
                                           features=features, 
                                           plot_perHM = plot)
  
  PctExprDF$max = apply(PctExprDF, 1, max)
  PctExprDF$min = apply(PctExprDF, 1, min)
  
  PctExprDF$maxminDelta = PctExprDF$max-PctExprDF$min
  
  
  Idents(serobj) = group.by
  
  vln <- VlnPlot(serobj, 
          features = rownames( PctExprDF[order(PctExprDF$maxminDelta, decreasing = T)[1:topN],]), 
          stack = TRUE, flip = TRUE, fill.by = "ident",
          cols = cols) + NoLegend() + xlab(xlab) #+
    # theme_classic(base_size = 14) 
  
  
  if(addDimplot){
    dim <- DimPlot(serobj, 
            label = T, label.size = 10, label.box = T, repel = T,
            cols=cols,
            group.by=group.by) + 
      coord_flip()  + scale_y_reverse() +
      theme_classic(base_size = 14) + 
      NoLegend() +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank()
      )
    if(returnLs) {
      list(dim=dim, vln=vln)
    } else {
      dim / vln
    }
  } else  {
    vln
  }
 
  
}





#' Calculate the percent expressed cells for each group in a Seurat object
#'
#' @param so A Seurat object
#' @param group.by A meta feature used to group the cells
#' @param features A vector of features (genes) to calculate the percent expressed cells for
#' @param plot_perHM A logical indicating whether to plot a heatmap of the results
#' @return A data frame with the percent expressed cells for each group and feature
#' @export
PercentExpressed <- function(so, group.by, features, plot_perHM = F){
  
  print("removing:")
  print(features[!features %in% rownames(so)])
  
  
  features = features[features %in% rownames(so)]
  MyDot = DotPlot(so, group.by = group.by, features = (features) )
  
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
#' @param seuratObj A Seurat object.
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
make_RNA_heatmap = function(seuratObj, markerVec, labelColumn, rowsplit=NULL, columnsplit=NULL,
                            size=6,
                            clus_cols = TRUE, show_column_dend = T,
                            clus_rows=TRUE, show_row_dend = T,
                            fontsize=10, titlefontsize=10, pairedList2=NULL, asGG = F,
                            column_names_side = "bottom",
                            column_dend_side = "bottom",
                            row_names_side = "left",
                            row_dend_side = "left",
                            assays = 'RNA', useCDnames = T){
  avgSeurat <- Seurat::AverageExpression(seuratObj, group.by = labelColumn,
                                         features = markerVec,
                                         slot = 'counts', return.seurat = T,
                                         assays = assays)
  avgSeurat <- Seurat::NormalizeData(avgSeurat)
  mat <- t(as.matrix(Seurat::GetAssayData(avgSeurat, slot = 'data')))
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
#' @param seuratObj A Seurat object.
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
make_RNA_heatmap2 = function(seuratObj, labelColumn = 'ClusterNames_0.2', 
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
  
  
  mat <- asinh(scale(t(AverageExpression(seuratObj, group.by = labelColumn, features = markerVec)[[assay]])))
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

#' Produce a combination heatmap
#' 
#' @param seuratObj A seurat object
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
ProduceComboHeatmap <- function(seuratObj, markerVec, pairedList, pairedList2=NULL, 
                                labelColumn, prefix, adtoutput = "unpaired", 
                                rowsplit = NULL, columnsplit = NULL, 
                                size, coldend = TRUE, rowdend = TRUE, 
                                coldendside = "bottom", rowdendside = "left",
                                row_names_side = "left",  column_names_side = "bottom",
                                fontsize = 12, titlefontsize = 20, gap = 0,
                                show_column_dend = T, show_row_dend = T){
  

  
  
  P1 = scCustFx::make_RNA_heatmap(seuratObj = seuratObj, markerVec = markerVec, labelColumn = labelColumn, 
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
    P2 <- scCustFx::CreateProteinHeatmap(seuratObj, proteinList = unname(unlist(pairedList)), labelColumn=labelColumn, prefix = prefix, size, fontsize = fontsize, titlefontsize = titlefontsize)
    col_ADT <- P2[[2]]
    P2 <- P2[[1]]
  } else {
    P2 <- scCustFx::CreatePairedHeatmap(seuratObj, pairedList, labelColumn=labelColumn, prefix = prefix)
    
  }
  
  
  P3 <- scCustFx::plotEnrichment(seuratObj, field1 = labelColumn, field2 = 'Tissue', 
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

#' Create a protein heatmap
#' 
#' @param seuratObj A Seurat object.
#' @param proteinList A list of proteins to plot.
#' @param labelColumn Column name to use for labeling rows.
#' @param prefix Prefix for plot title.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @return A heatmap of protein expression.
#' @export
CreateProteinHeatmap <- function(seuratObj, proteinList, labelColumn, 
                                 prefix, size, fontsize, titlefontsize) {
  print(proteinList)
  adt_columns <- proteinList # unname(unlist(pairedList))
  mat_ADT <- cbind("Label" = seuratObj@meta.data[[labelColumn]], seuratObj@meta.data[,adt_columns])
  colnames(mat_ADT) <- gsub(paste0(prefix, "."), "", colnames(mat_ADT))
  colnames(mat_ADT) <- gsub("_UCell_pos", "", colnames(mat_ADT)) %>% as.factor()
  colnames(mat_ADT) <- gsub("\\.", "-", colnames(mat_ADT)) %>% as.factor()
  ordering <- names(mat_ADT)
  mat_ADT <- mat_ADT %>% pivot_longer(cols = 2:length(mat_ADT), names_to = "ADT") %>% group_by(Label, ADT) %>% summarize(avg = mean(value)) %>% pivot_wider(id_cols = Label, names_from = "ADT", values_from = "avg") %>% as.data.frame()
  mat_ADT <- mat_ADT[, ordering]
  rownames(mat_ADT) <- mat_ADT$Label
  mat_ADT <- mat_ADT %>% select(-Label)
  ordering <- names(mat_ADT)
  
  avgSeuratADT <- Seurat::AverageExpression(seuratObj, group.by = labelColumn,
                                            slot = 'data', return.seurat = T,
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
#' @param seuratObj A Seurat object.
#' @param field1 Field 1 to compare.
#' @param field2 Field 2 to compare.
#' @param size Size of plot.
#' @param fontsize Font size for labels.
#' @param titlefontsize Font size for plot title.
#' @param column_title Column title for plot.
#' @return A plot of enrichment analysis between two fields.
#' @export
plotEnrichment <- function(seuratObj, field1, field2, size, fontsize, titlefontsize, column_title) {
  mat <- asinh(chisq.test(table(seuratObj[[field1]][[field1]], seuratObj[[field2]][[field2]]))$res) %>% as.matrix()  
  
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
#' @param so A Seurat object to subset
#' @param reduction_method Character input for the dimensionality reduction method to subset by. Accepts "tSNE", "UMAP", or "PCA". Default is "tSNE"
#' @param x_range A numeric vector of length 2 for the x-coordinate range to subset by
#' @param y_range A numeric vector of length 2 for the y-coordinate range to subset by
#' @param logfc.threshold A numeric threshold for the log fold change of genes to be considered differentially expressed. Default is 0.25
#' @param test.use A character input for the test to use for finding differentially expressed genes. Default is "wilcox"
#' @param slot Character input for which slot to use for finding differentially expressed genes. Default is "data"
#' @param min.pct A numeric threshold for the minimum percentage of cells expressing a gene to be considered. Default is 0.65
#' @param min.diff.pct A numeric threshold for the minimum difference in percentage of cells expressing a gene between groups to be considered differentially expressed. Default is 0.2
#' @param only.pos Logical input indicating whether to only consider positively expressed genes. Default is TRUE
#' @param savePathName A character input for the path to save the differentially expressed genes as an RDS file. Default is NULL
#' @return A list of differentially expressed genes
#' @export
#' @examples
#' DEgenes_unsupclusts <- SubsetSerAndFindDEGenes(so = pbmc, reduction_method = "tSNE", x_range = c(-10, 10), y_range = c(-20, 20), savePathName = "DEgenes_unsupclusts.rds")
FindAllMarkers_SubsetCoord <- function(so, reduction_method = "tSNE", x_range = c(-1, 1), y_range = c(-1, 1),
                                       logfc.threshold = 0.25, test.use = "wilcox", slot = "data",
                                       min.pct = 0.65, min.diff.pct = 0.2, only.pos = T,
                                       savePathName = NULL){
  so_sub = SubsetSerByCoordinates(so=so, reduction_method = reduction_method, x_range = x_range, y_range = y_range)
  
  DEgenes_unsupclusts = FindAllMarkers(ComboSerObj,
                            logfc.threshold = logfc.threshold,
                            test.use = test.use, slot = slot,
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
#' @param so A Seurat Object 
#' @param reduction_method an input of  = c("tSNE", "UMAP", "PCA") default tSNE
#' @param x_range an input of  two numbers defining the x range default c(-1, 1)
#' @param y_range aan input of  two numbers defining the y range default c(-1, 1)
#' @return A subset Seurat object.
#' @export
SubsetSerByCoordinates <- function(so, reduction_method = "tSNE", assay = "RNA",
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
  
  cells = (so@reductions[[reduction_method]]@cell.embeddings[,dim_1] > x_range[1] & 
             so@reductions[[reduction_method]]@cell.embeddings[,dim_1] < x_range[2]) & 
    (so@reductions[[reduction_method]]@cell.embeddings[,dim_2] > y_range[1] & 
       so@reductions[[reduction_method]]@cell.embeddings[,dim_2] < y_range[2])
  
  print(paste0("Selection included: ", length(cells), " cells"))
  subsetted_obj <- so[, cells]
  
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