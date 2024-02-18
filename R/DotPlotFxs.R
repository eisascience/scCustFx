#' Enhanced Dot Plot Visualization for Seurat Objects
#'
#' This function provides an enhanced dot plot visualization for Seurat objects,
#' allowing for flexible grouping, assay selection, and dendrogram inclusion options.
#' It averages expression data by specified groups and then utilizes a customized
#' plotting function to display the results.
#'
#' @param SerObj A Seurat object containing the data to be visualized.
#' @param group.by The name of the metadata column in `SerObj` to group the data by.
#' @param slot The slot from which to retrieve data. Common slots include "data", "scale.data", and "counts".
#' @param assay The name of the assay to use for data retrieval. Default is "RBA".
#' @param markerVec A vector of marker genes to include in the plot.
#' @param rowsplit Optional. A vector specifying how to split rows in the plot. NULL means no splitting.
#' @param columnsplit Optional. A vector specifying how to split columns in the plot. NULL means no splitting.
#' @param size The size of the dots in the plot.
#' @param coldend Logical, whether to include a column dendrogram.
#' @param rowdend Logical, whether to include a row dendrogram.
#' @param coldendside The side on which to draw the column dendrogram ("top" or "bottom").
#' @param rowdendside The side on which to draw the row dendrogram ("left" or "right").
#' @param fontsize The font size for text in the plot.
#' @param titlefontsize The font size for the plot title.
#' @param gap The gap between dots in the plot.
#'
#' @return Generates a dot plot visualization based on the provided Seurat object and parameters.
#'
#' @examples
#' DotPlot2(mySeuratObject, group.by = "ident", slot = "data", assay = "RNA",
#'          markerVec = c("GeneA", "GeneB", "GeneC"))
#'
#' @importFrom circlize circos.initialize circos.par circos.trackPlotRegion circos.text circos.points
#' @importFrom tidyr gather spread
#' @export
DotPlot2 <- function(SerObj, group.by, slot="data", assay="RBA", markerVec, 
                     rowsplit = NULL, columnsplit = NULL, 
                     size = 6, coldend = FALSE, rowdend = FALSE, 
                     coldendside = "bottom", rowdendside = "left", 
                     fontsize = 12, titlefontsize = 14, gap = 0) {
  
  library(circlize)
  library(tidyr)
  
  avgSeurat <- Seurat::AverageExpression(SerObj, group.by = group.by,
                                         slot = slot, return.seurat = TRUE,
                                         assays = assay)
  mat <- as.matrix(Seurat::GetAssayData(avgSeurat, slot = slot))
  
  JustDots(SerObj, mat, 
           markerVec = markerVec, 
           pairedList2 = NULL, 
           labelColumn = group.by,
           rowsplit = rowsplit, columnsplit = columnsplit, 
           size = size, coldend = coldend, rowdend = rowdend, 
           coldendside = coldendside, rowdendside = rowdendside, 
           fontsize = fontsize, titlefontsize = titlefontsize, gap = gap)
}

#' Custom Dot Plot Visualization for Seurat Objects
#'
#' This function creates a custom dot plot visualization for Seurat objects. It is designed
#' to allow for complex customizations such as splitting rows and columns, adjusting dot sizes,
#' and incorporating dendrograms. It scales the expression data, applies custom naming conventions,
#' and utilizes a heatmap to represent data visually.
#'
#' @param ComboSerObj A Seurat object containing the data to be visualized.
#' @param mat A matrix of data to be visualized, typically derived from a Seurat object.
#' @param markerVec A vector of marker genes to be included in the plot.
#' @param pairedList2 A named list where each name corresponds to a marker gene in `markerVec`
#' and each value is the name to replace it with in the plot.
#' @param labelColumn The name of the metadata column in `ComboSerObj` to use for labeling.
#' @param rowsplit Optional vector specifying how to split rows in the plot; NULL means no splitting.
#' @param columnsplit Optional vector specifying how to split columns in the plot; NULL means no splitting.
#' @param size Size of the dots in the plot.
#' @param coldend Logical indicating whether to include a column dendrogram.
#' @param rowdend Logical indicating whether to include a row dendrogram.
#' @param coldendside Specifies the side on which to draw the column dendrogram ("top" or "bottom").
#' @param rowdendside Specifies the side on which to draw the row dendrogram ("left" or "right").
#' @param fontsize Font size for text in the plot.
#' @param titlefontsize Font size for the plot title.
#' @param gap Gap between dots in the plot.
#'
#' @return A ComplexHeatmap object representing the custom dot plot. This object can be directly
#' plotted in R or further customized using the ComplexHeatmap and circlize packages.
#'
#' @examples
#' JustDots(mySeuratObject, mat = expressionMatrix, 
#'          markerVec = c("GeneA", "GeneB"), 
#'          pairedList2 = list(GeneA = "Gene A", GeneB = "Gene B"), 
#'          labelColumn = "ident")
#'
#' @importFrom circlize colorRamp2
#' @importFrom Seurat DotPlot AverageExpression GetAssayData BlueAndRed
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid grid.circle gpar unit
#' @export
JustDots <- function(ComboSerObj, mat, 
                     markerVec, 
                     #pairedList, 
                     pairedList2, 
                     labelColumn,
                     #prefix, 
                     #adtoutput = "unpaired",
                     rowsplit = NULL, columnsplit = NULL, 
                     size, coldend = TRUE, rowdend = TRUE, 
                     coldendside = "bottom", rowdendside = "left", 
                     fontsize = 12, titlefontsize = 14, gap = 0){
  
  
  mat <- mat  %>% as.data.frame() %>% filter(rownames(mat) %in% markerVec) %>% as.matrix() %>% pheatmap:::scale_mat(scale = 'row') %>% as.data.frame()
  
  
  colnames(mat) <- CellMembrane::RenameUsingCD(colnames(mat))
  
  mat <- t(as.matrix(mat)) %>% as.data.frame()
  mat1 <- mat
  
  for (nam in names(pairedList2)){
    foundnam_pos <- grep(nam, colnames(mat))
    rep <- pairedList2[[nam]]
    colnames(mat)[foundnam_pos] <- rep
  }
  
  col_RNA = circlize::colorRamp2(c(min(mat), 0, max(mat)), c(Seurat::BlueAndRed(20)[1], "gray85", Seurat::BlueAndRed(20)[20]), space = "sRGB")
  
  # col_RNA = c(Seurat::BlueAndRed(20)[c(1,3,5,7)], Seurat::BlueAndRed(20)[c(14,16,18,20)])
  #Above emulates Seurat's BlueAndRed color scheme.
  fullorder <- rownames(mat)
  fullorder1 <- rownames(mat1)
  
  plt <- Seurat::DotPlot(ComboSerObj, features = unique(markerVec), group.by = labelColumn, assay = "RNA")
  pct <- plt$data %>% select(pct.exp, id, features.plot) %>% pivot_wider(id_cols = features.plot, names_from = id, values_from = pct.exp) %>% as.data.frame()
  row.names(pct) <- pct$features.plot
  pct <- pct %>% select(-features.plot) %>% as.matrix() %>% t()
  
  pct <- pct[fullorder1,]
  colordering <- colnames(mat1)
  pct <- pct[,colordering]
  
  P1 <-
    ComplexHeatmap::Heatmap(mat,
                            width = ncol(mat)*unit(size, "mm"), 
                            
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              grid.circle(x = x, y = y, r = sqrt(pct[i,j])* min(unit.c(width)/25), #(mat_fin[i, j])/3 * min(unit.c(width, height)), 
                                          gp = gpar(fill = col_RNA(mat[i, j]), col = NA))
                            },
                            rect_gp = gpar(type="none"),
                            border_gp = gpar(col = "black", lty = 1), 
                            height = nrow(mat)*unit(size, "mm"),
                            row_names_side = "left",
                            row_dend_side = rowdendside,
                            column_names_rot = 45,
                            # col = col_RNA,
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
  return(P1)
}

#' Calculate Average Expression Matrix for Marker Genes in Seurat Object
#'
#' Computes the average expression for a specified set of marker genes across groups
#' defined by a label column within a Seurat object. This function facilitates the
#' exploration of gene expression patterns across different cell types or conditions
#' defined in the dataset.
#'
#' @param ComboSerObj A Seurat object containing single-cell RNA-seq data.
#' @param labelColumn A string specifying the name of the column in the metadata of
#' the Seurat object that defines the groups (e.g., cell types or conditions) for
#' which the average expression should be calculated.
#' @param markerVec A vector of strings representing the marker genes for which
#' the average expression is to be computed.
#'
#' @return A matrix with rows corresponding to marker genes and columns to the
#' groups defined in `labelColumn`. Each cell in the matrix represents the average
#' expression of a gene in a group.
#'
#' @examples
#' avgExprMatrix <- GetAvgExpressionMatrix(mySeuratObject, "cell_type", c("GeneA", "GeneB", "GeneC"))
#'
#' @importFrom Seurat AverageExpression GetAssayData
#' @export
GetAvgExpressionMatrix <- function(ComboSerObj, labelColumn, markerVec) {
  avgSeurat <- Seurat::AverageExpression(ComboSerObj, group.by = labelColumn,
                                         features = markerVec,
                                         slot = 'data', return.seurat = T,
                                         assays = 'RNA')
  # avgSeurat <- Seurat::NormalizeData(avgSeurat)
  mat <- t(as.matrix(Seurat::GetAssayData(avgSeurat, slot = 'data')))
  return(mat)
}
