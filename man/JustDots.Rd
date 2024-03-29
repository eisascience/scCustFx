% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DotPlotFxs.R
\name{JustDots}
\alias{JustDots}
\title{Custom Dot Plot Visualization for Seurat Objects}
\usage{
JustDots(
  SerObj,
  mat,
  markerVec,
  pairedList2,
  labelColumn,
  rowsplit = NULL,
  columnsplit = NULL,
  size,
  coldend = TRUE,
  rowdend = TRUE,
  coldendside = "bottom",
  rowdendside = "left",
  fontsize = 12,
  titlefontsize = 14,
  gap = 0,
  assay = "RNA"
)
}
\arguments{
\item{SerObj}{A Seurat object containing the data to be visualized.}

\item{mat}{A matrix of data to be visualized, typically derived from a Seurat object.}

\item{markerVec}{A vector of marker genes to be included in the plot.}

\item{pairedList2}{A named list where each name corresponds to a marker gene in `markerVec`
and each value is the name to replace it with in the plot.}

\item{labelColumn}{The name of the metadata column in `SerObj` to use for labeling.}

\item{rowsplit}{Optional vector specifying how to split rows in the plot; NULL means no splitting.}

\item{columnsplit}{Optional vector specifying how to split columns in the plot; NULL means no splitting.}

\item{size}{Size of the dots in the plot.}

\item{coldend}{Logical indicating whether to include a column dendrogram.}

\item{rowdend}{Logical indicating whether to include a row dendrogram.}

\item{coldendside}{Specifies the side on which to draw the column dendrogram ("top" or "bottom").}

\item{rowdendside}{Specifies the side on which to draw the row dendrogram ("left" or "right").}

\item{fontsize}{Font size for text in the plot.}

\item{titlefontsize}{Font size for the plot title.}

\item{gap}{Gap between dots in the plot.}
}
\value{
A ComplexHeatmap object representing the custom dot plot. This object can be directly
plotted in R or further customized using the ComplexHeatmap and circlize packages.
}
\description{
This function creates a custom dot plot visualization for Seurat objects. It is designed
to allow for complex customizations such as splitting rows and columns, adjusting dot sizes,
and incorporating dendrograms. It scales the expression data, applies custom naming conventions,
and utilizes a heatmap to represent data visually.
}
\examples{
JustDots(mySeuratObject, mat = expressionMatrix, 
         markerVec = c("GeneA", "GeneB"), 
         pairedList2 = list(GeneA = "Gene A", GeneB = "Gene B"), 
         labelColumn = "ident")

}
