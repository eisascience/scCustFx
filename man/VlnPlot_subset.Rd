% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations_Seurat_based.R
\name{VlnPlot_subset}
\alias{VlnPlot_subset}
\title{Calculate the percent expressed cells for each group in a Seurat object}
\usage{
VlnPlot_subset(
  SerObj,
  group.by,
  features,
  plot = F,
  topN = 10,
  xlab = "",
  cols = c("#8DD3C7", "#33A02C", "#F4CAE4", "#A6CEE3", "#FDCDAC", "#B2DF8A", "#E78AC3",
    "#1F78B4", "#FB9A99", "#E31A1C", "#66C2A5", "#FF7F00"),
  addDimplot = T,
  returnLs = T
)
}
\arguments{
\item{SerObj}{A Seurat object}

\item{group.by}{A meta feature used to group the cells}

\item{features}{A vector of features (genes) to calculate the percent expressed cells for}

\item{plot}{A logical indicating whether to plot}

\item{topN}{select top N genes}

\item{cols}{color vector}
}
\value{
A data frame with the percent expressed cells for each group and feature
}
\description{
Calculate the percent expressed cells for each group in a Seurat object
}
