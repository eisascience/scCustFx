% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_based.R
\name{BiSplit_DE}
\alias{BiSplit_DE}
\title{BiSplit_DE}
\usage{
BiSplit_DE(
  SerObj,
  gate_feat = NULL,
  gate_thr = 0,
  plot_DimPlot = F,
  doDE = T,
  cols = c("gray", "red"),
  raster = F,
  assay = "RNA",
  layer = "data",
  only.pos = F,
  min.pct = 0.65,
  min.diff.pct = 0.2,
  logfc.threshold = 0.25
)
}
\arguments{
\item{SerObj}{An object of class Seurat}

\item{gate_feat}{A character vector of feature names to use as gatekeepers for the BiSplit algorithm.}

\item{gate_thr}{The threshold value for the gatekeepers. Default is 0.}

\item{plot_DimPlot}{Logical value indicating whether to generate a DimPlot of the results. Default is FALSE.}

\item{doDE}{set doDE to F from T (default) when just visualizing the gating dimplot}

\item{cols}{A character vector of colors to use for plotting the results. Default is c("gray", "red").}

\item{raster}{Logical value indicating whether to rasterize the plot. Default is FALSE.}

\item{assay}{The assay type to use for the analysis. Default is "RNA".}

\item{layer}{The layer to use for the analysis. Default is "data".}

\item{only.pos}{Logical value indicating whether to only consider positive differential expression. Default is FALSE.}

\item{min.pct}{The minimum percentage of cells expressing a feature to consider it. Default is 0.65.}

\item{min.diff.pct}{The minimum percentage difference in expression between groups to consider a feature differentially expressed. Default is 0.2.}

\item{logfc.threshold}{The log-fold change threshold to use for identifying differentially expressed features. Default is 0.25.}
}
\value{
A DF containing the results of the differential expression analysis.
}
\description{
A function to perform differential expression analysis using BiSplit algorithm.
}
