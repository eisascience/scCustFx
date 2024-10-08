% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_based.R
\name{RunSilhouette}
\alias{RunSilhouette}
\title{RunSilhouette}
\usage{
RunSilhouette(SerObj = NULL, batch = "celltype", reduction = "pca", dims = 10)
}
\arguments{
\item{SerObj}{A seurat SerObj to correct batch effects.}

\item{batch}{Batch labels.}

\item{reduction}{Reduction to use.}

\item{dims}{Number of dimensions to use.}
}
\value{
Silhouette width score.
}
\description{
RunSilhouette
}
