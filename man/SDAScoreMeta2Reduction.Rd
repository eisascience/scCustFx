% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDA_Seurat.R
\name{SDAScoreMeta2Reduction}
\alias{SDAScoreMeta2Reduction}
\title{Create a new dimension reduction object using SDA scores}
\usage{
SDAScoreMeta2Reduction(
  SerObj,
  loadingMat = NULL,
  sdaComps = NULL,
  reduction.key = "SDA_",
  assayName = "RNA",
  reduction.name = "SDA",
  includeLoading = F
)
}
\arguments{
\item{SerObj}{A Seurat object containing SDA scores to use for the reduction.}

\item{loadingMat}{A matrix of loadings to use for the reduction. If not specified, the loadings will be calculated from the SDA scores.}

\item{sdaComps}{A vector of column names or indices from SerObj@meta.data that correspond to the SDA scores to use for the reduction.}

\item{reduction.key}{A prefix string to use for the key of the new dimension reduction object.}

\item{assayName}{The name of the assay in SerObj to use for the reduction.}

\item{reduction.name}{The name to use for the new reduction in the Seurat object.}
}
\value{
A modified Seurat object with a new dimension reduction object added.
}
\description{
This function creates a new dimension reduction object using SDA scores from a Seurat object, and adds it to the Seurat object as a new assay.
}
\examples{
## Create a new SDA dimension reduction object from a Seurat object
my_seurat <- SDAScoreMeta2Reduction(SerObj = my_seurat, sdaComps = c("SDA_1", "SDA_2"))

}
