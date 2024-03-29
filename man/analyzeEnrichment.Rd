% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations_Seurat_based.R
\name{analyzeEnrichment}
\alias{analyzeEnrichment}
\title{Analyze PC Enrichment}
\usage{
analyzeEnrichment(
  SerObj,
  CompN = 1,
  topN = 20,
  qrhub = "org.Hs.eg.db",
  invertWeights = F,
  labsize = 3.5,
  fontsize = 0.25,
  base_size = 20,
  max.overlaps = 20,
  MaxNsig = 30,
  reduction = "pca",
  redLab = "PC",
  topNGenes = 100,
  bottomNGenes = 100
)
}
\arguments{
\item{SerObj}{A Seurat object containing the data.}

\item{CompN}{The principal component number.}

\item{topN}{The number of top genes to consider for enrichment analysis.}

\item{qrhub}{The gene database used for enrichment analysis.}
}
\value{
A list containing two ggplot objects for positive and negative loadings.
}
\description{
This function performs enrichment analysis for the genes associated with positive
and negative loadings of a given principal component (PC) in a Seurat object.
}
\examples{
\dontrun{
result <- analyzeEnrichment(SerObj, CompN = 1, topN = 20, qrhub = "org.Hs.eg.db")
result$ggGO_pos | result$ggGO_neg
}
}
