% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes.R
\name{run_bidirectional_enrichment}
\alias{run_bidirectional_enrichment}
\title{Run Gene Set Enrichment Analysis on Loadings}
\usage{
run_bidirectional_enrichment(
  loadings,
  comp,
  topN,
  species = "Homo sapiens",
  category = "H",
  pvalueCutoff = 0.05,
  doPlot = F
)
}
\arguments{
\item{loadings}{A numeric matrix or data frame of gene loadings with genes as row names.
Rows represent genes and columns represent different components.}

\item{comp}{Numeric. The component number (row index) to analyze.}

\item{topN}{Numeric. The number of top genes (from both positive and negative directions) to consider.}

\item{species}{Character. The species to query from \code{msigdbr}. Default is \code{"Homo sapiens"}.}

\item{category}{Character. The gene set category to use (e.g., \code{"H"} for Hallmark). Default is \code{"H"}.}

\item{pvalueCutoff}{Numeric. The p-value cutoff to consider enrichment as significant. Default is 0.05.}
}
\value{
A list with two elements:
\describe{
  \item{enrichedPos}{An object of class \code{enrichResult} with enrichment results for positively loaded genes.}
  \item{enrichedNeg}{An object of class \code{enrichResult} with enrichment results for negatively loaded genes.}
}
}
\description{
This function performs gene set enrichment analysis for both positive and negative loadings 
using the \code{clusterProfiler} package with gene sets obtained via \code{msigdbr}.
It extracts the top \code{topN} genes in both directions for a specified component 
and then uses \code{enricher} for enrichment analysis. If any significant enrichment is detected 
(based on the adjusted p-value cutoff), a dotplot is displayed.
}
\examples{
\dontrun{
  # Assuming results$loadings[[1]] is a numeric matrix with gene symbols as row names:
  loadings <- results$loadings[[1]]
  enrichmentResults <- run_bidirectional_enrichment(loadings, comp = 2, topN = 100)
}

}
