% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDA_processing_based.R
\name{SDApostprocessing}
\alias{SDApostprocessing}
\title{Post-process Single-cell Data Analysis (SDA) Results}
\usage{
SDApostprocessing(
  resultsDir,
  outputFolder,
  normedDGE,
  orgDb = "org.Hs.eg.db",
  storeGoEnrichment = T
)
}
\arguments{
\item{resultsDir}{A character string specifying the directory where the SDA analysis results are stored.}

\item{outputFolder}{A character string specifying the path to the output folder where additional
processed results and analyses will be saved.}

\item{normedDGE}{A matrix or data frame containing the normalized gene expression data used in the
SDA analysis. Row names should correspond to features (genes), and column names should correspond to cell barcodes.}

\item{storeGoEnrichment}{A logical value indicating whether to perform and store Gene Ontology (GO)
enrichment analysis on the analysis components. Default is `TRUE`.}
}
\value{
A list containing the enriched SDA results, including cell barcodes, feature names, component
  statistics, and optionally GO enrichment analysis results.
}
\description{
This function loads the results from an SDA analysis, enriches them with additional metadata such
as cell barcodes and feature (gene) names, computes component statistics, and optionally performs
Gene Ontology (GO) enrichment analysis on the components. It is designed to work with results
generated by the SDAtools and provides a comprehensive overview of the analysis for further interpretation.
}
