% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalizations.R
\name{normaliseDGEv2}
\alias{normaliseDGEv2}
\title{normaliseDGEv2}
\usage{
normaliseDGEv2(
  dge,
  verbose = FALSE,
  center = TRUE,
  scale = TRUE,
  outliers = FALSE,
  threshold = 10,
  min_library_size = 0,
  gene_subset = 0.8,
  transformation = "sqrt"
)
}
\arguments{
\item{dge}{A matrix rows of genes and cols of cells}
}
\value{
A quantile-based normalized matrix
}
\description{
a sparse implementation of dropsim::normaliseDGE()
}
