% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_based.R
\name{RunComBatseq}
\alias{RunComBatseq}
\title{RunComBatseq}
\usage{
RunComBatseq(
  SerObj = NULL,
  batch = "batch",
  runningTime = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{SerObj}{A seurat SerObj to correct batch effects.}

\item{batch}{Batch labels.}

\item{runningTime}{Return the running time.}

\item{verbose}{Print verbose.}

\item{...}{Arguments passed to other methods.}
}
\value{
Corrected and normalized Seurat SerObj.
}
\description{
RunComBatseq
}
