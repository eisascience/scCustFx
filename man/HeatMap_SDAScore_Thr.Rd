% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations_SDA.R
\name{HeatMap_SDAScore_Thr}
\alias{HeatMap_SDAScore_Thr}
\title{This function generates a Heatmap of thresholding SDA score}
\usage{
HeatMap_SDAScore_Thr(
  SDAscore,
  Meta = NULL,
  GT = T,
  CutThresh = 0,
  clustering_method = "ward.D2"
)
}
\arguments{
\item{SDAscore}{vector of SDA score}

\item{Meta}{a meta vector or NULL}

\item{GT}{if T greater than}

\item{CutThresh}{Cut the score based on threshold .}
}
\value{
A ggplot object representing the loadings along genomic coordinates.
}
\description{
This function generates a Heatmap of thresholding SDA score
}
