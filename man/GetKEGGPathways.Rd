% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KEGG_Analysis.R
\name{GetKEGGPathways}
\alias{GetKEGGPathways}
\title{GetKEGGPathways}
\usage{
GetKEGGPathways(df = NULL, species = "human")
}
\arguments{
\item{df}{data frame containing differential gene expression results}

\item{species}{species to map . Available option: "human" and "mouse".}
}
\value{
Data frame containing KEGG pathways
}
\description{
GetKEGGPathways
}
