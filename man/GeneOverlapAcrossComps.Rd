% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDA_util.R
\name{GeneOverlapAcrossComps}
\alias{GeneOverlapAcrossComps}
\title{Gene overlap across components}
\usage{
GeneOverlapAcrossComps(GeneSet, compDF, Ngenes = 10, verbose = T)
}
\arguments{
\item{GeneSet}{a character vector of gene names to compare across the components}

\item{compDF}{a data.frame object containing gene expression values across components, where each column corresponds to a component and each row to a gene}

\item{Ngenes}{the number of genes to consider in the top overlap}

\item{verbose}{a boolean indicating whether to print the top found genes for each component}
}
\value{
A named list of length equal to the length of the input \code{GeneSet}, where each element of the list is a numeric vector of length equal to the number of components in \code{compDF} representing the number of overlapping genes in the top \code{Ngenes} of each component with the input gene set.
}
\description{
Given a set of genes and a dataframe of gene expression values across multiple components, this function computes the number of overlapping genes in the top \code{Ngenes} of each component with the input gene set.
}
