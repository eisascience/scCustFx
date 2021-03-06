% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Normalizations.R
\name{avgExpr_norm}
\alias{avgExpr_norm}
\title{avgExpr_norm}
\usage{
avgExpr_norm(mat, doGlobal = T, doPerFeat = T, rm.NA = T)
}
\arguments{
\item{mat}{A matrix where the features are the rows, cells/measured units are columns}
}
\value{
A normalized matrix
}
\description{
Normalize a matrix by average expression (total and or per feature)
}
