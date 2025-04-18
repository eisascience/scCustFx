% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils_tables.R
\name{table2DF}
\alias{table2DF}
\title{Covert a table to a dataframe}
\usage{
table2DF(tbl = NULL, do.prop = F, prop.marg = 1)
}
\arguments{
\item{tbl}{table}

\item{do.prop}{if T, does prop and need to set prop.marg}

\item{prop.marg}{either 1 or 2 passed to prop.table margin}
}
\value{
a dataframe version of the input table
}
\description{
Covert a table to a dataframe
}
