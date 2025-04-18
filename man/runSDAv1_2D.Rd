% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDA_processing_based.R
\name{runSDAv1_2D}
\alias{runSDAv1_2D}
\title{Run SDA V1.1 in 2D i.e., matrix form}
\usage{
runSDAv1_2D(
  path.sda,
  resultsDir,
  rawDataFile,
  numComps = 30,
  max_iter = 100,
  save_freq = 20,
  set_seed = 1234,
  n_cells,
  nThreads = 6,
  num_blocks = 6
)
}
\arguments{
\item{path.sda}{A character string specifying the path to the SDA software executable.}

\item{resultsDir}{A character string specifying the directory where the SDA results will be saved.}

\item{rawDataFile}{A character string specifying the path to the preprocessed and normalized gene
expression data file to be analyzed.}

\item{numComps}{Integer, the number of components to be used in the analysis, default is 30.}

\item{max_iter}{Integer, the maximum number of iterations for the SDA algorithm, default is 100.}

\item{save_freq}{Integer, frequency of saving interim results, default is 20 iterations.}

\item{set_seed}{Integer, the seed for random number generation to ensure reproducibility, default is 1234.}

\item{n_cells}{Integer, the number of cells in the dataset to be analyzed.}

\item{nThreads}{Integer, the number of threads for parallel processing, default is 6.
If greater than 1, eigen decomposition will be done in parallel.}

\item{num_blocks}{Integer, the number of blocks to divide the dataset into for parallel processing, default is 6.}
}
\value{
Invisible NULL. The function is called for its side effects, including the execution
  of the SDA analysis and saving of results to the specified directory.
}
\description{
Executes the SDA algorithm on preprocessed and normalized single-cell gene expression data,
allowing for dimensionality reduction and analysis. This function is designed to work with
the specified version of SDA, taking in several parameters to customize the analysis process,
including the number of components, iteration settings, and parallel processing options.
}
