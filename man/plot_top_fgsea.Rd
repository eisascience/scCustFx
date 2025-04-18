% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genes.R
\name{plot_top_fgsea}
\alias{plot_top_fgsea}
\title{Plot Top Enriched Pathways from fgsea Results}
\usage{
plot_top_fgsea(
  fgsea_res,
  ranked_genes,
  gene_sets,
  top_n = 3,
  absNES = 1.5,
  padj = 0.05,
  filterSig = T
)
}
\arguments{
\item{fgsea_res}{Data frame output from `fgsea()`.}

\item{ranked_genes}{Named vector of ranked log2 fold changes.}

\item{gene_sets}{List of gene sets used in fgsea.}

\item{top_n}{Integer. Number of pathways to plot.}
}
\value{
List of ggplot enrichment plots.
}
\description{
Plot Top Enriched Pathways from fgsea Results
}
