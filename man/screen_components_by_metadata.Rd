% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations_SDA.R
\name{screen_components_by_metadata}
\alias{screen_components_by_metadata}
\title{Screen All Components for Differences by Metadata Feature}
\usage{
screen_components_by_metadata(
  scores_matrix,
  meta_df,
  meta_feature,
  top_n = 5,
  method = "auto"
)
}
\arguments{
\item{scores_matrix}{A numeric matrix where rows are samples and columns are components.}

\item{meta_df}{A data frame with rownames matching `scores_matrix`, containing metadata.}

\item{meta_feature}{Character. The metadata column to group by (e.g., `"seqrun"`).}

\item{top_n}{Integer. Number of top components to visualize (default: 5).}

\item{method}{Character. Statistical test to use (`"anova"` for normal data, `"kruskal"` for non-parametric, default: `"auto"`).}
}
\value{
A list with:
\item{ranked_results}{A dataframe ranking components by significance (p-values & effect sizes).}
\item{plots}{A list of ggplot objects for top components.}
}
\description{
This function systematically analyzes all components in a score matrix, testing their differences across 
a specified metadata feature. Components are ranked by statistical significance.
}
\examples{
\dontrun{
screen_results <- screen_components_by_metadata(scores_matrix = results$scores, 
                                                meta_df = MetaDF, 
                                                meta_feature = "seqrun")
print(screen_results$ranked_results) # Ranked table of components
screen_results$plots[[1]]  # View first plot
}
}
