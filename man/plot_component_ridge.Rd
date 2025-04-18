% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualizations_SDA.R
\name{plot_component_ridge}
\alias{plot_component_ridge}
\title{Generate Ridge Plot of Component Scores by Metadata Feature}
\usage{
plot_component_ridge(
  scores_matrix,
  meta_df,
  component,
  meta_feature,
  title = "Ridge Plot of Component Scores"
)
}
\arguments{
\item{scores_matrix}{A numeric matrix where rows are samples (with rownames) and columns are components.}

\item{meta_df}{A data frame with rownames matching those in `scores_matrix` and containing metadata features.}

\item{component}{Integer or character. The component column (e.g., `"Component_1"`) to visualize.}

\item{meta_feature}{Character. The metadata column to split the ridge plot by (e.g., `"seqrun"`).}

\item{title}{Character. Plot title (default: `"Ridge Plot of Component Scores"`).}
}
\value{
A ggplot object displaying the ridge plot.
}
\description{
This function creates a ridge plot to visualize the distribution of component scores
split by a specified metadata feature.
}
\examples{
\dontrun{
if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggridges", quietly = TRUE)) {
    # Example scores matrix
    scores <- matrix(rnorm(500), nrow = 100, ncol = 5)
    rownames(scores) <- paste0("Sample", 1:100)
    colnames(scores) <- paste0("Component_", 1:5)

    # Example metadata dataframe
    meta <- data.frame(
        seqrun = sample(c("Run1", "Run2", "Run3"), 100, replace = TRUE),
        InfectionGroup = sample(c("Control", "Infected"), 100, replace = TRUE),
        MonkeyName = sample(LETTERS[1:5], 100, replace = TRUE),
        TreatmentCode = sample(c("A", "B", "C"), 100, replace = TRUE),
        PD1status = sample(c("High", "Low"), 100, replace = TRUE),
        row.names = rownames(scores)
    )

    # Generate ridge plot for Component 1 split by seqrun
    plot_component_ridge(scores, meta, component = "Component_1", meta_feature = "seqrun")
}
}
}
