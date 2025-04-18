% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nih_grant_scraper.R
\name{write_nih_deadlines_excel}
\alias{write_nih_deadlines_excel}
\title{Write NIH Deadline Data to Excel}
\usage{
write_nih_deadlines_excel(urls, output_path)
}
\arguments{
\item{urls}{A character vector of NIH grant opportunity page URLs.}

\item{output_path}{A character string specifying the full file path (including filename)
for the output Excel file.}
}
\value{
Invisibly returns the final data frame written to Excel.
}
\description{
This function takes a character vector of NIH grant page URLs and an output file path, calls
\code{parse_nih_due_dates} to extract key dates and metadata from each grant page, transforms
the result into a transposed matrix (so that deadlines appear as rows and grants as columns),
and writes the resulting data frame to an Excel file using the \code{xlsx} package.
}
\details{
The transformation process involves:
\enumerate{
  \item Converting the resulting wide tibble (with one row per grant) into a matrix.
  \item Setting the URL values as row names.
  \item Replacing \code{NA} values with empty strings for aesthetics.
  \item Transposing the matrix so that deadline dates become row labels and metadata columns become columns.
  \item Converting the transposed matrix back to a data frame and writing to an Excel file.
}
}
\examples{
\dontrun{
  urls <- c("https://grants.nih.gov/grants/guide/pa-files/PAR-25-362.html", 
            "https://grants.nih.gov/grants/guide/pa-files/PAR-23-198.html")
  write_nih_deadlines_excel(urls, "./NIH_deadlines.xlsx")
}

}
