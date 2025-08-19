#' Prepare Gene Expression (GEX) Matrix with Full Cell Identifiers and Optional Prefix
#'
#' This function constructs a cleaned GEX count matrix from raw CosMx SMI transcript data,
#' ensuring that full cell IDs (c_<Slide>_<FOV>_<CellId>) are preserved. If the raw data
#' already contains a valid `cell` column, it will be used directly; otherwise the function
#' will reconstruct full cell IDs from `Slide`, `fov`, and `CellId` columns. You can also
#' provide a `file_prefix` to save the matrix under a sub-directory or filename prefix.
#'
#' @param raw_transcripts A data frame or a character string path to a CSV file containing
#'   raw transcript data. Must include either a `cell` column in the format
#'   `c_<Slide>_<FOV>_<CellId>` or the separate columns `Slide`, `fov`, and `CellId`.
#' @param z_stacks An integer vector of Z‑slice indices to retain in the analysis.
#' @param output_dir Optional character string specifying a directory to save the
#'   cleaned expression matrix when `save_matrix = TRUE`.
#' @param save_matrix Logical; if `TRUE`, the resulting count matrix is written as
#'   `exprMat_zclean.csv` under `output_dir` (and `file_prefix` if provided). Defaults to `FALSE`.
#' @param file_prefix Optional character string to prefix the output file path. Can be used as a
#'   sub-directory name or filename prefix. Defaults to an empty string (no prefix).
#' @param slide_id Single integer.  Which slide are we on?  Defaults to 1.
#' 
#' @return A data frame (matrix) of counts where row names are full cell IDs and column
#'   names are transcript targets.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' # Using an in-memory data frame and save without prefix
#' raw_df <- read.csv("raw_transcripts.csv", stringsAsFactors = FALSE)
#' gex_mat <- PrepareGEXMatrix(raw_df, z_stacks = 1:5, save_matrix = TRUE, output_dir = "results/")
#'
#' # Using a file prefix for saving under a sub-directory
#' gex_mat <- PrepareGEXMatrix(
#'   raw_transcripts = "raw_transcripts.csv",
#'   z_stacks       = c(1,2,3),
#'   output_dir     = "results/",
#'   save_matrix    = TRUE,
#'   file_prefix    = "exp1"
#' )
#' }
PrepareGEXMatrix <- function(
    raw_transcripts,
    z_stacks,
    output_dir     = NULL,
    save_matrix    = FALSE,
    file_prefix    = "",
    slide_id    = 1
) {
  library(dplyr)
  
  
  # --- load data ---
  if (is.character(raw_transcripts)) {
    raw_df <- read.csv(raw_transcripts, stringsAsFactors = FALSE)
  } else if (is.data.frame(raw_transcripts)) {
    raw_df <- raw_transcripts
  } else {
    stop("'raw_transcripts' must be a file path or data frame.")
  }
  
  # --- filter ---
  raw_df <- raw_df %>%
    filter(CellId != 0, z %in% z_stacks) %>%
    filter(!grepl("SystemControl", target))
  
  # --- construct full cell IDs from fov + CellId (and slide_id) ---
  if (! all(c("fov","CellId") %in% colnames(raw_df))) {
    stop("Cannot determine cell IDs: need both 'fov' and 'CellId' columns.")
  }
  raw_df$cell <- paste0("c_", slide_id, "_", raw_df$fov, "_", raw_df$CellId)
  
  # --- build the GEX matrix ---
  GEX_clean_z_stack <- as.data.frame.matrix(
    table(raw_df$cell, raw_df$target),
    stringsAsFactors = FALSE
  )
  
  # --- optionally save ---
  if (save_matrix) {
    if (is.null(output_dir)) {
      stop("Please specify output_dir when save_matrix = TRUE.")
    }
    csv_path <- file.path(output_dir, file_prefix, "exprMat_zclean.csv")
    dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
    write.csv(GEX_clean_z_stack, csv_path, row.names = TRUE)
    message("GEX matrix saved to: ", csv_path)
  }

  
  return(GEX_clean_z_stack)
}
