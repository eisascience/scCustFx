#' Counts the overlap of items in rows per column feature
#'
#' This function takes in a list of ggplots and makes a GIF
#'
#' @param df a dataframe
#' @return rturns a count matrix
#' @export
count_overlap <- function(df) {
  overlap <- unique(unlist(df))
  
  count_matrix <- matrix(0, nrow = length(overlap), ncol = ncol(df), dimnames = list(overlap, colnames(df)))
  
  for (i in 1:ncol(df)) {
    overlapping_genes <- intersect(overlap, df[, i])
    count_matrix[overlapping_genes, i] <- rep(1, length(overlapping_genes))
  }
  return(count_matrix)
}

#' @title transposedt
#'
#' @description Transpose a data.table
#' @return A transposed data.table
#' @keywords transpose, t
#' @param dt The datatable
#' @import data.table
#' @export
TransposeDT <- function(dt, varlabel="myVar") {
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}




#' Generate and Save/Display a Styled Table
#'
#' This function creates a styled table from a given dataframe and either displays it 
#' or saves it as an image.
#'
#' @param df A dataframe to be formatted into a table.
#' @param file_name A string specifying the name of the file if saving the table.
#' @param save A logical value indicating whether to save the table (`TRUE`) or just display it (`FALSE`). Default is `TRUE`.
#'
#' @return If `save = FALSE`, the table is printed in the console. If `save = TRUE`, it is saved to the specified file.
#' @export
#'
#' @examples
#' SaveMy_table(mtcars, "table.png", save = TRUE)
#' SaveMy_table(iris, "table.html", save = FALSE)
SaveMy_table <- function(df, file_name, save = TRUE) {
  library("knitr")
  library("kableExtra")
  # Generate the table
  kable_obj <- kable(df, format = "html") %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
      full_width = FALSE, 
      position = "center",
      font_size = 14
    ) %>%
    row_spec(0, bold = TRUE, color = "white", background = "#0072B2")
  
  # Save or print the table
  if (save) {
    save_kable(kable_obj, file_name)
  } else {
    print(kable_obj)
  }
}
