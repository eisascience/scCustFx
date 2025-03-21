#' Save a List of DE Results to an Excel File
#'
#' Saves a list of  data frames into an Excel file, with each list item as a separate sheet.
#'
#' @param DFLS A named list of data frames 
#' @param file_path File path for the output Excel file
#'
#' @return Saves an Excel file with multiple sheets.
#' @export
#'
#' @examples
#' save_de_results_to_excel(DEgeneLS, "DE_results.xlsx")
Save_LS_Excel <- function(DFLS, file_path = "DE_results.xlsx") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx: install.packages('openxlsx')")
  
  wb <- openxlsx::createWorkbook()
  
  for (name in names(DFLS)) {
    openxlsx::addWorksheet(wb, sheetName = name)
    openxlsx::writeData(wb, sheet = name, DFLS[[name]])
  }
  
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message("Excel file saved: ", file_path)
}