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



#' Save DE Results to Excel, Splitting by Up/Down Regulation
#'
#' This function takes a named list of DESeq2 results data frames (`DEgeneLS`),
#' splits each into **Upregulated (log2FoldChange > 0) and Downregulated (log2FoldChange < 0)**
#' genes, and saves them in an Excel workbook, creating separate sheets for each.
#'
#' @param DFLS Named list of data frames containing DE results.
#' @param file_path File path to save the Excel file (default: "DE_results.xlsx").
#' @param logfc_col Name of the column containing log2 fold changes (default: "log2FoldChange").
#'
#' @return Saves an Excel file with organized sheets.
#' @export
#'
#' @examples
#' Save_LS_ExcelV2(DEgeneLS, file_path = "DE_results.xlsx")
Save_LS_ExcelV2 <- function(DFLS, file_path = "DE_results.xlsx", logfc_col = "log2FoldChange") {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Please install openxlsx: install.packages('openxlsx')")
  
  wb <- openxlsx::createWorkbook()
  
  for (name in names(DFLS)) {
    df <- DFLS[[name]]
    
    # Check if the column exists
    if (!(logfc_col %in% colnames(df))) {
      warning(paste("Skipping", name, "—", logfc_col, "column not found."))
      next
    }
    
    # Split into upregulated and downregulated
    df_up <- df[df[[logfc_col]] > 0, ]
    df_down <- df[df[[logfc_col]] < 0, ]
    
    # Add both as separate sheets
    if (nrow(df_up) > 0) {
      sheet_name_up <- paste0(name, "_Up")
      openxlsx::addWorksheet(wb, sheetName = sheet_name_up)
      openxlsx::writeData(wb, sheet = sheet_name_up, df_up)
    }
    
    if (nrow(df_down) > 0) {
      sheet_name_down <- paste0(name, "_Down")
      openxlsx::addWorksheet(wb, sheetName = sheet_name_down)
      openxlsx::writeData(wb, sheet = sheet_name_down, df_down)
    }
  }
  
  openxlsx::saveWorkbook(wb, file_path, overwrite = TRUE)
  message("Excel file saved: ", file_path)
}