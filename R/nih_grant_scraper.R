

#' Parse a NIH Grant Page for Key Dates and Metadata
#'
#' This function reads an NIH grant opportunity webpage from a given URL and extracts the key deadlines
#' from the 'Key Dates' table as well as relevant metadata such as the Funding Opportunity Number (FON),
#' grant title, activity code (with a fallback from the title if missing), and an indicator for early investigator
#' status (including checks for "early investigator", "early stage", "early-stage", or "ESI").
#'
#' @param url A character string containing the URL of the NIH grant opportunity page.
#'
#' @return A tibble with one or more rows containing deadline information and associated metadata. The returned
#' tibble contains the following columns:
#'   \item{url}{The URL of the grant page.}
#'   \item{date}{The parsed deadline date (Date object).}
#'   \item{type}{The deadline type (e.g. "NEW" or "AIDS").}
#'   \item{fon}{The Funding Opportunity Number.}
#'   \item{grant_title}{The title of the funding opportunity.}
#'   \item{activity_code}{The activity code (e.g. R01, R21); extracted from the usual location or, if missing, from the grant title.}
#'   \item{early_investigator}{Indicator ("Yes"/"No") showing if the opportunity is targeted to early investigators.}
#'
#' @examples
#' \dontrun{
#'   # Example URL:
#'   url <- "https://grants.nih.gov/grants/guide/pa-files/PAR-25-362.html"
#'   parse_single_url(url)
#' }
#'
#' @export
parse_single_url <- function(url) {
  
  library(rvest)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  
  page <- tryCatch(
    read_html(url),
    error = function(e) {
      message("Error reading URL: ", url, " -- ", e$message)
      return(NULL)
    }
  )
  if (is.null(page)) return(NULL)
  
  # ---- Extract Metadata ----
  
  # Funding Opportunity Number (FON)
  fon <- page %>% 
    html_element("span.noticenum") %>% 
    html_text(trim = TRUE)
  
  # Grant Title (Funding Opportunity Title)
  grant_title <- page %>% 
    html_element("span.title") %>% 
    html_text(trim = TRUE)
  
  # Activity Code / Type (e.g., R01); usually found in row with data-index="7"
  activity_code <- page %>% 
    html_element("div[data-index='7'] a") %>% 
    html_text(trim = TRUE)
  
  # If the activity_code is missing or blank, try to extract it from grant_title.
  if (is.na(activity_code) || str_trim(activity_code) == "") {
    # Look for a pattern like "R01", "R21", etc. (case-insensitive)
    extracted <- str_extract(grant_title, "(?i)R\\d+")
    if (!is.na(extracted)) {
      activity_code <- toupper(extracted)
    } else {
      activity_code <- ""  # or "Not Specified"
    }
  }
  
  # Determine if this is for early investigators by scanning the entire page text.
  page_text <- page %>% html_text()
  
  # early_investigator <- ifelse(str_detect(tolower(page_text), "early investigator"), "Yes", "No")
  early_investigator <- ifelse(
    str_detect(tolower(grant_title), "early[ -]?investigator|early[ -]?stage|esi"),
    "Yes",
    "No"
  )  

  
  # ---- Extract Key Dates Table ----
  table_node <- page %>% html_element("table#keyDatesContentTable")
  if (is.null(table_node)) {
    message("No matching table found for: ", url)
    # Return metadata only if deadlines table is missing.
    return(tibble(
      url = url,
      fon = fon,
      grant_title = grant_title,
      activity_code = activity_code,
      early_investigator = early_investigator
    ))
  }
  
  # Parse the table without treating any row as header.
  raw_tbl <- table_node %>% html_table(header = FALSE, fill = TRUE)
  if (nrow(raw_tbl) < 2) {
    message("Table doesn't contain enough rows: ", url)
    return(NULL)
  }
  
  # Assume that the second row contains the header (due types),
  # then drop the first two rows from the raw table.
  hdr <- raw_tbl[2, ]
  data_tbl <- raw_tbl[-1, ][-1, ]
  names(data_tbl) <- as.character(unlist(hdr))
  
  # Pivot the deadlines table so that each cell is paired with its header.
  long_data <- data_tbl %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = -row_id,
      names_to = "due_type",
      values_to = "due_date"
    ) %>%
    filter(!is.na(due_date) & due_date != "")
  
  # Remove rows containing "Not Applicable"
  long_data <- long_data %>%
    filter(!str_detect(due_date, regex("not applicable", ignore_case = TRUE)))
  
  # Extract the date from cells (expected format e.g. "May 07, 2025 *")
  long_data <- long_data %>%
    mutate(
      date_str    = str_extract(due_date, "([A-Za-z]+\\s+\\d{1,2},\\s+\\d{4})"),
      date_parsed = as.Date(date_str, format = "%B %d, %Y")
    ) %>%
    filter(!is.na(date_parsed))
  
  # Use the header field to set a deadline type label.
  long_data <- long_data %>%
    mutate(
      deadline_type = case_when(
        str_to_lower(due_type) == "new" ~ "NEW",
        str_detect(str_to_lower(due_type), "aids") ~ "AIDS",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(deadline_type))
  
  # Remove duplicates and add the URL column.
  cleaned_deadlines <- long_data %>%
    distinct(date_parsed, deadline_type, .keep_all = TRUE) %>%
    mutate(url = url) %>%
    select(url, date = date_parsed, type = deadline_type)
  
  # Attach metadata to each deadline row.
  enriched_deadlines <- cleaned_deadlines %>%
    mutate(
      fon = fon,
      grant_title = grant_title,
      activity_code = activity_code,
      early_investigator = early_investigator
    )
  
  return(enriched_deadlines)
}

#' Parse Multiple NIH Grant Pages and Combine Key Dates Data
#'
#' This function takes a vector of NIH grant opportunity page URLs, calls \code{parse_single_url}
#' for each URL to extract both deadline information and associated metadata,
#' and then combines the results into a single wide-format tibble.
#' The resulting tibble contains the following metadata columns:
#' \item{url}{The grant page URL.}
#' \item{fon}{Funding Opportunity Number.}
#' \item{grant_title}{The title of the funding opportunity.}
#' \item{activity_code}{Activity code (e.g., R01, R21).}
#' \item{early_investigator}{Indicator ("Yes"/"No") for early investigator targeting.}
#'
#' In addition, the tibble includes one column per unique deadline date (formatted as YYYY-MM-DD)
#' with cell values indicating the deadline type ("NEW" or "AIDS").
#'
#' @param urls A character vector of NIH grant opportunity page URLs.
#'
#' @return A tibble in wide format with one row per grant. The left-hand side columns contain metadata,
#' and the remaining columns represent the deadlines for that grant.
#'
#' @examples
#' \dontrun{
#'   urls <- c("https://grants.nih.gov/grants/guide/pa-files/PAR-25-362.html",
#'             "https://grants.nih.gov/grants/guide/pa-files/PAR-23-198.html")
#'   parse_nih_due_dates(urls)
#' }
#'
#' @export
parse_nih_due_dates <- function(urls) {
  results_list <- lapply(urls, parse_single_url)
  combined <- bind_rows(results_list)
  
  if (nrow(combined) == 0) {
    warning("No data parsed from any URLs.")
    return(NULL)
  }
  
  # Pivot the deadlines information to wide format.
  # The metadata (url, fon, grant_title, activity_code, early_investigator) become
  # the left-hand side columns; then each unique date becomes a column.
  final_df <- combined %>%
    group_by(url, fon, grant_title, activity_code, early_investigator, date) %>%
    summarize(type = first(type), .groups = "drop") %>%
    mutate(date_col = format(date, "%Y-%m-%d")) %>%
    select(url, fon, grant_title, activity_code, early_investigator, date_col, type) %>%
    pivot_wider(names_from = date_col, values_from = type)
  
  # Reorder columns so metadata appear on the left.
  meta_cols <- c("url", "fon", "grant_title", "activity_code", "early_investigator")
  other_cols <- setdiff(names(final_df), meta_cols)
  final_df <- final_df[, c(meta_cols, sort(other_cols))]
  
  return(final_df)
}

# Example usage:
urls <- c(
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-362.html", 
  "https://grants.nih.gov/grants/guide/pa-files/PAR-23-198.html", 
  "https://grants.nih.gov/grants/guide/pa-files/PAR-24-075.html", 
  "https://grants.nih.gov/grants/guide/pa-files/PA-25-304.html", 
  "https://grants.nih.gov/grants/guide/pa-files/PAR-23-297.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-165.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-228.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-229.html",
  "https://grants.nih.gov/grants/guide/pa-files/PA-25-169.html",
  "https://grants.nih.gov/grants/guide/pa-files/PA-25-301.html",
  "https://grants.nih.gov/grants/guide/pa-files/PA-25-305.html",
  "https://grants.nih.gov/grants/guide/pa-files/PA-25-304.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-238.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-131.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-25-235.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-23-236.html",
  "https://grants.nih.gov/grants/guide/pa-files/PAR-23-169.html"
)


#' Write NIH Deadline Data to Excel
#'
#' This function takes a character vector of NIH grant page URLs and an output file path, calls
#' \code{parse_nih_due_dates} to extract key dates and metadata from each grant page, transforms
#' the result into a transposed matrix (so that deadlines appear as rows and grants as columns),
#' and writes the resulting data frame to an Excel file using the \code{xlsx} package.
#'
#' The transformation process involves:
#' \enumerate{
#'   \item Converting the resulting wide tibble (with one row per grant) into a matrix.
#'   \item Setting the URL values as row names.
#'   \item Replacing \code{NA} values with empty strings for aesthetics.
#'   \item Transposing the matrix so that deadline dates become row labels and metadata columns become columns.
#'   \item Converting the transposed matrix back to a data frame and writing to an Excel file.
#' }
#'
#' @param urls A character vector of NIH grant opportunity page URLs.
#' @param output_path A character string specifying the full file path (including filename)
#'   for the output Excel file.
#'
#' @return Invisibly returns the final data frame written to Excel.
#'
#' @examples
#' \dontrun{
#'   urls <- c("https://grants.nih.gov/grants/guide/pa-files/PAR-25-362.html", 
#'             "https://grants.nih.gov/grants/guide/pa-files/PAR-23-198.html")
#'   write_nih_deadlines_excel(urls, "./NIH_deadlines.xlsx")
#' }
#'
#' @export
write_nih_deadlines_excel <- function(urls, output_path) {
  # Parse the NIH grant pages to get metadata and deadlines in wide format.
  result <- parse_nih_due_dates(urls)
  
  if (is.null(result) || nrow(result) == 0) {
    warning("No data was parsed from the provided URLs.")
    return(invisible(NULL))
  }
  
  # Create a matrix for output: Remove the 'url' column, assign URL values as row names.
  result_matrix <- as.matrix(result[,-1])
  rownames(result_matrix) <- result$url
  
  # Replace NA with an empty string for a cleaner appearance.
  result_matrix[is.na(result_matrix)] <- ""
  
  # Transpose the matrix if you prefer deadlines as rows and grants as columns.
  result_matrix <- t(result_matrix)
  
  # Convert the transposed matrix back to a data frame.
  result_df <- as.data.frame(result_matrix, check.names = FALSE)
  
  # Print the head of the final data frame for confirmation.
  print(head(result_df))
  
  # Write the final data frame to an Excel file.
  xlsx::write.xlsx(result_df, output_path, row.names = TRUE)
  
  return(result_df)
}




# urls <- sort(urls, decreasing = TRUE)
# result <- parse_nih_due_dates(urls)
# 
# # Optionally, if you wish to create a matrix for output:
# # Remove the 'url' column (or keep if desired) and assign row names.
# result_matrix <- as.matrix(result[,-1])
# rownames(result_matrix) <- result$url
# 
# # Replace NA with an empty string for aesthetics.
# result_matrix[is.na(result_matrix)] <- ""
# 
# # If you want to transpose so deadlines appear as rows and grants as columns:
# result_matrix <- t(result_matrix)
# result_df <- as.data.frame(result_matrix, check.names = FALSE)
# 
# # Print the final data frame.
# head(result_df)
# 
# # Write the result to an Excel file using the xlsx package.
# xlsx::write.xlsx(result_df, "./NIH_deadlines.xlsx", row.names = TRUE)
