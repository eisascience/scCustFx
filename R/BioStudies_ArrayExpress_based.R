




#' Scrape ArrayExpress metadata and fastq.gz links
#'
#' This function scrapes metadata and fastq.gz links from ArrayExpress for a given accession.
#'
#' @param accession Character string specifying the ArrayExpress accession number.
#' @return No direct return value. The function saves metadata as a CSV file and fastq.gz links as a text file.
#' @seealso [rvest](https://cran.r-project.org/web/packages/rvest/index.html), [RSelenium](https://cran.r-project.org/web/packages/RSelenium/index.html)
#' @examples
#' \dontrun{
#' # Example usage:
#' scrape_arrayexpress("E-MTAB-1234")
#' }
#' @importFrom rvest html_nodes html_attr read_html
#' @importFrom RSelenium rsDriver
#' @export
scrape_arrayexpress <- function(accession, port = 4562L, browser = "firefox") {
  
  
  
  # Install and load necessary packages
  if (!requireNamespace("rvest", quietly = TRUE)) {
    install.packages("rvest")
  }
  if (!requireNamespace("RSelenium", quietly = TRUE)) {
    install.packages("RSelenium")
  }
  if (!requireNamespace("XML", quietly = TRUE)) {
    install.packages("XML")
  }
  
  library(rvest)
  library(RSelenium)
  library(XML)
  
  
  # Start Selenium server (make sure you have Java installed)
  driver <- rsDriver(browser = browser, port = port)
  remote_driver <- driver[["client"]]
  
  # Construct the URL
  url <- paste0("https://www.ebi.ac.uk/biostudies/arrayexpress/studies/", accession, "/sdrf")
  
  print(paste0("Trying URL: ", url))
  
  # Open the URL using Selenium
  remote_driver$navigate(url)
  
  dropdown <- remote_driver$findElement(using = 'name', value = 'sdrf_length')
  dropdown$clickElement()
  
  # Locate and click the option with the value 'All'
  all_option <- remote_driver$findElement(using = 'xpath', value = '//option[@value="-1"]')
  all_option$clickElement()
  
  
  # Get the page source after the JavaScript has executed
  page_source <- remote_driver$getPageSource()[[1]]
  
 

  page <- htmlParse(page_source, asText = TRUE)
  tables <- XML::readHTMLTable(page)
  
  # Close the Selenium server
  remote_driver$close()
  
  metadata_table = tables[["sdrf"]]

  # Save metadata as CSV
  metadata_csv <- paste0(accession, "_metadata.csv")
  write.csv(metadata_table, file = metadata_csv, row.names = FALSE)
  cat("Metadata saved as:", metadata_csv, "\n")
  
  
  fastq_links <- read_html(page_source) %>%
    html_nodes("td.FASTQ a") %>%
    html_attr("href")
  
  
  # Save fastq.gz links as a text file
  fastq_txt <- paste0(accession, "_fastq_links.txt")
  writeLines(fastq_links, fastq_txt)
  cat("Fastq.gz links saved as:", fastq_txt, "\n")
}


# Example usage
scrape_arrayexpress(accession = "E-MTAB-11069")


