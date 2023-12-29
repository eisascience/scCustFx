




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
scrape_arrayexpress <- function(accession, 
                                port = 4562L, browser = "firefox",
                                LoadWaitSec=5,
                                R1R2only=T, splitData = F) {
  
  
  
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
  library(dplyr)
  
  
  # Start Selenium server (make sure you have Java installed)
  driver <- rsDriver(browser = browser, port = port)
  
  tryCatch({
    
    remote_driver <- driver[["client"]]
    
    # Construct the URL
    url <- paste0("https://www.ebi.ac.uk/biostudies/arrayexpress/studies/", accession, "/sdrf")
    
    print(paste0("Trying URL: ", url))
    
    # Open the URL using Selenium
    remote_driver$navigate(url)
    
    # Wait for the page to update
    Sys.sleep(LoadWaitSec)
    
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
    #remote_driver$close()
    
    metadata_table = tables[["sdrf"]]
    
    # Save metadata as CSV
    metadata_csv <- paste0(accession, "_metadata.csv")
    write.csv(metadata_table, file = metadata_csv, row.names = FALSE)
    cat("Metadata saved as:", metadata_csv, "\n")
    
    
    fastq_links <- read_html(page_source) %>%
      html_nodes("td.FASTQ a") %>%
      html_attr("href")
    
    assay_links <- read_html(page_source) %>%
      html_nodes("td.ASSAY a") %>%
      html_attr("href")
    
    metadata_table$FASTQ_URL = fastq_links
    
    
    
    
    if(R1R2only){
      
      if(length(grep("R1", fastq_links))>0){
        
        metadata_table =  metadata_table[sort(c(grep("_1", fastq_links), 
                                                grep("_2", fastq_links))),] 
        
        fastq_links = metadata_table$FASTQ_URL
        
      } else {
        
        if(length(grep("_1", fastq_links))>0){
          
          metadata_table =  metadata_table[sort(c(grep("_1", fastq_links), 
                                grep("_2", fastq_links))),] 
          
        } else {
          print("R1 R2 or _1 _2 not found in fastq files")
        }
      }
      
    }
    
    
    fastq_links = metadata_table$FASTQ_URL
    
   
    
    
    if (splitData) {
      # Check if "organism" and "organism part" are in the columns
      if ("organism" %in% colnames(metadata_table) && "organism part" %in% colnames(metadata_table)) {
        
        colnames(metadata_table) <- make.unique(colnames(metadata_table))
        
        ordered_table <- metadata_table %>%
          arrange(organism, `organism part`)
        
        
        # Get unique values of "organism" and "organism part"
        unique_organisms <- unique(ordered_table$organism)
        unique_organism_parts <- unique(ordered_table$`organism part`)
        
        # Create a list to store sub-tables
        sub_tables <- list()
        
        # Iterate over unique combinations of "organism" and "organism part"
        for (org in unique_organisms) {
          for (org_part in unique_organism_parts) {
            # Create a sub-table for the current combination
            sub_table <- ordered_table[ordered_table$organism == org & ordered_table$`organism part` == org_part, ]
            
            # Order the sub-table by individual and FASTQ_URL
            sub_table <- sub_table[order(sub_table$individual, sub_table$FASTQ_URL), ]
            
            # Add the sub-table to the list
            sub_tables[[paste(org, org_part, sep = "_")]] <- sub_table
          }
        }
        
        # Save sub-tables
        for (i in seq_along(sub_tables)) {
          sub_table <- sub_tables[[i]]
          sub_table_name <- names(sub_tables)[i]
          sub_table_csv <- paste0(accession, "_", sub_table_name, "_metadata.csv")
          write.csv(sub_table, file = sub_table_csv, row.names = FALSE)
          cat("Sub-table saved as:", sub_table_csv, "\n")
          
          # Save FASTQ links as a text file
          fastq_links <- sub_table$FASTQ_URL
          fastq_txt <- paste0(accession, "_", sub_table_name, "_fastq_links.txt")
          writeLines(fastq_links, fastq_txt)
          cat("Fastq.gz links saved as:", fastq_txt, "\n")
          
        }
      }
      
      
    } else {
      # Save fastq.gz links as a text file
      fastq_txt <- paste0(accession, "_fastq_links.txt")
      writeLines(fastq_links, fastq_txt)
      cat("Fastq.gz links saved as:", fastq_txt, "\n")
    }
    
    
  }, error = function(e) {
    message("Error: ", e$message)
  }, finally = {
    # Close the Selenium server
    remote_driver$close()
  })
  
  
}


# Example usage
# scCustFx:::scrape_arrayexpress(accession = "E-MTAB-8122", port = 4562L)

# scrape_arrayexpress(accession = "E-MTAB-8122", port = 4562L, splitData = T)

