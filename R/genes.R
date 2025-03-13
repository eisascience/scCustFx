
#' Get.Gene_Info
#'
#' This function get gene info via biomaRt
#'
#' @param gene_names  gene set
#' @param dataset  default hsapiens_gene_ensembl, try mmusculus_gene_ensembl
#' @param attributes attributes to grab
#' @return converted human gene set
#' @export
Get.Gene_Info <- function(gene_names,
                          dataset = "hsapiens_gene_ensembl",
                          attributes = c("external_gene_name", "chromosome_name", 
                                         "start_position", "end_position", "description")) {
  library(biomaRt)
  # Connect to the Ensembl database
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Define the attributes we want to retrieve
  
  
  # Retrieve the information for the given gene names
  results <- getBM(attributes = attributes, 
                   filters = "external_gene_name", 
                   values = gene_names, 
                   mart = ensembl)
  
  # Return the results as a dataframe
  return(results)
}



#' convert.Mouse2Human.biomaRt
#'
#' This function converts mouse genes to human genes via biomaRt
#'
#' @param x mouse gene set
#' @return converted human gene set
#' @export
convert.Mouse2Human.biomaRt <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#' convert.Human2Mouse.biomaRt
#'
#' This function converts human genes to mouse genes via biomaRt
#'
#' @param x mouse gene set
#' @return converted human gene set
#' @export
convert.Human2Mouse.biomaRt <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

#' convert.Human2Mouse
#'
#' This function converts human genes to mouse genes via text manipulation e.g. CD8A -> Cd8as
#'
#' @param geneV mouse gene set
#' @return converted human gene set
#' @export
convert.Human2Mouse <- function(geneV){
  sub('^(\\w?)', '\\U\\1', tolower(geneV), perl=T)
}


#' clean gene set of commonly unwanted genes
#'
#' This function cleans the input gene set
#'
#' @param geneV mouse gene set
#' @return cleaned mouse gene set
#' @export
cleanmygenes.mouse = function(geneV, keepV=NULL){
  geneV = geneV[!grepl("Rps", geneV)]
  geneV = geneV[!grepl("Rpl", geneV)]
  # geneV = geneV[!grepl("LOC", geneV)]
  geneV = geneV[!grepl("^Mt-", geneV)]
  geneV = geneV[!grepl("^Hgb", geneV)]
  geneV = geneV[!grepl("^Hb", geneV)]
  geneV = geneV[!grepl("^Trav", geneV)]
  geneV = geneV[!grepl("^Igk", geneV)]
  
  
  
  
  if(!is.null(keepV)){
    geneV = geneV[geneV %in% keepV]
  }
  return(geneV)
}


#' Retrieve Homology Information from Ensembl
#'
#' This function queries Ensembl using the `biomaRt` package to retrieve ortholog information 
#' for multiple species relative to human genes. It returns a list where each element contains 
#' the homology data for a given species.
#'
#' @param mirror A character string specifying the Ensembl mirror to use. Default is `"https://uswest.ensembl.org"`.
#' @param species_list A character vector of species names (in Ensembl format) for which homology data should be retrieved. 
#'                     Default includes `"mmusculus"`, `"mmulatta"`, `"ppaniscus"`, `"ptroglodytes"`, `"ggorilla"`, `"nleucogenys"`, `"hsapiens"`, and `"cjacchus"`.
#'
#' @return A named list where each element contains a data frame with homology information for a given species.
#'
#' @import biomaRt dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' homology_data <- getGeneHomology()
#' head(homology_data$mmusculus)
#' }
getGeneHomology <- function(mirror = "https://uswest.ensembl.org", 
                            species_list = c(#"mmusculus", 
                                             "mmulatta"#, 
                                             # "ppaniscus", 
                                             # "ptroglodytes", 
                                             # "ggorilla", 
                                             # "nleucogenys", 
                                             #"cjacchus"
                                             )) {
  
  # Initialize HomologyLS list
  HomologyLS <- list()
  
  # Initialize Ensembl Mart
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = mirror)
  
  # Function to retrieve ortholog information for a given species
  get_orthologs <- function(species) {
    attributes <- c('ensembl_gene_id', 'chromosome_name', 'description',
                    paste0(species, '_homolog_ensembl_gene'),
                    paste0(species, '_homolog_orthology_type'),
                    paste0(species, '_homolog_orthology_confidence'))
    getBM(attributes = attributes, mart = ensembl)
  }
  
  # Retrieve gene symbols
  gene_symbols <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), mart = ensembl)
  
  # Loop through species list to get ortholog information and merge with gene symbols
  for (Species in species_list) {
    message("Processing: ", Species)
    
    # Retrieve ortholog information
    orthologs <- get_orthologs(Species)
    
    # Merge with gene symbols
    combined_data <- merge(orthologs, gene_symbols, by = "ensembl_gene_id", all.x = TRUE)
    
    # Store in list
    HomologyLS[[Species]] <- combined_data
    
    Sys.sleep(2) # Sleep for 2 seconds to avoid spamming the query system
    message(Species, " done.")
  }
  
  # Return the list
  return(HomologyLS)
}
