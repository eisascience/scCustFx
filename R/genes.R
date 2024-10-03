
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

