#' Run GSEA from a DE Results Data Frame and a Gene Signature
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using a differential 
#' expression (DE) data frame and a gene signature data frame. The gene signature 
#' is split into up-regulated and down-regulated gene sets based on the signcon 
#' column. The function returns both the GSEA results object and a sorted, named 
#' vector of log2FoldChange values.
#'
#' @param de_df A data frame containing DE results with at least two columns: 
#'   \code{gene} (gene symbol) and \code{log2FoldChange} (log-fold change values).
#' @param signature_df A data frame containing the gene signature. It must include 
#'   at least the columns \code{SYMBOL} (gene symbol) and \code{signcon} (a numeric 
#'   value where positive indicates up-regulation and negative indicates down-regulation).
#' @param minGSSize An integer specifying the minimum gene set size to test in GSEA. 
#'   Default is \code{5}.
#' @param maxGSSize An integer specifying the maximum gene set size to test in GSEA. 
#'   Default is \code{500}.
#' @param pvalueCutoff A numeric value defining the p-value cutoff for GSEA significance. 
#'   Default is \code{0.05}.
#'
#' @details
#' The function splits the provided \code{signature_df} into two groups based on 
#' the \code{signcon} value: genes with \code{signcon > 0} (Signature_Up) and 
#' genes with \code{signcon < 0} (Signature_Down). It creates a TERM2GENE data frame 
#' with these two groups, then prepares a ranked gene vector from the \code{de_df} 
#' (ensuring there are no NA or duplicated entries for gene symbols). Finally, it 
#' performs GSEA using the \pkg{clusterProfiler} package.
#'
#' @return A list with two components:
#' \item{gsea_result}{An object containing the GSEA results from \code{clusterProfiler::GSEA}.}
#' \item{ranked_genes}{A named vector of log2 fold-change values sorted in decreasing order.}
#'
#' @examples
#' \dontrun{
#' # Read in your DE results and signature file
#' de_df <- read.csv("path/to/DE_results.csv")
#' signature_df <- read.csv("path/to/signature.csv")
#'
#' # Run GSEA
#' results <- run_GSEA_from_signature(de_df, signature_df)
#'
#' # Access the GSEA result and ranked genes
#' gsea_result <- results$gsea_result
#' ranked_genes <- results$ranked_genes
#'
#' # Plot the enrichment for one gene set
#' library(enrichplot)
#' gseaplot2(gsea_result, geneSetID = "Signature_Down", title = "GSEA Plot: Signature_Down")
#' }
#'
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr filter pull
#'
#' @export
run_GSEA_from_signature <- function(de_df, signature_df,
                                    minGSSize = 5,
                                    maxGSSize = 500,
                                    pvalueCutoff = 0.05) {
  # Ensure required packages are installed
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("Install clusterProfiler")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")
  
  # ---------------------------
  # Build TERM2GENE from the signature
  # ---------------------------
  # Split the signature based on signcon: positive for up-regulated, negative for down-regulated
  signature_up <- signature_df %>% dplyr::filter(signcon > 0) %>% dplyr::pull(SYMBOL)
  signature_down <- signature_df %>% dplyr::filter(signcon < 0) %>% dplyr::pull(SYMBOL)
  
  term2gene <- data.frame(
    term = rep(c("Signature_Up", "Signature_Down"), 
               times = c(length(signature_up), length(signature_down))),
    gene = c(signature_up, signature_down)
  )
  
  # ---------------------------
  # Prepare ranked genes from DE results
  # ---------------------------
  if (!all(c("gene", "log2FoldChange") %in% colnames(de_df))) {
    stop("The DE dataframe must include 'gene' and 'log2FoldChange' columns.")
  }
  
  # Remove NA values and duplicate gene entries based on gene symbol
  de_df <- de_df[!is.na(de_df$log2FoldChange), ]
  de_df <- de_df[!duplicated(de_df$gene), ]
  
  # Create a named vector with the log2FoldChange; gene names are used as names
  ranked_genes <- de_df$log2FoldChange
  names(ranked_genes) <- de_df$gene
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # ---------------------------
  # Run GSEA using clusterProfiler
  # ---------------------------
  gsea_res <- clusterProfiler::GSEA(
    geneList = ranked_genes,
    TERM2GENE = term2gene,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff
  )
  
  # Return a list containing both the GSEA results and the ranked gene vector
  return(list(gsea_result = gsea_res, ranked_genes = ranked_genes))
}



# calc_running_score <- function(ranked_genes, gene_set, weight = 1) {
#   # Create a binary vector: 1 if gene is in gene_set, 0 otherwise
#   hits <- as.integer(names(ranked_genes) %in% gene_set)
#   
#   # Total genes, and number of hits and misses
#   N <- length(ranked_genes)
#   Nh <- sum(hits)
#   Nmiss <- N - Nh
#   
#   # Check to avoid division by zero in case no gene hits
#   if (Nh == 0) {
#     stop("No genes from gene_set are in the ranked list.")
#   }
#   
#   # Compute weighted increments for hits
#   # This follows the original GSEA definition: hit increments are proportional to |ranks|^weight.
#   # Here, weight = 1 (as in the standard GSEA method).
#   hit_scores <- abs(ranked_genes) ^ weight
#   P_hit <- cumsum(hits * hit_scores) / sum(hit_scores[hits == 1])
#   
#   # For misses, each miss contributes equally.
#   P_miss <- cumsum(1 - hits) / Nmiss
#   
#   # Running enrichment score: difference between cumulative hits and misses
#   running_score <- P_hit - P_miss
#   
#   return(running_score)
# }