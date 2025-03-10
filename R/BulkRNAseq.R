#' Create a Model Matrix for Experimental Design
#'
#' This function generates a model matrix for statistical analysis based on metadata. It creates a grouping variable
#' from the specified contrast columns, ensuring proper formatting for downstream analysis.
#'
#' @param meta_df A data frame containing sample metadata.
#' @param contrast_columns Character vector. Column names in `meta_df` to be used for defining experimental groups.
#' @param sampleIdCol Character. Name of the column in `meta_df` that contains sample identifiers (default: `"cDNA_ID"`).
#'
#' @return A model matrix where rows correspond to samples and columns correspond to experimental conditions.
#'         The matrix includes attributes `"contrast_columns"` and `"sampleIdCol"` for reference.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Merges (`unites`) values from `contrast_columns` into a single `group` column.
#'   \item Ensures valid R variable names using `make.names()`.
#'   \item Constructs a design matrix using `stats::model.matrix()`, encoding the experimental conditions.
#'   \item Attaches metadata attributes (`contrast_columns`, `sampleIdCol`) to the resulting matrix.
#' }
#'
#' @importFrom dplyr mutate select across
#' @importFrom tidyr unite all_of
#' @importFrom stats model.matrix
#' @importFrom magrittr set_rownames set_colnames
#'
#'
#' @examples
#' \dontrun{
#' if(requireNamespace("dplyr", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)) {
#'     metadata <- data.frame(
#'         cDNA_ID = c("Sample1", "Sample2", "Sample3"),
#'         Treatment = c("A", "B", "A"),
#'         Batch = c("X", "X", "Y")
#'     )
#'     design_matrix <- DesignModelMatrix(metadata, contrast_columns = c("Treatment", "Batch"))
#'     print(design_matrix)
#' }
#' }
#' @export
DesignModelMatrix <- function (meta_df, contrast_columns, sampleIdCol = "cDNA_ID") 
{
  
  # contrast_columns = "TreatmentCode"
  # sampleIdCol = "MonkeyNumber"
  
  colData_intermediate <- tidyr::unite(mutate(select(as.data.frame(meta_df), 
                                                     all_of(contrast_columns)), dplyr::across(everything(), 
                                                                                              function(x) gsub(x, pattern = "_", replacement = "."))), 
                                       "group", tidyr::all_of(contrast_columns))
  colData_intermediate$group <- make.names(colData_intermediate$group)
  meta_df$group <- colData_intermediate$group
  experiment_information <- data.frame(meta_df, 
                                       row.names = NULL) %>% dplyr::select(tidyr::all_of(c(sampleIdCol, 
                                                                                           "group")))
  design <- stats::model.matrix(~0 + experiment_information$group) %>% 
    magrittr::set_rownames(experiment_information[[sampleIdCol]]) %>% 
    magrittr::set_colnames(levels(factor(experiment_information$group)))
  attr(design, "contrast_columns") <- contrast_columns
  attr(design, "sampleIdCol") <- sampleIdCol
  return(design)
}




#' Map Gene Identifiers in DESeq2 Results
#'
#' This function maps Ensembl gene IDs in DESeq2 differential expression results to target gene names
#' using a provided gene mapping data frame.
#'
#' @param DE_Seq_Result A data frame containing DESeq2 results, where row names are Ensembl gene IDs.
#' @param gene_mapping A data frame with at least two columns: `"input"` (Ensembl IDs) and `"target"` (gene symbols).
#'                      Default is `gene_mapping_MMUL10`.
#'
#' @return A data frame with an additional `"ensembl_ID"` column and a `"gene"` column containing mapped gene names.
#'         If no match is found in `gene_mapping`, the Ensembl ID is retained as the gene name.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Extracts Ensembl gene IDs from `DE_Seq_Result` row names.
#'   \item Merges DESeq2 results with `gene_mapping` to retrieve corresponding gene symbols.
#'   \item Assigns `"gene"` column values based on mapped names or retains Ensembl IDs if no match is found.
#' }
#'
#' @importFrom dplyr left_join
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example DESeq2 results (mock data)
#' DE_Seq_Result <- data.frame(log2FoldChange = c(1.5, -2.3, 0.8),
#'                             padj = c(0.01, 0.05, 0.2))
#' rownames(DE_Seq_Result) <- c("ENSG000001", "ENSG000002", "ENSG000003")
#'
#' # Example gene mapping table
#' gene_mapping <- data.frame(input = c("ENSG000001", "ENSG000002"),
#'                            target = c("GeneA", "GeneB"))
#'
#' # Run gene mapping
#' mapped_results <- map_genes(DE_Seq_Result, gene_mapping)
#' print(mapped_results)
#' }
map_genes <- function(DE_Seq_Result, gene_mapping = gene_mapping_MMUL10) {
  DE_Seq_Result$ensembl_ID <- rownames(DE_Seq_Result)
  Mapped_DE_Seq_Result <- merge(as.data.frame(DE_Seq_Result), 
                                gene_mapping[,c("input", "target")], 
                                by.x = "ensembl_ID", 
                                by.y = "input")
  Mapped_DE_Seq_Result$gene <- ifelse(is.na(Mapped_DE_Seq_Result$target), 
                                      Mapped_DE_Seq_Result$ensembl_ID, 
                                      Mapped_DE_Seq_Result$target)
  return(Mapped_DE_Seq_Result)
}
