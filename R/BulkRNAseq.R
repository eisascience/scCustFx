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




#' Draw a volcano plot for one contrast derived from DEseq
#' @param contrast_df Data frame with columns: log2FoldChange, padj, gene, contrast_name
#' @param Title Title text or NULL
#' @return ggplot object
#' @export
draw_volcano_DEseq <- function(contrast_df, Title=NULL) {
  
  if(is.null(Title)) Title = unique(contrast_df$contrast_name)
  
  de_df <- contrast_df %>%
    mutate(Differential_Expression = case_when(
      log2FoldChange >= 0.58 & padj <= 0.05  ~ "Upregulated",
      log2FoldChange <= -0.58 & padj <= 0.05 ~ "Downregulated",
      TRUE                                   ~ "Not DE"
    )) %>%
    filter(!is.na(log2FoldChange), !is.na(padj), !is.na(gene))
  
  ggplot(de_df, aes(log2FoldChange, -log10(padj), color = Differential_Expression)) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.58, 0.58),     linetype = "dashed") +
    geom_label_repel(data = filter(de_df, Differential_Expression!="Not DE"),
                     aes(label = gene),
                     size = 3, box.padding = 0.5,
                     show.legend = FALSE,
                     seed = 1234, max.overlaps = 25) +
    xlim(c(-max(abs(de_df$log2FoldChange)), max(abs(de_df$log2FoldChange)))) +
    scale_color_manual(values = c(
      "Upregulated"   = "orange",
      "Downregulated" = "cadetblue2",
      "Not DE"        = "black"
    )) +
    egg::theme_article() +
    ggtitle(unique(contrast_df$contrast_name))
}

#' Plot a heatmap of differentially expressed genes for a single contrast
#'
#' This function selects DEGs from a contrast result, optionally
#' applies voom transformation and batch correction, then displays a
#' z‑scored heatmap of expression across samples.
#'
#' @param contrast_df A data.frame of DE results with columns:
#'   \describe{
#'     \item{ensembl_ID}{gene IDs matching rownames of \code{count_mat}}
#'     \item{log2FoldChange}{numeric log2 fold changes}
#'     \item{padj}{numeric adjusted p‑values}
#'     \item{contrast_name}{character string used for the plot title}
#'   }
#' @param count_mat Numeric matrix of raw counts (rows=genes, cols=samples).
#'   Row names must match \code{contrast_df$ensembl_ID}, column names match rownames of \code{meta_df}.
#' @param meta_df Data.frame of sample metadata; rows must align with cols of \code{count_mat}.
#'   Must include the column named by \code{BatchFeat} for batch removal and optionally \code{PD1status}, \code{seqrun} for annotations.
#' @param cleanBatch Logical; if \code{TRUE}, calls \code{limma::removeBatchEffect} on \code{BatchFeat}. Default \code{TRUE}.
#' @param VoomTrans Logical; if \code{TRUE}, applies \code{limma::voom} before batch correction. Default \code{FALSE}.
#' @param BatchFeat Character; the column name in \code{meta_df} indicating batch. Required when \code{cleanBatch = TRUE}.
#' @param gene_meta Data.frame mapping \code{ensembl_ID} (rownames) to final gene names. Must contain column \code{final_GeneName}.
#' @param filter Logical; if \code{TRUE}, DEGs are filtered by \code{padj} and \code{absLog2FC}. Default \code{FALSE}.
#' @param padj_thr Numeric cutoff for adjusted p‑value when \code{filter=TRUE}. Default \code{0.05}.
#' @param absLog2FC_thr Numeric cutoff for absolute log2 fold‑change when \code{filter=TRUE}. Default \code{0.58}.
#' @return Invisibly returns the \code{pheatmap} object.
#' @importFrom limma voom removeBatchEffect
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_heatmap_DEseq <- function(contrast_df,
                         count_mat,
                         meta_df,
                         cleanBatch = TRUE,
                         VoomTrans  = FALSE,
                         BatchFeat  = NULL,
                         gene_meta,
                         filter     = FALSE,
                         padj_thr       = 0.05,
                         absLog2FC_thr  = 0.58) {
  # — align count matrix to metadata —
  if (!all(rownames(meta_df) %in% colnames(count_mat))) {
    stop("Column names of count_mat must include all rownames(meta_df).")
  }
  count_mat <- count_mat[, rownames(meta_df), drop = FALSE]
  
  # — select DEGs —
  if (filter) {
    DEGs <- contrast_df %>%
      filter(!is.na(ensembl_ID),
             padj <= padj_thr,
             abs(log2FoldChange) >= absLog2FC_thr) %>%
      pull(ensembl_ID)
  } else {
    DEGs <- contrast_df$ensembl_ID
  }
  if (length(DEGs) < 3) {
    message("Too few DEGs for heatmap: ", unique(contrast_df$contrast_name))
    return(invisible(NULL))
  }
  
  # — optional voom transform —
  if (VoomTrans) {
    mat_trans <- voom(count_mat, plot = FALSE)$E
  } else {
    mat_trans <- count_mat
  }
  
  # — optional batch correction —
  if (cleanBatch) {
    if (is.null(BatchFeat)) 
      stop("`BatchFeat` must be provided when `cleanBatch = TRUE`")
    mat_trans <- removeBatchEffect(mat_trans, batch = meta_df[[BatchFeat]])
  }
  
  # — subset to DEGs and z‑score —
  mat <- mat_trans[DEGs, , drop = FALSE]
  rownames(mat) <- gene_meta[rownames(mat), "final_GeneName"]
  mat_z <- t(scale(t(mat)))
  
  # — annotation setup —
  show_rows <- length(DEGs) <= 100
  ann_cols <- meta_df[, c("PD1status", "seqrun"), drop = FALSE]
  ann_colors <- list(
    PD1status = c(IgG = "#1b9e77", `anti-PD1` = "#d95f02"),
    seqrun     = setNames(
      brewer.pal(length(unique(meta_df$seqrun)), "Set2"),
      unique(meta_df$seqrun)
    )
  )
  
  # — draw heatmap —
  pheatmap::pheatmap(
    mat_z,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    show_rownames     = show_rows,
    show_colnames     = FALSE,
    annotation_col    = ann_cols,
    annotation_colors = ann_colors,
    main              = unique(contrast_df$contrast_name)
  )
}