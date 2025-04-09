
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



#'  STRING Network from Gene List
#'
#' Maps a vector of gene symbols to STRINGdb and plots the resulting protein-protein interaction network.
#'
#' @param gene_vector A character vector of gene symbols (e.g., DE genes).
#' @param species Numeric: NCBI taxonomy ID for the species (default: 9606 for Homo sapiens).
#' @param version Character: STRING database version (default: "11.5").
#' @param score_threshold Numeric: Minimum interaction confidence score (default: 400).
#' @param remove_unmapped Logical: Whether to remove unmapped genes (default: TRUE).
#' @param plot Logical: Whether to display the STRING network plot (default: TRUE).
#'
#' @return A data frame of mapped genes with STRING IDs. Optionally displays the network plot.
#' @export
#'
#' @examples
#' de_genes <- c("CD3D", "IL2", "FOXP3", "CD28")
#' result <- plot_string_network(de_genes, species = 9606)
STRINGdb_network <- function(gene_vector,
                                species = 9606,
                                version = "11.5",
                                score_threshold = 400,
                                remove_unmapped = TRUE,
                                plot = TRUE) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Please install the STRINGdb package: install.packages('STRINGdb')")
  }
  
  # Load package and initialize
  library(STRINGdb)
  string_db <- STRINGdb$new(version = version, species = species, score_threshold = score_threshold)
  
  # Prepare input and map to STRING IDs
  input_df <- data.frame(gene = gene_vector)
  mapped <- string_db$map(input_df, "gene", removeUnmappedRows = remove_unmapped)
  
  # Plot if requested
  if (plot && nrow(mapped) > 0) {
    string_db$plot_network(mapped$STRING_id)
  } else if (plot && nrow(mapped) == 0) {
    message("No genes mapped to STRING. Plot not generated.")
  }
  
  return(mapped)
}



#' Run GSEA Using fgsea
#'
#' Performs preranked GSEA and returns results for plotting.
#'
#' @param de_df A data frame with `gene` and `log2FoldChange` columns.
#' @param species Character. Species name for msigdbr (default: "Homo sapiens").
#' @param category Character. MSigDB category (default: "H").
#' @param subcategory Character. Optional subcategory (e.g., "CP:KEGG").
#' @param method Character. "multilevel" (default) or "simple".
#' @param top_n Integer. Number of top pathways to return.
#'
#' @return A data frame of GSEA results sorted by p-value.
#' @export
run_fgsea_from_de <- function(de_df,
                              species = "Homo sapiens",
                              category = "H",
                              subcategory = NULL,
                              method = "multilevel",
                              top_n = 10) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Install fgsea")
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr")
  
  if (!all(c("gene", "log2FoldChange") %in% colnames(de_df))) {
    stop("Input must include 'gene' and 'log2FoldChange' columns.")
  }
  
  # "H"
  # Hallmark
  # Most coherent & interpretable great first choice
  # "C2"
  # Curated gene sets
  # Useful; includes KEGG (CP:KEGG), Reactome (CP:REACTOME), BioCarta
  # "C7"
  # Immunologic signatures
  # Excellent for immune-focused datasets (e.g. colitis, T-cell responses)
  # "C8"
  # Cell type signatures
  # Great for deconvolution of mixed samples
  # "C5"
  # Gene Ontology (BP, MF, CC)
  # More redundant; use carefully
  # "C3"
  # Motif-based sets (TF & miRNA targets)
  # Good for upstream regulatory hypotheses
  # "C6"
  # Oncogenic signatures
  # More cancer-specific
  # "C1"
  # Chromosomal locations
  # Rarely useful for expression data
  
  # Prepare ranked genes
  de_df <- de_df[!is.na(de_df$log2FoldChange), ]
  de_df <- de_df[order(abs(de_df$log2FoldChange), decreasing = TRUE), ]
  de_df <- de_df[!duplicated(de_df$gene), ]
  ranked_genes <- de_df$log2FoldChange
  names(ranked_genes) <- de_df$gene
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Get gene sets
  m_df <- msigdbr::msigdbr(species = species, category = category)
  if (!is.null(subcategory)) {
    m_df <- subset(m_df, gs_subcat == subcategory)
  }
  gene_sets <- split(m_df$gene_symbol, m_df$gs_name)
  
  # Run fgsea
  res <- if (method == "multilevel") {
    fgsea::fgseaMultilevel(pathways = gene_sets, stats = ranked_genes)
  } else {
    fgsea::fgsea(pathways = gene_sets, stats = ranked_genes, nperm = 10000)
  }
  
  res <- res[order(res$pval), ]
  return(res)
}

#' Plot Top Enriched Pathways from fgsea Results
#'
#' @param fgsea_res Data frame output from `fgsea()`.
#' @param ranked_genes Named vector of ranked log2 fold changes.
#' @param gene_sets List of gene sets used in fgsea.
#' @param top_n Integer. Number of pathways to plot.
#'
#' @return List of ggplot enrichment plots.
#' @export
plot_top_fgsea <- function(fgsea_res, ranked_genes, gene_sets, top_n = 3, absNES  = 1.5, padj = 0.05, filterSig = T) {
  if (!requireNamespace("fgsea", quietly = TRUE)) stop("Install fgsea")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2")
  
  if(filterSig){
    # Filter significant pathways first
    fgsea_res <- fgsea_res[fgsea_res$padj < padj & abs(fgsea_res$NES) > absNES, ]
  }
  
  
  # Sort first by p-value (padj), then by absolute NES
  fgsea_res <- fgsea_res[order(fgsea_res$padj, -abs(fgsea_res$NES)), ]
  
  if (nrow(fgsea_res) == 0) {
    message("No significant enriched pathways found.")
    return(NULL)
  }
  
  # Select the top N pathways
  top_paths <- fgsea_res$pathway[1:min(top_n, nrow(fgsea_res))]
  
  # Plot each pathway
  plots <- list()
  for (path in top_paths) {
    p <- fgsea::plotEnrichment(gene_sets[[path]], ranked_genes) +
      ggplot2::ggtitle(path)
    print(p)  # Ensure the plot appears in RMarkdown
    plots[[path]] <- p
  }
  
  return(plots)
}


#' Volcano Plot for GSEA Results with Improved Labeling
#'
#' Creates a volcano plot for fgsea results, ensuring all significant pathways are labeled with minimal overlap.
#'
#' @param gsea_results Data frame output from fgsea.
#' @param ranked_genes Named vector of ranked log2 fold changes.
#' @param top_n Number of top enriched pathways to label.
#' @param pval_cutoff Adjusted p-value threshold for significance (default: 0.05).
#' @param nes_cutoff NES threshold for strong enrichment (default: 1.5).
#' @param repel Logical. Use `ggrepel` for labels (default: TRUE).
#' @param title Plot title (optional).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' gsea_plot <- plot_gsea_volcano(gsea_results, ranked_genes)
#' print(gsea_plot)
plot_gsea_volcano <- function(gsea_results,
                              # ranked_genes,
                              top_n = 20,
                              pval_cutoff = 0.05,
                              nes_cutoff = 1,
                              repel = TRUE,
                              title = "GSEA Enrichment Volcano Plot") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Please install ggrepel")
  
  # Prepare data
  gsea_results$log_padj <- -log10(gsea_results$padj)
  gsea_results$direction <- ifelse(gsea_results$NES > 0, "Up", "Down")
  
  # Compute % overlap of my genes in each pathway
  gsea_results$percent_overlap <- (gsea_results$leadingEdge %>% lengths()) / gsea_results$size * 100
  
  # Define colors
  color_palette <- c("Up" = "gold", "Down" = "dodgerblue")
  
  # Highlight significant pathways
  gsea_results$significant <- gsea_results$padj < pval_cutoff & abs(gsea_results$NES) > nes_cutoff
  
  # Select pathways to label
  labeled_paths <- gsea_results[gsea_results$significant, ]
  if (nrow(labeled_paths) > top_n) {
    labeled_paths <- labeled_paths[order(labeled_paths$padj), ][1:top_n, ]
  }
  
  # Make the plot
  p <- ggplot(gsea_results, aes(x = NES, y = log_padj, color = direction, size = percent_overlap)) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = c(-nes_cutoff, nes_cutoff), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black") +
    scale_color_manual(values = color_palette) +
    scale_size_continuous(name = "% Overlap", range = c(2, 6)) +
    labs(x = "Normalized Enrichment Score (NES)", 
         y = "-log10(Adjusted P-value)", 
         color = "Direction",
         title = title) +
    theme_minimal()
  
  # Improved labeling
  if (repel) {
    p <- p + ggrepel::geom_text_repel(data = labeled_paths, 
                                      aes(label = pathway),
                                      size = 3.5, 
                                      max.overlaps = Inf,  # Ensure all labels are shown
                                      min.segment.length = 0,  # Ensure segments are always drawn
                                      box.padding = 0.5,  # Better spacing around text
                                      point.padding = 0.3, 
                                      segment.alpha = 0.7, 
                                      max.iter = 10000)  # Avoid overlapping labels
  } else {
    p <- p + geom_text(data = labeled_paths, aes(label = pathway), size = 3, vjust = -0.8)
  }
  
  return(p)
}

#' Run Gene Set Enrichment Analysis on Loadings
#'
#' This function performs gene set enrichment analysis for both positive and negative loadings 
#' using the \code{clusterProfiler} package with gene sets obtained via \code{msigdbr}.
#' It extracts the top \code{topN} genes in both directions for a specified component 
#' and then uses \code{enricher} for enrichment analysis. If any significant enrichment is detected 
#' (based on the adjusted p-value cutoff), a dotplot is displayed.
#'
#' @param loadings A numeric matrix or data frame of gene loadings with genes as row names.
#'   Rows represent genes and columns represent different components.
#' @param comp Numeric. The component number (row index) to analyze.
#' @param topN Numeric. The number of top genes (from both positive and negative directions) to consider.
#' @param species Character. The species to query from \code{msigdbr}. Default is \code{"Homo sapiens"}.
#' @param category Character. The gene set category to use (e.g., \code{"H"} for Hallmark). Default is \code{"H"}.
#' @param pvalueCutoff Numeric. The p-value cutoff to consider enrichment as significant. Default is 0.05.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{enrichedPos}{An object of class \code{enrichResult} with enrichment results for positively loaded genes.}
#'   \item{enrichedNeg}{An object of class \code{enrichResult} with enrichment results for negatively loaded genes.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming results$loadings[[1]] is a numeric matrix with gene symbols as row names:
#'   loadings <- results$loadings[[1]]
#'   enrichmentResults <- run_bidirectional_enrichment(loadings, comp = 2, topN = 100)
#' }
#'
#' @importFrom clusterProfiler enricher dotplot
#' @import msigdbr
#' @export
run_bidirectional_enrichment <- function(loadings, comp, topN, 
                           species = "Homo sapiens", category = "H", 
                           pvalueCutoff = 0.05, doPlot = F) {
  
  # Check that loadings has row names representing gene symbols
  if(is.null(rownames(loadings))) {
    stop("loadings must have row names representing gene symbols")
  }
  
  # Load gene sets from msigdbr
  m_df <- msigdbr::msigdbr(species = species, category = category)
  
  # Extract top positively and negatively loaded genes for the specified component
  topPos <- names(sort(loadings[comp, ], decreasing = TRUE))[1:topN]
  topNeg <- names(sort(loadings[comp, ], decreasing = FALSE))[1:topN]
  
  # Enrichment analysis using clusterProfiler's enricher function
  enrichedPos <- clusterProfiler::enricher(topPos,
                                           TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
                                           pvalueCutoff = pvalueCutoff)
  
  enrichedNeg <- clusterProfiler::enricher(topNeg,
                                           TERM2GENE = m_df[, c("gs_name", "gene_symbol")],
                                           pvalueCutoff = pvalueCutoff)
  if(doPlot){
    # Display dotplots if there are significant enrichment results
    if(nrow(subset(enrichedPos@result, p.adjust < pvalueCutoff)) > 0) {
      # print(clusterProfiler::dotplot(enrichedPos, showCategory = 10, title = paste0("Comp n = ", comp, " Top Pos Genes")))
      print(clusterProfiler::dotplot(enrichedPos, showCategory = 10) + ggtitle(paste0("Comp n = ", comp, " Top Pos Genes")))
      
    }
    
    if(nrow(subset(enrichedNeg@result, p.adjust < pvalueCutoff)) > 0) {
      # print(clusterProfiler::dotplot(enrichedNeg, showCategory = 10), title = paste0("Comp n = ", comp, " Top Neg Genes"))
      print(clusterProfiler::dotplot(enrichedNeg, showCategory = 10) + ggtitle(paste0("Comp n = ", comp, " Top Neg Genes")))
    }
  } 
  
  
  # Return the enrichment objects
  return(list(enrichedPos = enrichedPos, enrichedNeg = enrichedNeg))
}




#' Run Bidirectional fgsea Analysis from Gene Loadings
#'
#' This function performs a bidirectional Gene Set Enrichment Analysis (GSEA) using the fgsea package
#' on a ranked gene loadings vector. It retrieves gene sets using msigdbr for a specified species and category,
#' executes fgsea (either using fgseaMultilevel or fgsea based on the selected method), and optionally generates
#' plots for the top enriched pathways.
#'
#' @param loadings_vector A named numeric vector of gene loadings. The names must correspond to gene symbols.
#' @param topN Integer specifying the number of top pathways to plot for each direction (default: 6).
#' @param species Character string indicating the species for which to retrieve gene sets (default: "Homo sapiens").
#' @param category Character string representing the gene set category from msigdbr (default: "H").
#' @param pvalueCutoff Numeric value for the adjusted p-value cutoff to determine significant pathways (default: 0.05).
#' @param doPlot Logical flag indicating whether to generate plots for the top pathways (default: FALSE).
#' @param nperm Integer specifying the number of permutations for the fgsea function (default: 1000).
#' @param method Character string specifying the fgsea method to use; "multilevel" for fgseaMultilevel,
#'   or any other value to use the standard fgsea function (default: "multilevel").
#'
#' @return A list containing:
#'   \describe{
#'     \item{gsea_results}{A data frame with the full results from the GSEA analysis.}
#'     \item{pos_results}{A data frame with significantly enriched pathways having positive Normalized Enrichment Scores (NES).}
#'     \item{neg_results}{A data frame with significantly enriched pathways having negative NES.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Create a named vector of gene loadings
#'   loadings <- c(GeneA = 1.2, GeneB = -0.8, GeneC = 2.3)
#'   
#'   # Run the bidirectional fgsea analysis
#'   result <- run_bidirectional_fgsea_from_loadings(loadings)
#' }
#'
#' @export
run_bidirectional_fgsea_from_loadings <- function(loadings,
                                                  comp = 1,
                                                  topN = 6, 
                                                  species = "Homo sapiens", 
                                                  category = "H",
                                                  subcategory = NULL,
                                                  pvalueCutoff = 0.05, 
                                                  doPlot = FALSE,
                                                  nperm = 1000,
                                                  method = "multilevel") {
  
  loadings_vector = loadings[comp, ]
  
  if(sd(loadings_vector) > 0){ # 0 var vectors are useless right!!?!
    # Check that loadings_vector has names representing gene symbols
    if (is.null(names(loadings_vector))) {
      stop("loadings_vector must have gene names as names")
    }
    
    if(sum(duplicated(names(loadings_vector))) > 0){
      print("removing duplicates")
      loadings_vector = loadings_vector[!duplicated(names(loadings_vector))]
      
    }
    
    # Create a ranked gene vector (sorted in decreasing order)
    ranked_genes <- sort(loadings_vector, decreasing = TRUE)
    
    # Load gene sets from msigdbr and split by pathway (gs_name)
    m_df <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory) #%>%
      # dplyr::select(gs_name, gene_symbol) %>%
      # split(x = .$gene_symbol, f = .$gs_name)
    if (!is.null(subcategory)) {
      m_df <- subset(m_df, gs_subcat == subcategory)
    }
    gene_sets <- split(m_df$gene_symbol, m_df$gs_name)
    
    # Run fgsea
    gsea_results <-  if (method == "multilevel") {
      fgsea::fgseaMultilevel(pathways = gene_sets, stats = ranked_genes)
    } else {
      # fgsea::fgsea(pathways = gene_sets, stats = ranked_genes, nperm = 10000)
      # Run fgsea with the ranked gene vector and gene sets
      fgsea::fgsea(pathways = gene_sets, 
                   stats = ranked_genes, 
                   nperm = nperm)
    }
    
    
    
    # Filter results based on the adjusted p-value cutoff
    # sig_results <- subset(gsea_results, padj < pvalueCutoff)
    
    # Split significant results into positively and negatively enriched pathways
    pos_results <- subset(gsea_results, NES > 0)
    neg_results <- subset(gsea_results, NES < 0)
    
    
    
    
    if (doPlot) {
      # Plot top pathways for positive enrichment if any exist
      if (nrow(pos_results) > 0) {
        pos_plots <- plot_top_fgsea(pos_results, ranked_genes, gene_sets, top_n = topN, filterSig = F)
        
        # Create the bottom label plot
        bottom_label <- ggdraw() + 
          draw_label(paste0("Comp n=", comp, " - top pos"), fontface = "italic", size = 14)
        
        pos_plot_grid <- cowplot::plot_grid(plotlist = pos_plots, ncol = 2)
        
        # Combine the grouped plot with the bottom label vertically
        combined_plot <- cowplot::plot_grid(pos_plot_grid, bottom_label, ncol = 1, rel_heights = c(1, 0.1))
        
        print(combined_plot)
        
        
        
      }
      
      # Plot top pathways for negative enrichment if any exist
      if (nrow(neg_results) > 0) {
        neg_plots <- plot_top_fgsea(neg_results, ranked_genes, gene_sets, top_n = topN, filterSig = F)
        
        
        # Create the bottom label plot
        bottom_label <- ggdraw() + 
          draw_label(paste0("Comp n=", comp, " - top neg"), fontface = "italic", size = 14)
        
        neg_plot_grid <- cowplot::plot_grid(plotlist = neg_plots, ncol = 2)
        
        # Combine the grouped plot with the bottom label vertically
        combined_plot <- cowplot::plot_grid(neg_plot_grid, bottom_label, ncol = 1, rel_heights = c(1, 0.1))
        
        print(combined_plot)
        
        
      }
    }
    
    # Return the full GSEA results and the split significant results
    return(list(gsea_results = gsea_results,
                pos_results = pos_results,
                neg_results = neg_results))
    
  } else {
    return(list(gsea_results = "",
                pos_results = "",
                neg_results = ""))
  }
  
}
