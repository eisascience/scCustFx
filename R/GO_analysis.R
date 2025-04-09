#' Perform GO Enrichment Analysis and Visualization
#'
#' This function performs Gene Ontology (GO) enrichment analysis using `clusterProfiler::enrichGO()`
#' and optionally generates a ggplot-based visualization of the results.
#'
#' @param gene_set A character vector of gene symbols to be analyzed.
#' @param title A character string specifying the plot title.
#' @param title_suffix A character string to be appended to the title.
#' @param gene_universe A character vector of all possible genes for enrichment analysis.
#' @param p.adjust.thr Numeric; threshold for adjusted p-values to highlight significant terms in the plot. Default is `0.05`.
#' @param returnPlot Logical; if `TRUE`, returns both the enrichment result and a ggplot visualization. Default is `TRUE`.
#' @param qrhub The OrgDb object specifying the organism database (e.g., `"org.Hs.eg.db"` for human, `"org.Mmu.eg.db"` for mouse).
#' @param labsize Numeric; text label size for GO terms in the plot. Default is `3`.
#' @param fontsize Numeric; font size scaling factor for ggplot. Default is `0.25`.
#' @param base_size Numeric; base font size for ggplot theme. Default is `18`.
#' @param max.overlaps Integer; controls the maximum number of overlapping labels in `geom_label_repel()`. Default is `20`.
#' @param MaxNsig Integer; maximum number of significant GO terms to label in the plot. Default is `30`.
#'
#' @return If `returnPlot = TRUE`, returns a list containing:
#'   \itemize{
#'     \item `ggGO` - A ggplot object visualizing the GO enrichment results.
#'     \item `ego` - The `enrichGO()` result object.
#'   }
#'   If `returnPlot = FALSE`, only the `enrichGO()` result object is returned.
#'
#' @import clusterProfiler ggplot2 ggrepel data.table
#' @export
#'
#' @examples
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' gene_set <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN") # Example genes
#' gene_universe <- unique(keys(org.Hs.eg.db, keytype = "SYMBOL"))
#' results <- enrichGOFunction(gene_set, title = "GO Enrichment", title_suffix = " - Example",
#'                             gene_universe = gene_universe, qrhub = org.Hs.eg.db)
#' print(results[[1]]) # Plot
#' }
enrichGOFunction <- function(gene_set, title, title_suffix, gene_universe, p.adjust.thr = 0.05, returnPlot=T, qrhub, 
                             labsize=3, fontsize=.25, base_size=18, max.overlaps = 20, MaxNsig=30,
                             legend.position = "bottom") {
  ego <- clusterProfiler::enrichGO(gene = gene_set,
                  universe = gene_universe,
                  OrgDb = qrhub, #"org.Hs.eg.db" "org.Mmu.eg.db
                  keyType = 'SYMBOL',
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1)
  
  ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio) / frac_to_numeric(ego@result$BgRatio)
  
  ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, "/", fixed = TRUE), function(x) {
    x <- as.numeric(x)
    x[1] / (x[2] - x[1])
  }))
  
  ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, "/", fixed = TRUE), function(x) {
    x <- as.numeric(x)
    x[1] / (x[2] - x[1])
  }))
  
  ego@result$GeneOdds[is.infinite(ego@result$GeneOdds)] <- 1
  
  
  if(returnPlot){
    ggGO <- ggplot(data.table(ego@result), aes(asinh(GeneOdds/BgOdds), -log10(pvalue), size = Count)) +
      geom_point(aes(colour = p.adjust < p.adjust.thr)) +
      scale_size_area() +
      geom_label_repel(data = data.table(ego@result)[order(p.adjust)][1:MaxNsig][p.adjust < 0.7],
                       aes(label = Description, size = fontsize), size = labsize, force = 2, max.overlaps = max.overlaps) +
      ggtitle(paste0(title, "\n", title_suffix)) +
      xlab("asinh(Odds Ratio)") +
      scale_x_log10(limits = c(1, NA), breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
      theme_classic(base_size = base_size) #+  # Automatically adjust x and y limits with extra room
      # xlim(scales::expand_range(range(ego@result$GeneOdds / ego@result$BgOdds, na.rm = TRUE), mul = 0.05)) +
      # ylim(scales::expand_range(range(-log10(ego@result$pvalue), na.rm = TRUE), mul = 0.05))
      
      # Move legend to the bottom
      theme(legend.position = legend.position )
      
      
      top_terms <- ego@result %>%
        arrange(p.adjust) %>%
        head(20) 
      
      # Create bar plot
      ggGObar <-  ggplot(top_terms, aes(x = reorder(Description, -log10(p.adjust)),
                                        y = -log10(p.adjust), fill = FoldEnrichment)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste0("Top GO Enriched Terms", "\n", title, "\n", title_suffix),
             x = "GO Terms",
             y = "-log10(Adjusted P-Value)",
             fill = "FoldEnrichment") +
        theme_minimal(base_size = 14) + theme(legend.position = legend.position )

      # 
      # ggGObubble <- ggplot(ego@result, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment), size = Count, color = p.adjust)) +
      #   geom_point(alpha = 0.7) +
      #   scale_color_gradient(low = "blue", high = "red") +
      #   labs(title = "GO Term Enrichment Bubble Plot",
      #        x = "Fold Enrichment",
      #        y = "GO Terms",
      #        size = "Gene Count",
      #        color = "Adjusted P-Value") +
      #   theme_minimal(base_size = 14)
      
      # Convert Count to factor for discrete sizes
      # ego@result$Count <- as.factor(ego@result$Count)
        # ego@result$Count <- as.numeric(ego@result$Count)
      
      # # Bubble plot with repelled gene labels
      # ggGObubble2 <- ggplot(ego@result, aes(x = FoldEnrichment, 
      #                        y = reorder(Description, FoldEnrichment), 
      #                        size = Count, 
      #                        color = p.adjust)) +
      #   geom_point(alpha = 0.7) +
      #   scale_color_gradient(low = "blue", high = "red") +
      #   scale_size_discrete(range = c(2, 10)) +  # Ensure discrete sizes
      #   geom_text_repel(aes(label = geneID), 
      #                   size = 3, 
      #                   max.overlaps = 10, 
      #                   force = 5, 
      #                   box.padding = 0.5) +  # Adjust text spacing
      #   labs(title = paste0("GO Term Enrichment Bubble Plot with Gene Labels", "\n", title, "\n", title_suffix),
      #        x = "Fold Enrichment",
      #        y = "GO Terms",
      #        size = "Gene Count",
      #        color = "Adjusted P-Value") +
      #   theme_minimal(base_size = 14) + theme(legend.position = legend.position )
      
      if (nrow(ego@result) > 0) {
        
        # Remove NA values
        ego@result <- ego@result %>% filter(!is.na(FoldEnrichment) & !is.na(Description))
        
        # Keep only top 30 enriched GO terms
        ego_top <- ego@result %>%
          arrange(p.adjust) %>%
          head(30)
        
        # Fix Count for size scaling
        if (length(unique(ego_top$Count)) == 1) {
          ego_top$Count <- as.numeric(ego_top$Count)  # Convert to numeric if all values are the same
        } else {
          ego_top$Count <- as.factor(ego_top$Count)  # Otherwise, treat as factor
        }
        
        # Bubble plot with repelled gene labels
        ggGObubble2 <- ggplot(ego_top, aes(x = FoldEnrichment, 
                                           y = reorder(Description, FoldEnrichment), 
                                           size = Count, 
                                           color = p.adjust)) +
          geom_point(alpha = 0.7) +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_discrete(range = c(2, 10)) +  # Ensure discrete sizes
          geom_text_repel(aes(label = geneID), 
                          size = 3, 
                          max.overlaps = 10, 
                          force = 5, 
                          box.padding = 0.5) +  # Adjust text spacing
          labs(title = paste0("GO Term Enrichment Bubble Plot with Gene Labels", "\n", title, "\n", title_suffix),
               x = "Fold Enrichment",
               y = "GO Terms",
               size = "Gene Count",
               color = "Adjusted P-Value") +
          theme_minimal(base_size = 14) + 
          theme(legend.position = legend.position )
        
      } else {
        ggGObubble2 <- NULL  # Prevents function from breaking if no terms are enriched
      }
      
      # library(enrichplot)
      # emapplot(ego, showCategory = 30)
      
      # library(ComplexUpset)
      
      # Convert geneID column to list format
      # ego@result$geneID_list <- strsplit(ego@result$geneID, "/")
      
      # Generate upset plot
      # ggUpset = upset(
      #   fromList(setNames(ego@result$geneID_list, ego@result$Description)),
      #   sets = ego@result$Description[1:10],
      #   order.by = "freq",
      #   mainbar.y.label = "GO Term Overlaps",
      #   sets.x.label = "Genes Per GO Term"
      # )
    
    return(list(ggGO, 
                ggGObar, 
                # ggGObubble,
                cnetplot(ego, showCategory = 10, foldChange = NULL),
                # ggUpset, 
                ggGObubble2))
    
  } else{
    return(ego)
  }
  
}


#' GO Volcano Plot
#'
#' Generates a volcano plot for Gene Ontology (GO) enrichment results from a specified component of a GO data object.
#' Highlights significant GO terms (based on adjusted p-value) and labels the top 30 enriched terms with p.adjust < 0.7.
#'
#' @param x A list-like object containing GO enrichment results for different components. Default is `GO_data`.
#' @param component A character string indicating the name of the component within `x` to plot. Default is `"V5N"`.
#' @param extraTitle An optional string to append to the plot title. If empty, it defaults to `"Component : <component>"`.
#'
#' @return A `ggplot` object displaying the volcano plot of GO terms for the specified component.
#' 
#' @import ggplot2
#' @import data.table
#' @import ggrepel
#' @export
#'
#' @examples
#' # Assuming GO_data is a properly formatted object with a "V5N" component
#' go_volcano_plot(GO_data, component = "V5N")
go_volcano_plot <- function(x=GO_data, component="V5N", extraTitle=""){
  if(extraTitle=="") extraTitle = paste("Component : ", component, sep="")
  
  #print(
  ggplot(data.table(x[[component]]), aes(GeneOdds/BgOdds, -log(pvalue), size=Count)) +
    geom_point(aes(colour=p.adjust<0.05)) +
    scale_size_area() +
    geom_label_repel(data = data.table(x[[component]])[order(p.adjust)][1:30][p.adjust<0.7], aes(label = Description, size=0.25), size = 3, force=2) + 
    ggtitle(paste("",extraTitle, sep="\n") ) +
    xlab("Odds Ratio") +
    scale_x_log10(limits=c(1,NA), breaks=c(1,2,3,4,5,6,7,8))
  #)
}


