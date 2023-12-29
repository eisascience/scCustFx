#' Enrich GO Function
#'
#' This function performs Gene Ontology (GO) enrichment analysis and generates a
#' ggplot object for visualization.
#'
#' @param gene_set A vector of gene symbols for the enrichment analysis.
#' @param title The title to be added to the plot title.
#' @param title_suffix The suffix to be added to the plot title.
#' @param gene_universe The vector of gene symbols as the universe.
#'
#' @return A ggplot object displaying the GO enrichment analysis results.
#'
#' @export
enrichGOFunction <- function(gene_set, title, title_suffix, gene_universe, p.adjust.thr = 0.05, returnPlot=T, qrhub, 
                             labsize=3, fontsize=.25, base_size=18, max.overlaps = 20, MaxNsig=30) {
  ego <- enrichGO(gene = gene_set,
                  universe = gene_universe,
                  OrgDb = qrhub,
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
  
  if(returnPlot){
    ggGO <- ggplot(data.table(ego@result), aes(GeneOdds/BgOdds, -log(pvalue), size = Count)) +
      geom_point(aes(colour = p.adjust < p.adjust.thr)) +
      scale_size_area() +
      geom_label_repel(data = data.table(ego@result)[order(p.adjust)][1:MaxNsig][p.adjust < 0.7],
                       aes(label = Description, size = fontsize), size = labsize, force = 2, max.overlaps = max.overlaps) +
      ggtitle(paste0(title, title_suffix)) +
      xlab("Odds Ratio") +
      scale_x_log10(limits = c(1, NA), breaks = c(1, 2, 3, 4, 5, 6, 7, 8)) +
      theme_classic(base_size = base_size)
    
    return(ggGO)
  } else{
    return(ego)
  }
  
}


