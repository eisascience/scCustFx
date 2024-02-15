#' @title SDA.AddCompStats
#'
#' @description SDA.AddCompStats
#' @param SDAres An SDA results list/object 
#' @return An SDA results list/object with Comp statistics
#' @import data.table
#' @export
SDA.AddCompStats <- function(SDAres, sdThrsh = 0.04, maxscoreThrsh=20, 
                             maxloadThrsh = 1, redoCalc=T){
  if (redoCalc)   SDAres$component_statistics <- NULL
  if (is.null(SDAres$component_statistics) ) {
    
    SDAres$component_statistics <- as.data.frame(data.table(
      Component = 1:SDAres$n$components, 
      Component_name = dimnames(SDAres$scores)[[2]], 
      max_score = apply(abs(SDAres$scores),  2, max),
      max_loading = apply(abs(SDAres$loadings[[1]]), 1, max),
      mean_score = apply(abs(SDAres$scores),  2, mean),
      mean_loading = apply(abs(SDAres$loadings[[1]]), 1, mean),
      sd_score = apply(abs(SDAres$scores),  2, sd),
      sd_loading = apply(abs(SDAres$loadings[[1]]), 1, sd),
      ssqrd_score = apply(SDAres$scores,  2, function(x) sum(x^2)),
      ssqrd_loading = apply(SDAres$loadings[[1]], 1, function(x) sum(x^2))
    )[order(-Component)])
    
  }
  
  SDAres$component_statistics$Component_namev2    <- SDAres$component_statistics$Component_name
  SDAres$component_statistics$Component_name_plot <- SDAres$component_statistics$Component_name
  
  SDAres$component_statistics[which(SDAres$component_statistics$sd_loading<sdThrsh),]$Component_namev2         <- rep("", length(which(SDAres$component_statistics$sd_loading<sdThrsh)))
  
  SDAres$component_statistics[which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh),]$Component_name_plot <- rep("", length(which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh)))
  
  SDAres$component_statistics <- data.table(SDAres$component_statistics)
  return(SDAres)
  
}



#' Perform Gene Ontology Enrichment Analysis on SDA Components
#'
#' This function conducts Gene Ontology (GO) enrichment analysis on specified components
#' from Single-cell Data Analysis (SDA) results. It supports analysis in positive, negative,
#' or both directions and allows specifying the number of top genes to consider. The function
#' utilizes `clusterProfiler` for the enrichment analysis and is capable of handling different
#' organism databases.
#'
#' @param sdaResults A list containing SDA results, expected to include a 'loadings' element
#'   with components and gene loadings.
#' @param components A numeric vector indicating which components from the SDA results to analyze.
#' @param orgDb A character string specifying the Bioconductor organism annotation package
#'   to use for the enrichment analysis. Default is "org.Hs.eg.db" for human. Other examples
#'   include "org.Mm.eg.db" for mouse and "org.Mmu.eg.db" for rhesus macaque.
#' @param geneNumber An integer specifying the number of top genes to consider for enrichment
#'   analysis in each component and direction. Default is 100.
#' @param direction A character string indicating the direction of gene loading to consider
#'   for the analysis. Can be "Pos" for positive, "Neg" for negative, or "Both" for both directions.
#'   Default is "Both".
#'
#' @return A list where each element corresponds to a component and direction analyzed, containing
#'   the results of the GO enrichment analysis. Results include enrichment scores, gene ratios,
#'   background ratios, and odds ratios among others.
#'   
#' @export
SDA.GO_Enrichment <- function (sdaResults, components, 
                               orgDb = "org.Hs.eg.db", 
                               # mouse org.Mm.eg.db
                               # rhesus macaque org.Mmu.eg.db
                               geneNumber = 100, 
          direction = "Both") 
{
  ret <- list()
  if (direction == "Both") {
    directions <- c("Pos", "Neg")
  }
  else {
    directions <- direction
  }
  for (comp in components) {
    for (direction in directions) {
      print(paste0("Loading component ", comp, ", direction: ", 
                   direction))
      if (direction == "Neg") {
        top_genes <- names(sort(sdaResults$loadings[[1]][comp, 
        ]))[1:geneNumber]
      }
      else if (direction == "Pos") {
        top_genes <- names(sort(sdaResults$loadings[[1]][comp, 
        ], decreasing = T))[1:geneNumber]
      }
      gene_universe <- names(sdaResults$loadings[[1]][comp, 
      ])
      ego <- clusterProfiler::enrichGO(gene = top_genes, 
                                       universe = gene_universe, OrgDb = orgDb, keyType = "SYMBOL", 
                                       ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1, 
                                       qvalueCutoff = 1)
      frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text = x)))
      ego@result$Enrichment <- frac_to_numeric(ego@result$GeneRatio)/frac_to_numeric(ego@result$BgRatio)
      ego@result$GeneOdds <- unlist(lapply(strsplit(ego@result$GeneRatio, 
                                                    "/", fixed = TRUE), function(x) {
                                                      x <- as.numeric(x)
                                                      x[1]/(x[2] - x[1])
                                                    }))
      ego@result$BgOdds <- unlist(lapply(strsplit(ego@result$BgRatio, 
                                                  "/", fixed = TRUE), function(x) {
                                                    x <- as.numeric(x)
                                                    x[1]/(x[2] - x[1])
                                                  }))
      ret[[paste0(comp, "-", direction)]] <- ego@result
    }
  }
  return(ret)
}