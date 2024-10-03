

#' This function generates a plot of gene loadings along genomic coordinates based on a Seurat object.
#'
#' @param SDARedDataLS A list of SDA data
#' @param reduction the name of the SDA run
#' @param dimN the reduction 
#' @param highlight_genes A character vector of gene symbols to highlight in the plot (default is NULL).
#' @param TopNpos The number of top positive loadings to display (default is 10).
#' @param TopNneg The number of top negative loadings to display (default is 10).
#' @param data_set passed to SDAtools:::get.location data_set default "hsapiens_gene_ensembl" can be mmulatta_gene_ensembl or mmusculus_gene_ensembl .. see ??useMart or do: mart = useMart('ensembl'), followed by listDatasets(mart).
#' @param invertWeights to invert the loading weight i.e., -loadings
#' @param includeUnMapped to include un mapped genes e.g. LOCs default = T
#' 
#' 
#' @return A ggplot object representing the loadings along genomic coordinates.
#' @export
plot_loadings_coordinates <- function(SDARedDataLS,
                                      reduction,
                                      mart = NULL,
                                      genes = NULL,
                                      dimN,
                                      highlight_genes = NULL, 
                                      TopNpos = 10, TopNneg=10,
                                      data_set = "hsapiens_gene_ensembl", #mmulatta_gene_ensembl
                                      invertWeights=F, includeUnMapped = T, geneLocPath=NULL ) {
  
  library(biomaRt)
  library(ggplot2)
  library(ggrepel)
  
  # Use the Ensembl Mart for human genes
  if(is.null(mart)) mart <- useMart("ensembl", data_set)
  
  # Get gene coordinates
  if(is.null(genes)) genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position"), mart = mart)
  
  # Calculate chromosome lengths
  chromosomes <- unique(genes$chromosome_name)
  
  chromosome_length_red <- sapply(chromosomes, function(chr) {
    chr_genes <- subset(genes, chromosome_name == chr)
    max_end <- max(chr_genes$end_position)
    return(max_end)
  })
  
  chromosome_length_red = chromosome_length_red[names(chromosome_length_red) %in% c(1:22, "X", "Y", "MT")]
  chromosome_lengthsDF = data.frame(chromosome = names(chromosome_length_red),
                                    length = as.numeric(chromosome_length_red))
  
  
  
  Genes2Map = colnames(SDARedDataLS$loadings[[reduction]]$loadings)
  
  
  
  if(!is.null(geneLocPath)) {
    
    if(file.exists(geneLocPath)){
      gene_locations = readRDS(geneLocPath)
      
    } 
    if(!file.exists(geneLocPath)){
      print("file does not exist, downloading new results")
      
      gene_locations <- SDAtools:::get.location(
        gene.symbols = Genes2Map,
        data_set = data_set,
        gene_name = "external_gene_name"
      )
      
      saveRDS(gene_locations, geneLocPath)
    }
  } else{
    # if(is.null(geneLocPath)){
    gene_locations <- SDAtools:::get.location(
      gene.symbols = Genes2Map,
      data_set = data_set,
      gene_name = "external_gene_name"
    )
    # } 
  }
  
  
  gene_locations[is.na(gene_locations$start_position), ]$chromosome_name = "?"
  gene_locations[is.na(gene_locations$start_position), ]$start_position = 1
  
  cl = lapply(unique(gene_locations$chromosome_name), function(xN){
    data.frame(xN, round( (max(subset(gene_locations, chromosome_name %in% xN)$start_position) - min(subset(gene_locations, chromosome_name %in% xN)$start_position))/2))
  }) %>% data.table::rbindlist() %>% as.data.frame()
  colnames(cl) = colnames=c("chr", "center")
  
  if(includeUnMapped){
    if(sum(grepl("^LOC", Genes2Map)) > 0){
      LOCgenes = Genes2Map[grepl("^LOC", Genes2Map)]
      
      geneLocPath_fix = gsub(".rds", "_LOCfix.rds", geneLocPath)
      
    
      
      gene_locations = rbind(gene_locations,
                             data.frame(gene_symbol = LOCgenes, 
                                        chromosome_name = rep("?", length(LOCgenes)),
                                        start_position = rep(1, length(LOCgenes))))
      
    }
  }
  
  
  temp <- merge(gene_locations, chromosome_lengthsDF, by.x = "chromosome_name", by.y = "chromosome", all.x = TRUE) %>% 
    as.data.frame()
  temp$genomic_position <- temp$start_position + temp$length / 2
  
  
  component = SDARedDataLS$loadings[[reduction]]$loadings[dimN,]
  
  if(invertWeights){
    component = component * -1
  }
  
  # Subset genes with weights
  label_data <- data.frame(
    gene_symbol = names(component),
    loading = component
  )
  
  # Merge with genomic positions
  label_data <- merge(label_data, temp, by = "gene_symbol", all.x = TRUE)
  
  if(sum(is.na(label_data$start_position))>0) {
    
    label_data[is.na(label_data$start_position), ]$chromosome_name = "?"
    label_data[is.na(label_data$start_position), ]$length = 3
      
    if(includeUnMapped & sum(is.na(label_data$start_position))>0) label_data[is.na(label_data$start_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
    label_data[is.na(label_data$start_position), ]$start_position = 1
  }
  
  if(sum(is.na(label_data$length))>0) label_data[is.na(label_data$length), ]$length = 3
  if(includeUnMapped & sum(is.na(label_data$genomic_position))>0) label_data[is.na(label_data$genomic_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
  
  
  label_data$gene_symbol_show = ""
  label_data$gene_symbol_show[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] = label_data$gene_symbol[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] 
  label_data$gene_symbol_show[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] = label_data$gene_symbol[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] 
  
  
  
  
  P <- ggplot(label_data, aes(genomic_position, loading, size = abs(loading)^2)) + 
    geom_point(stroke = 0, aes(alpha = (abs(loading))^0.7, 
                               color = ifelse(abs(loading)>.05, "cornflowerblue", "grey") ) ) + 
    scale_color_manual(values =  c("cornflowerblue", "grey"))  + 
    # scale_x_continuous(breaks = cl$center, 
    #                    labels = cl$chr, minor_breaks = NULL) +
    # scale_color_manual(values = rep_len("black", length(unique(label_data$chromosome_name))+1))+
    # scale_colour_manual(values = c(rep_len(c("black", "cornflowerblue"), 
    #                                        length(unique(label_data$chromosome_name))), "grey")) +
    xlab("Genomic Coordinate") + ylab("Weight")  +
    geom_label_repel(aes(label = gene_symbol_show), max.overlaps = max(c(TopNneg,TopNpos ))*2+1,
                     size = 3, box.padding = unit(0.5, "lines"), 
                     point.padding = unit(0.1, "lines"), force = 10, segment.size = 0.2, segment.color = "blue") +
    theme_minimal() + theme(legend.position = "none"); P
  
  # Highlight genes if necessary
  if (!is.null(highlight_genes)) {
    P <- P + geom_point(data = label_data[label_data$gene_symbol %in% highlight_genes, ], color = "red")
  }
  
  return(P)
}


#' This function generates a Heatmap of thresholding SDA score
#'
#' @param SDAscore vector of SDA score
#' @param Meta a meta vector or NULL 
#' @param GT if T greater than 
#' @param CutThresh Cut the score based on threshold .
#' 
#' 
#' @return A ggplot object representing the loadings along genomic coordinates.
#' @export
HeatMap_SDAScore_Thr <- function(SDAscore, Meta = NULL,  GT = T, CutThresh = 0, clustering_method = "ward.D2"){
  
  if(is.null(Meta)) Meta = rep(0, length(SDAscore))
  
  tempDF = data.frame(SDAscore = SDAscore, 
                      Meta = Meta)
  
  
  
  
  tempDF$GtThr = "BelowThreshold"
  
  
  if(GT) tempDF$GtThr[tempDF[,"SDAscore"]>=CutThresh] = "AboveThreshold"
  if(!GT) tempDF$GtThr[tempDF[,"SDAscore"]<CutThresh] = "AboveThreshold"
  
  
  pheatmap::pheatmap(asinh(chisq.test(table(tempDF$Meta, tempDF$GtThr))$res), clustering_method = clustering_method)
  
  
  
}