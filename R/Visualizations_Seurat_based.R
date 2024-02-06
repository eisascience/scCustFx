#' Generate a boxplot with jittered points for the percentage of values exceeding a threshold.
#'
#' @param SerObj The input data or object.
#' @param FeatName The name of the feature variable.
#' @param CutThresh The threshold value for categorizing values.
#' @param MetaDataName The name of the metadata variable for grouping (default: "Population").
#' @param TitleExtra Additional title text for the plot.
#' @param col The color for the boxplot and jittered points (default: "#A6CEE3").
#' @param MetaDataName2 The name of the second metadata variable (default: "SubjectId").
#' @param ylab_text The label for the y-axis.
#' @param xlab_text The label for the x-axis.
#' @param base_size The base font size for the plot (default: 16).
#' @param my_comparisons A vector specifying pairwise comparisons for statistical testing (default: NULL).
#' @param limits A vector to define the limits of the percent axis (default: c(0, 110)).
#' @param GT T/F to select greater than and equal to threshold if T else less than (default: T).
#' @param StatMethod passed to stat_compare_means only in t.test and in wilcox.test (default: wilcox.test).
#' @param StatLab passed to stat_compare_means  p.signif" (shows the significance levels), "p.format" (shows the formatted p value), (default: p.signif).
#'
#'
#' @return A ggplot object representing the boxplot with jittered points.
#'
#' @examples
#' # Example usage:
#' FeatThrTab_boxplot(SerObj, FeatName = "Feature1", CutThresh = 0.5)
#'
#'
#' @export
FeatThrTab_boxplot <- function(SerObj, 
                               FeatName = NULL, CutThresh = NULL,
                               MetaDataName = "Population", TitleExtra = "", 
                               col = "#A6CEE3", 
                               MetaDataName2 = "SubjectId",
                               ylab_text = "", xlab_text = "",
                               base_size = 16, 
                               my_comparisons = NULL,
                               limits = c(0, 110), GT=T, StatMethod = "wilcox.test",
                               StatLab = "p.signif"){
  
  if(is.null(FeatName)) stop("FeatName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(FeatName) != 1) stop ("FeatName Length != 1")
  
  
  tempDF <- FetchData(SerObj, vars = c(FeatName, MetaDataName, MetaDataName2))
  # colnames(tempDF) <- c("var1", "var2", MetaDataName2, MetaDataName3)
  
  tempDF$GtThr = "FALSE"
  
  if(GT) tempDF$GtThr[tempDF[,FeatName]>=CutThresh] = "TRUE"
  if(!GT) tempDF$GtThr[tempDF[,FeatName]<CutThresh] = "TRUE"
  
  tempDF$percent_gt = tempDF %>%
    group_by(!!sym(MetaDataName), !!sym(MetaDataName2)) %>%
    mutate(num_gt_thr = sum(GtThr == "TRUE"),
           percent_gt = num_gt_thr / n()) %>%
    ungroup() %>% pull(percent_gt) #%>% hist(breaks = 100)
  
  
  
  
  subj_tissue_data <- aggregate(as.formula(paste("percent_gt ~", MetaDataName, "+", MetaDataName2)), 
                                data = tempDF, FUN = mean)
  
  subj_tissue_data$percent_gt = subj_tissue_data$percent_gt *100
  
  # print(table(subj_tissue_data[,MetaDataName]))
  LevLeng = length(levels(factor(as.character(subj_tissue_data[,MetaDataName]))))
  
  library(ggpubr)
  
  
  gg2 = ggboxplot(subj_tissue_data, x = MetaDataName, y = "percent_gt", 
                  palette = c(rep(col, LevLeng)),
                  # add = "jitter", 
                  fill=MetaDataName, outlier.shape = NA) +
    geom_jitter(alpha=.5, size=.8) +
    theme_classic(base_size = base_size) +
    ggtitle(TitleExtra) + 
    xlab(xlab_text) + ylab(ylab_text)  +
    theme(legend.position="none",
          # legend.direction="horizontal",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    scale_y_continuous(limits = limits);
  
  if(!is.null(my_comparisons)){
    gg2 = gg2 +
      stat_compare_means(comparisons = my_comparisons, method = StatMethod , label = StatLab)
  }
  
  
  
  return(gg2)
  
  
}

#' Plot violin plots with Wilcoxon p-values
#'
#' This function generates a plot of violin plots with Wilcoxon p-values
#' for the given feature and groupings in a Seurat object.
#'
#' @param SerObj A Seurat object containing the data.
#' @param feature The feature to plot.
#' @param group.by The grouping variable to use for the plot (default: "ident").
#' @param layer The layer to use for the plot (default: "data").
#' @param GrpFacLevels Optional vector of grouping factor levels to plot.
#' @param assay The assay to use for the plot (default: "RNA").
#' @param colors A vector of colors to use for the plot.
#' @param comparisons A list of pairwise comparisons to perform.
#' @param title A string to use as the plot title.
#' @param doJitter Whether to add jitter to the points (default: TRUE).
#'
#' @return A ggplot object.
#' @export
plot_violin_wpvalue <- function(SerObj, 
                                feature, 
                                group.by = "ident",
                                layer="data", 
                                GrpFacLevels=NULL,
                                assay="RNA", 
                                colors =col_vector, 
                                comparisons = NULL, title="",
                                doJitter = T){
  
  
  
  
  library(ggpubr)
  
  v1 <- VlnPlot(SerObj, features =feature, assay = assay, group.by=group.by)
  
  v1 <- v1$data
  colnames(v1) = c("Feat", "Grp")
  
  if (!is.null(GrpFacLevels)) {
    v1$Grp <- factor(as.character(v1$Grp), levels = GrpFacLevels)
  }
  
  
  if(!is.null(comparisons)){
    
    if(class(comparisons) != "list") stop("GrpNames needs to be a list")
    
    if(doJitter){
      ggviolin(v1, x = "Grp", y = "Feat", 
               fill = "Grp",
               palette=col_vector[1:length(levels(factor(v1$Grp)))],
               add = "jitter",
               add.params = list(size = .005, alpha = 0.05),
               title = paste0("Wilcox p.value. :: ", title)) + 
        stat_compare_means(comparisons = comparisons,
                           label = "p.signif")
    } else {
      ggviolin(v1, x = "Grp", y = "Feat", fill = "Grp",
               palette=col_vector[1:length(levels(factor(v1$Grp)))],
               title = paste0("Wilcox p.value. :: ", title)) + 
        stat_compare_means(comparisons = comparisons , label = "p.signif")
    }
    
  } else {
    if(doJitter){
      
      ggviolin(v1, x = "Grp", y = "Feat", 
               fill = "Grp",
               palette=col_vector[1:length(levels(factor(v1$Grp)))],
               add = "jitter",
               title = paste0("Wilcox p.value. :: ", title)) 
    } else {
      ggviolin(v1, x = "Grp", y = "Feat", fill = "Grp",
               palette=col_vector[1:length(levels(factor(v1$Grp)))],
               title = paste0("Wilcox p.value. :: ", title)) 
    }
  }
  
  
  
}



#' This function generates a plot of gene loadings along genomic coordinates based on a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param reduction The dimensionality reduction method used in Seurat (default is "pca").
#' @param dimN the reduction dim default = 1 
#' @param highlight_genes A character vector of gene symbols to highlight in the plot (default is NULL).
#' @param TopNpos The number of top positive loadings to display (default is 10).
#' @param TopNneg The number of top negative loadings to display (default is 10).
#' @param data_set passed to SDAtools:::get.location data_set default "hsapiens_gene_ensembl" can be  #mmulatta_gene_ensembl .. see ??useMart or do: mart = useMart('ensembl'), followed by listDatasets(mart).
#' @param invertWeights to invert the loading weight i.e., -loadings
#' @param includeUnMapped to include un mapped genes e.g. LOCs default = T
#' 
#' 
#' @return A ggplot object representing the loadings along genomic coordinates.
#' @export
plot_loadings_coordinates <- function(SerObj, reduction = "pca", redLab = "PC",
                                      highlight_genes = NULL, 
                                      TopNpos = 10, TopNneg=10,
                                      data_set = "hsapiens_gene_ensembl", #mmulatta_gene_ensembl
                                      dimN = 1, invertWeights=F, includeUnMapped = T, geneLocPath=NULL ) {
  
  library(biomaRt)
  library(ggplot2)
  library(ggrepel)
  
  # Use the Ensembl Mart for human genes
  mart <- useMart("ensembl", "hsapiens_gene_ensembl")
  
  # Get gene coordinates
  genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position"), mart = mart)
  
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
  

    Genes2Map = rownames(SerObj)


  
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
  
  if(includeUnMapped){
    if(sum(grepl("^LOC", Genes2Map)) > 0){
      LOCgenes = Genes2Map[grepl("^LOC", Genes2Map)]
      
      geneLocPath_fix = gsub(".rds", "_LOCfix.rds", geneLocPath)
      
      # if(!is.null(geneLocPath)) {
      #   
      #   if(file.exists(geneLocPath_fix)){
      #     LOCgenesFix = readRDS(geneLocPath_fix)
      #     
      #   } else {
      #     LOCgenesFix = CellMembrane:::ResolveLocGenes(LOCgenes)
      #     saveRDS(LOCgenesFix, geneLocPath_fix)
      #   }
      #   
      # }
      
      
      
      # LOCgenes %in% gene_locations$gene_symbol
      
      gene_locations = rbind(gene_locations,
                             data.frame(gene_symbol = LOCgenes, 
                                        chromosome_name = rep("?", length(LOCgenes)),
                                        start_position = rep(1, length(LOCgenes))))
      
    }
  }
  
  
  temp <- merge(gene_locations, chromosome_lengthsDF, by.x = "chromosome_name", by.y = "chromosome", all.x = TRUE) %>% 
    as.data.frame()
  temp$genomic_position <- temp$start_position + temp$length / 2
  
  # component = Seurat::Loadings(SerObj, reduction = reduction)[, paste0(redLab, "_", dimN)]
  component = Seurat::Loadings(SerObj, reduction = reduction)[, dimN]
  
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
  
 
  label_data[is.na(label_data$start_position), ]$chromosome_name = "?"
  label_data[is.na(label_data$start_position), ]$length = 3
  if(includeUnMapped) label_data[is.na(label_data$start_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
  label_data[is.na(label_data$start_position), ]$start_position = 1
  
  label_data[is.na(label_data$length), ]$length = 3
  if(includeUnMapped) label_data[is.na(label_data$genomic_position), ]$genomic_position = max(label_data$genomic_position, na.rm = T) + 2
  
  
  label_data$gene_symbol_show = ""
  label_data$gene_symbol_show[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] = label_data$gene_symbol[order(label_data$loading, decreasing = TRUE)[1:TopNpos]] 
  label_data$gene_symbol_show[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] = label_data$gene_symbol[order(label_data$loading, decreasing = FALSE)[1:TopNneg]] 
  
  
  
  
  P <- ggplot(label_data, aes(genomic_position, loading, size = abs(loading)^2)) + 
    geom_point(stroke = 0, aes(alpha = (abs(loading))^0.7, color = ifelse(abs(loading)>.05, "cornflowerblue", "grey"))) + 
    scale_color_manual(values =  c("cornflowerblue", "grey"))+
    # scale_color_manual(values = rep_len("black", length(unique(label_data$chromosome_name))+1))+
    # scale_colour_manual(values = c(rep_len(c("black", "cornflowerblue"), length(unique(label_data$chromosome_name))), "grey")) +
    xlab("Genomic Coordinate") + ylab("Weight")  +
    geom_label_repel(aes(label = gene_symbol_show), max.overlaps = max(c(TopNneg,TopNpos ))*2+1,
                     size = 3, box.padding = unit(0.5, "lines"), 
                     point.padding = unit(0.1, "lines"), force = 10, segment.size = 0.2, segment.color = "blue") +
    theme_minimal() + theme(legend.position = "none")
  
  # Highlight genes if necessary
  if (!is.null(highlight_genes)) {
    P <- P + geom_point(data = label_data[label_data$gene_symbol %in% highlight_genes, ], color = "red")
  }
  
  return(P)
}






#' Analyze PC Enrichment
#'
#' This function performs enrichment analysis for the genes associated with positive
#' and negative loadings of a given principal component (PC) in a Seurat object.
#'
#' @param SerObj A Seurat object containing the data.
#' @param CompN The principal component number.
#' @param topN The number of top genes to consider for enrichment analysis.
#' @param qrhub The gene database used for enrichment analysis.
#'
#' @return A list containing two ggplot objects for positive and negative loadings.
#'
#' @examples
#' \dontrun{
#' result <- analyzeEnrichment(SerObj, CompN = 1, topN = 20, qrhub = "org.Hs.eg.db")
#' result$ggGO_pos | result$ggGO_neg
#' }
analyzeEnrichment <- function(SerObj, CompN = 1, topN = 20, qrhub = "org.Hs.eg.db", invertWeights=F,
                                labsize=3.5, fontsize=.25, base_size = 20, max.overlaps = 20, MaxNsig=30,
                                reduction = "pca", redLab = "PC",
                              topNGenes = 100, bottomNGenes = 100) {
  
  library(AnnotationHub)
  library(clusterProfiler)
  library(data.table)
  
  # hub <- AnnotationHub()
  
  # gene_universe <- rownames(Seurat::Loadings(SerObj, reduction = "pca"))
  gene_universe = rownames(SerObj)
  PCloadingDF <- RIRA:::ExtractGeneWeights(seuratObj = SerObj, componentNum = CompN,
                                           topNGenes = topNGenes, bottomNGenes = bottomNGenes,
                                           reduction = reduction)
  
  if(invertWeights){
    PCloadingDF$weight = PCloadingDF$weight * -1
  }
  
  NegLoadingsDF <- subset(PCloadingDF, weight < 0)
  PosLoadingsDF <- subset(PCloadingDF, weight > 0)
  
  NegLoadingsDF <- NegLoadingsDF[order(abs(NegLoadingsDF$weight), decreasing = TRUE),]
  PosLoadingsDF <- PosLoadingsDF[order(abs(PosLoadingsDF$weight), decreasing = TRUE),]
  
  top_genes_neg <- NegLoadingsDF$feature[1:min(topN, nrow(NegLoadingsDF))]
  top_genes_pos <- PosLoadingsDF$feature[1:min(topN, nrow(PosLoadingsDF))]
  
  # loc_genes <- CellMembrane::ResolveLocGenes(gene_universe[grep("^LOC", gene_universe)])
  
  
  
  ggGO_pos <- enrichGOFunction(top_genes_pos, title_suffix =  " Pos", gene_universe=gene_universe, title = paste0(redLab, CompN), qrhub = qrhub, 
                               labsize=labsize, fontsize=fontsize, base_size = base_size, max.overlaps = max.overlaps, MaxNsig=MaxNsig)
  ggGO_neg <- enrichGOFunction(top_genes_neg, title_suffix =  " Neg", gene_universe=gene_universe, title = paste0(redLab, CompN), qrhub = qrhub, 
                               labsize=labsize, fontsize=fontsize, base_size = base_size, max.overlaps = max.overlaps, MaxNsig=MaxNsig)
  
  return(list(ggGO_pos = ggGO_pos, ggGO_neg = ggGO_neg))
}



#' PercExprThrTab
#'
#' This function calculates the percentage of cells expressing each gene above a given threshold in each metadata group of a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param GeneNames A character vector specifying the gene names to calculate the percentage for.
#' @param CutThresh A numeric value specifying the threshold for gene expression. If NULL, quantile thresholding will be used.
#' @param MetaDataName A character string specifying the metadata column name to group the cells.
#' @param slot The slot name to retrieve the data from. Default is "data".
#' @param Quant A numeric value specifying the quantile threshold to use if CutThresh is NULL.
#' @param assay The assay name to retrieve the data from. If NULL, the default assay of the Seurat object is used.
#'
#' @export
PercExprThrTab <- function(SerObj, GeneNames = NULL, CutThresh = NULL,
                           MetaDataName = NULL, slot="data", Quant=0.1, assay="RNA" ){
  # SerObj = RhAtlas452.Tcells.DS.CD8s; GeneNames= IDgenes; CutThresh = 1; MetaDataName = "CustClus"
  
  
  if(is.null(GeneNames)) stop("GeneNames is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  # if(is.null(CutThresh)) stop("CutThresh is NULL")
  
  if(is.null(CutThresh)) print("using quantile thresholding")
  
  tempAssayData = GetAssayData(object = SerObj, slot = slot, assay = assay)
  
  print(paste0("number of genes not in seurat = ", length(GeneNames[!GeneNames %in% rownames(tempAssayData) ])))
  
  GeneNames <- GeneNames[GeneNames %in% rownames(tempAssayData)]
  
  tempDGE <- tempAssayData[GeneNames, ]
  remove(tempAssayData)
  
  tempDFout <- as.data.frame( lapply(GeneNames, function(GeneName){
    
    if(is.null(CutThresh)) CutThresh = quantile(tempDGE[GeneName, ], Quant=Quant)
    
    tempTab <- table(tempDGE[GeneName, ]>CutThresh, 
                     SerObj@meta.data[,MetaDataName])
    
    if(nrow(tempTab) == 1) {
      tempTab <- as.table(rbind(tempTab, rep(0, length(tempTab))))
      rownames(tempTab) <- c("FALSE", "TRUE")
    }
    
    tempTab <- addmargins(tempTab)
    # prcs <- round(tempTab[2,]/tempTab[3,]*100, 4)[1:(ncol(tempTab)-1)]
    prcs <- round(tempTab["TRUE",]/tempTab["Sum",]*100, 4)[1:(ncol(tempTab)-1)]
    prcs[is.na(prcs)] <- 0
    prcs
  }))
  
  colnames(tempDFout) <- GeneNames
  
  
  
  return(tempDFout)
  
}

#' Plot Explained Variance of Principal Components
#'
#' This function generates a bar plot illustrating the explained variance of principal components (PCs) in a Seurat object.
#'
#' @param SerObj A Seurat object.
#' @param pcs Numeric vector specifying which PCs to plot. Default is 1:10.
#'
#' @return A ggplot object representing the explained variance of the specified principal components.
#'
#' @export
plot_pca_variance <- function(SerObj, pcs = 1:10) {
  
  # Get the explained variance by the PCs
  variance_explained <- SerObj@reductions$pca@stdev^2 / sum(SerObj@reductions$pca@stdev^2)
  
  # Subset data for the specified PCs
  variance_data <- data.frame(PC = pcs, Variance_Explained = variance_explained[pcs])
  
  # Plot the explained variance
  gg <- ggplot(variance_data, aes(x = factor(PC), y = Variance_Explained, fill = ifelse(PC > 3, "Pink", "Lightgray"))) + 
    geom_col() + 
    xlab("Principal Component") + 
    ylab("Explained Variance (%)") + 
    theme_classic(base_size = 14)+
    ggtitle(paste("Explained Variance of PCs", min(pcs), "to", max(pcs))) +
    scale_x_discrete(labels = paste("PC", pcs)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Pink", "Lightgray")) +
    geom_hline(yintercept = 0.05, linetype = "dashed") + 
    theme(legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, vjust = .8, hjust=.9))
  
  return(gg)
}


#' Create a grouped plot with percentage labels in the legend.
#'
#' This function generates a grouped plot for the specified feature in a Seurat object, 
#' adding percentage labels to the legend to represent the proportion of total points 
#' in each category.
#'
#' @param SerObj A Seurat object containing the data.
#' @param group_feature The feature based on which the data will be grouped and plotted.
#' @param color_values A vector of color values for the groups.
#' @param alpha The alpha (transparency) value for the points. Default is 0.3.
#' @param raster Logical, indicating whether raster graphics should be used. Default is FALSE.
#' @param raster.dpi A numeric vector of length 2 indicating the reSerObjlution of the raster 
#'   graphics. Default is c(2000, 2000).
#' @param pt.size The size of the points in the plot. Default is 0.5.
#' @param base_size The base font size for the plot. Default is 10.
#' @param theme_b The theme to be applied to the plot. Default is theme_bw with base_size parameter.
#' @return A ggplot object with grouped points and legend labels representing percentages.
#' @examples
#' \dontrun{
#' # Example usage:
#' # gg <- create_grouped_plot(SerObj, group_feature = "IsAlphaBeta", 
#' #                           color_values = c("dodgerblue", "maroon"), alpha = 0.3,
#' #                           raster = FALSE, raster.dpi = c(2000, 2000), pt.size = 0.5,
#' #                           base_size = 10, theme_b = theme_classic(base_size = 10))
#' }
#'
#' @export
create_grouped_plot <- function(SerObj, group_feature, color_values = c("dodgerblue", "maroon"), 
                                alpha = .3, raster=F, raster.dpi = c(2000, 2000), 
                                pt.size = .5,base_size = 10, shuffle = F,
                                theme_b = theme_classic(base_size = base_size)) {
  # Extract the specified group_feature column from the Seurat object
  group_column <- SerObj[[group_feature]]
  
  # Calculate percentage of total points for the specified group_feature
  percentages <- round(table(group_column) / nrow(group_column) * 100 , 2)
  
  gg <- DimPlot(SerObj, shuffle = shuffle, 
                group.by = group_feature, 
                raster = raster, raster.dpi = raster.dpi, 
                cols = alpha(color_values, alpha),
                pt.size = pt.size, order = T) +
    scale_color_manual(values = color_values, 
                       breaks = names(percentages), 
                       labels = paste(names(percentages), 
                                      " (", round(percentages, 1), "%)")) +
    theme_b +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    labs(color = group_feature)
  
  return(gg)
}


#' Plot 3D visualization of Seurat object
#'
#' @param SerObj A Seurat object
#' @param feature1 The name of the first feature to plot along the x-axis
#' @param feature2 The name of the second feature to plot along the y-axis
#' @param feature3 The name of the third feature to plot along the z-axis
#' @param color_feature The name of the feature to use for coloring the points (optional)
#' @param cols The color palette to use for the points (optional)
#' @param plot_title The title of the plot (optional)
#'
#' @return A 3D plotly visualization of the Seurat object
#' @export
plot3d_seurat <- function(SerObj, 
                          feature1, feature2, feature3, 
                          color_feature=NULL, 
                          cols = colorRamp(c("navy", "red")), 
                          plot_title="", slot = "data", markerSize=2) {
  library(plotly)
  
  
  
  # Check if feature1 is a gene name or metadata name
  # if (feature1 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   x <- SerObj@assays$RNA@data[feature1, ]
  # } else if (feature1 %in% colnames(SerObj@meta.data)) {
  #   x <- SerObj@meta.data[, feature1]
  # } else {
  #   stop(paste("Feature", feature1, "not found."))
  # }
  x = FetchData(SerObj, feature1, slot = slot)[,1]
  
  # Check if feature2 is a gene name or metadata name
  # if (feature2 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   y <- SerObj@assays$RNA@data[feature2, ]
  # } else if (feature2 %in% colnames(SerObj@meta.data)) {
  #   y <- SerObj@meta.data[, feature2]
  # } else {
  #   stop(paste("Feature", feature2, "not found."))
  # }
  y = FetchData(SerObj, feature2, slot = slot)[,1]
  
  # Check if feature3 is a gene name or metadata name
  # if (feature3 %in% rownames(SerObj@assays$RNA@scale.data)) {
  #   z <- SerObj@assays$RNA@data[feature3, ]
  # } else if (feature3 %in% colnames(SerObj@meta.data)) {
  #   z <- SerObj@meta.data[, feature3]
  # } else {
  #   stop(paste("Feature", feature3, "not found."))
  # }
  z = FetchData(SerObj, feature3, slot = slot)[,1]
  
  
  tempDF = data.frame(x=x, y=y, z=z)
  
  # Set point colors if color_feature is specified
  if (!is.null(color_feature)) {
    
    if (color_feature %in% colnames(SerObj@meta.data)) {
      color_vals <- SerObj@meta.data[, color_feature]
      
      if(!is.factor(color_vals)) tempDF$Cols = naturalsort::naturalfactor(color_vals) else tempDF$Cols = color_vals
      
      cols_red = cols[1:length(levels(tempDF$Cols))]
      names(cols_red) = levels(tempDF$Cols)
      
      p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, color = ~Cols, 
                   colors = cols_red, type = "scatter3d", mode = "markers", 
                   marker = list(size = markerSize))
    } else {
      GeneExpr = FetchData(SerObj, color_feature, slot = slot)[,1]
      p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, 
                   color = ~GeneExpr, colors = cols, 
                   type = "scatter3d", mode = "markers", marker = list(size = markerSize))
      
    }
  } else {
    # Set default color
    p <- plot_ly(tempDF, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", marker = list(size = 2))
  }
  
  # Set plot title and axis labels
  # plot_title <- paste(feature1, "vs", feature2, "vs", feature3)
  x_axis_label <- feature1
  y_axis_label <- feature2
  z_axis_label <- feature3
  
  p <- layout(p, title = plot_title, scene = list(xaxis = list(title = x_axis_label), 
                                                  yaxis = list(title = y_axis_label), 
                                                  zaxis = list(title = z_axis_label)))
  
  return(p)
}


#' Calculate the percent expressed cells for each group in a Seurat object
#'
#' @param SerObj A Seurat object
#' @param group.by A meta feature used to group the cells
#' @param features A vector of features (genes) to calculate the percent expressed cells for
#' @param plot A logical indicating whether to plot 
#' @param topN select top N genes  
#' @param cols color vector
#' @return A data frame with the percent expressed cells for each group and feature
#' @export
VlnPlot_subset <- function(SerObj, group.by, features, plot=F, topN = 10, xlab = "",
                           cols=c("#8DD3C7", "#33A02C", "#F4CAE4", "#A6CEE3", "#FDCDAC",
                                  "#B2DF8A","#E78AC3", "#1F78B4", "#FB9A99", "#E31A1C", 
                                  "#66C2A5", "#FF7F00"), addDimplot = T, returnLs = T){
  
  
  #fix factor
  SerObj@meta.data[,group.by] = naturalsort::naturalfactor(as.character(SerObj@meta.data[,group.by]))
  
  features = unique(features)
  
  
  PctExprDF =  scCustFx:::PercentExpressed(SerObj, 
                                           group.by=group.by, 
                                           features=features, 
                                           plot_perHM = plot)
  
  PctExprDF$max = apply(PctExprDF, 1, max)
  PctExprDF$min = apply(PctExprDF, 1, min)
  
  PctExprDF$maxminDelta = PctExprDF$max-PctExprDF$min
  
  
  Idents(SerObj) = group.by
  
  vln <- VlnPlot(SerObj, 
                 features = rownames( PctExprDF[order(PctExprDF$maxminDelta, decreasing = T)[1:topN],]), 
                 stack = TRUE, flip = TRUE, fill.by = "ident",
                 cols = cols) + NoLegend() + xlab(xlab) #+
  # theme_classic(base_size = 14) 
  
  
  if(addDimplot){
    dim <- DimPlot(SerObj, 
                   label = T, label.size = 10, label.box = T, repel = T,
                   cols=cols,
                   group.by=group.by) + 
      coord_flip()  + scale_y_reverse() +
      theme_classic(base_size = 14) + 
      NoLegend() +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank()
      )
    if(returnLs) {
      list(dim=dim, vln=vln)
    } else {
      dim / vln
    }
  } else  {
    vln
  }
  
  
}

#' @title PlotFeatThrTab
#' @description This function from a seurat object, takes the feature selected and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SerObj, A Seurat object.
#' @param FeatName, the Feature name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotFeatThrTab <- function(SerObj, FeatName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(FeatName)) stop("FeatName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(FeatName) != 1) stop ("FeatName Length != 1")
  
  tempDF <- FetchData(SerObj, vars = c(FeatName, MetaDataName))
  colnames(tempDF) <- c("var1", "var2")
  tempTab <- table(tempDF$var1>CutThresh, tempDF$var2)
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    
    
    
    vp <- VlnPlot(SerObj, features = FeatName, 
                  group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=var1)) + 
        geom_histogram(aes(y=after_stat(density)), colour="black", fill="white", binwidth = .03) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(FeatName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(data.table::melt(tempTab)) +
        geom_bar(aes(x=factor(Var2), y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
  }
  
  
  
}

#' @title PlotExprThrTab
#' @description This function from a seurat object, takes the gene expression of a single gene and provides a vizualization for thresholding and corresponding descritized tabulation.
#' @param SerObj, A Seurat object.
#' @param GeneName, the Gene name of interest. 
#' @param CutThresh, cutting threshold.
#' @param TitleExtra, title from user.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotExprThrTab <- function(SerObj, GeneName = NULL, CutThresh = NULL, PlotBar = T,
                           MetaDataName = NULL, TitleExtra = "", PlotHist = F, 
                           col_vector = col_vector){
  
  if(is.null(GeneName)) stop("GeneName is NULL")
  if(is.null(MetaDataName)) {
    warning("MetaDataName is NULL, setting MetaDataName  = 'orig.ident' ")
    MetaDataName  = "orig.ident"
  }
  if(is.null(CutThresh)) stop("CutThresh is NULL")
  if(length(GeneName) != 1) stop ("GeneName Length != 1")
  
  
  # tempTab <- table(SerObj@assays$RNA@data[GeneName, ]>CutThresh, SerObj@meta.data[,MetaDataName])
  tempTab <- table(GetAssayData(object = SerObj, slot = 'data')[GeneName, ]>CutThresh, 
                   SerObj@meta.data[,MetaDataName])
  
  
  
  print(addmargins(tempTab))
  
  if(!PlotHist){
    if(PlotBar) p2 <- ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells")
    
    print(p2)
    
  } else {
    
    tempDF  <- as.data.frame(cbind(Expr = SerObj@assays$RNA@data[GeneName, ], CellN = colnames(SerObj@assays$RNA@data)), stringsAsFactors = F)
    tempDF$Expr <- as.numeric(tempDF$Expr)
    
    vp <- VlnPlot(SerObj, features = GeneName, group.by = MetaDataName, cols = col_vector) + theme(legend.position = "none") +
      geom_hline(aes(yintercept=CutThresh),
                 color="blue", linetype="dashed", size=1)
    
    p1 <- cowplot::plot_grid(
      ggplot(tempDF, aes(x=Expr)) + 
        geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = .3) +
        geom_density(alpha=.2, fill="#E69F00") + theme_bw() + 
        geom_vline(aes(xintercept=CutThresh ),
                   color="blue", linetype="dashed", size=1) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)),
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90)) +
        ggtitle(paste0(GeneName, ">", CutThresh, " :: count", "\n", TitleExtra)) + ylab("Number of cells"),
      
      vp,
      
      ggplot(melt(tempTab)) +
        geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw()  + scale_fill_manual(values=col_vector) +
        theme(legend.position="bottom",
              legend.direction="horizontal",
              legend.title = element_blank(),
              axis.text.x = element_text(angle = 90))+
        scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                           labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
        ggtitle("\n Relative Contribution") + ylab("Relative % cells"),
      
      ncol = 2)
    
    print(p1)
    
    
    
  }
  
  
  
}

#' @title QuantileADTcounter
#'
#' @description Plots a histogram of ADT/protein/CITE-seq markers optional: identify percent of cells cut by a threshold
#' @param SerObj A Seurat Object 
#' @param min.cells a horizontal line is plotted at this value
#' @param features Features
#' @param slot which ADT slot data or count
#' @param thresh a vertical line is plotted at this value and perc. of cells below it is computed and plotted
#' @return A ggplot object
#' @export
QuantileADTcounter <- function(SerObj, features=NULL, min.cells=10, 
                               slot="data", thresh = NULL){
  print("Setting default assay to ADT")
  DefaultAssay(SerObj) = "ADT"
  
  ADT_range <- (GetAssayData(SerObj, slot=slot))
  # max(ADT_range)
  # ADT_range[1:10, 1:10]
  
  ADT_range.melt <- reshape2::melt(quantile(ADT_range, 
                                            SerObjrt(c(1:9.9/100, 
                                                   1:9/10, .75, .85, .995, .996, .997, 
                                                   .998, .999, .9993, .9995, .9998, .9999,
                                                   round(seq(from=90, to=100, length.out = 50), 2)/100)) ))
  
  ADT_range.melt$value <- round(ADT_range.melt$value, 3)
  ADT_range.melt <- unique(ADT_range.melt)
  
  # ADT_range.melt$Quant  <- as.numeric(gsub("%", "", rownames(ADT_range.melt)))/100
  
  if(is.null(features)) {
    warning("features was NULL, using all features of ADT")
    features = rownames(SerObj)
  }
  
  ggls <- lapply(features, function(adt_f, ADT.DF){
    
    # adt_f = "CD8a"
    
    ADT_range.melt.sub <- ADT_range.melt
    
    ADT_range.melt.sub$Bin  <- as.numeric(table(cut(GetAssayData(SerObj, slot = slot)[adt_f, ], 
                                                    breaks = c(-5, ADT_range.melt.sub$value), 
                                                    include.lowest = T)))
    
    
    ADT_range.melt.sub$Percent = ADT_range.melt.sub$Bin/sum(ADT_range.melt.sub$Bin)
    
    ADT_range.melt.sub$BinFac = factor(paste0("bin-", ADT_range.melt.sub$value), labels  = ADT_range.melt.sub$value)
    
    if(is.null(thresh)) thresh = as.numeric(quantile(ADT_range.melt.sub$value, .7))
    
    
    
    # sum(ADT_range.melt.sub$Percent[ADT_range.melt.sub$value < thresh])
    
    ThreshPerc  = sum(ADT_range.melt.sub[ADT_range.melt.sub$value < thresh,]$Bin)/sum(ADT_range.melt.sub$Bin)
    
    BinThresh <- which(ADT_range.melt.sub$BinFac ==
                         min(as.numeric(as.character(ADT_range.melt.sub$BinFac[as.numeric(as.character(ADT_range.melt.sub$BinFac))>thresh]))))
    
    ggplot(ADT_range.melt.sub) + 
      geom_bar(aes(x=BinFac, y=Bin, fill=Percent), position = 'dodge',stat='identity') + theme_bw()  +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      labs(title = paste0(adt_f, '\n No. of cells binned per global ADT quantiles\nPercent below Threshold (',round(thresh,2),') = %' , round(ThreshPerc*100, 3))) + 
      geom_hline(yintercept = min.cells) + geom_vline(xintercept = BinThresh, col="red", lty=2)
    
  })
  return(ggls)
}



#' @title plot_libSize
#' 
#' @description Plots a histogram of ADT/protein/CITE-seq markers (optional: identify percent of cells cut by a threshold)
#' @param SerObj A Seurat Object 
#' @param assay Assay within the seurat object
#' @param slot which ADT slot data or count
#' @param featureSplit a feature of metadata; if NULL no split id done. 
#' @return A ggplot object
#' @export
plot_libSize <- function(SerObj, assay="ADT", slot = "counts", 
                         featureSplit=NULL, col_vector = NULL){
  
  if(is.null) col_vector = scCustFx:::ColorTheme()$col_vector
  
  if(is.null(featureSplit)) stop("featureSplit is null")
    
  if(!is.null(featureSplit)) {
    
    if(!featureSplit %in% colnames(SerObj@meta.data)) {
      stop("featureSplit not in meta.data")
      
    } else {
      
      SerObj$tempLibSize = colSums(GetAssayData(SerObj, assay = assay, slot=slot))
      
      SerObj$tempGrouping = factor(SerObj@meta.data[,featureSplit])
      
      tempOut = SerObj@meta.data %>% group_by(tempGrouping)  %>% 
        summarize(mean_LibSize = sum(tempLibSize, na.rm = TRUE))
      
      
    }} 
  
  
  tempOut %>% 
    reshape2::melt() %>% 
    ggplot(aes(y=log10(value), x=tempGrouping, fill=tempGrouping)) + 
    geom_bar( position = 'dodge',stat='identity') +
    theme_classic() +
    theme(legend.position="bottom",
          legend.direction="horizontal",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90)) + 
    ggtitle(paste0(Title, "\nSum log10 library size per barcode Prefix")) +
    scale_fill_manual(values=col_vector)
  
}
                               
                               
