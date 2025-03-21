

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


#' Generate Ridge Plot of Component Scores by Metadata Feature
#'
#' This function creates a ridge plot to visualize the distribution of component scores
#' split by a specified metadata feature.
#'
#' @param scores_matrix A numeric matrix where rows are samples (with rownames) and columns are components.
#' @param meta_df A data frame with rownames matching those in `scores_matrix` and containing metadata features.
#' @param component Integer or character. The component column (e.g., `"Component_1"`) to visualize.
#' @param meta_feature Character. The metadata column to split the ridge plot by (e.g., `"seqrun"`).
#' @param title Character. Plot title (default: `"Ridge Plot of Component Scores"`).
#'
#' @return A ggplot object displaying the ridge plot.
#'
#'
#'
#' @examples
#' \dontrun{
#' if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggridges", quietly = TRUE)) {
#'     # Example scores matrix
#'     scores <- matrix(rnorm(500), nrow = 100, ncol = 5)
#'     rownames(scores) <- paste0("Sample", 1:100)
#'     colnames(scores) <- paste0("Component_", 1:5)
#'
#'     # Example metadata dataframe
#'     meta <- data.frame(
#'         seqrun = sample(c("Run1", "Run2", "Run3"), 100, replace = TRUE),
#'         InfectionGroup = sample(c("Control", "Infected"), 100, replace = TRUE),
#'         MonkeyName = sample(LETTERS[1:5], 100, replace = TRUE),
#'         TreatmentCode = sample(c("A", "B", "C"), 100, replace = TRUE),
#'         PD1status = sample(c("High", "Low"), 100, replace = TRUE),
#'         row.names = rownames(scores)
#'     )
#'
#'     # Generate ridge plot for Component 1 split by seqrun
#'     plot_component_ridge(scores, meta, component = "Component_1", meta_feature = "seqrun")
#' }
#' }
#' @export
plot_component_ridge <- function(scores_matrix, meta_df, component, meta_feature, 
                                 title = "Ridge Plot of Component Scores") {
  # Ensure required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggridges", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'ggridges' are required. Install them using install.packages('ggplot2', 'ggridges').")
  }
  
  # Convert input to data frame and ensure rownames match
  scores_df <- as.data.frame(scores_matrix)
  scores_df$SampleID <- rownames(scores_df)
  
  # Check if component exists in the matrix
  if (!(component %in% colnames(scores_df))) {
    stop(paste("Component", component, "not found in scores matrix."))
  }
  
  # Merge scores with metadata
  merged_df <- merge(scores_df, meta_df, by = "row.names", all.x = TRUE)
  merged_df <- merged_df[, c("Row.names", component, meta_feature)]  # Keep relevant columns
  colnames(merged_df) <- c("SampleID", "ComponentScore", "MetaFeature")  # Standardized names
  
  # Generate ridge plot
  p <- ggplot(merged_df, aes(x = ComponentScore, y = as.factor(MetaFeature), fill = MetaFeature)) +
    ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2) +
    theme_minimal() +
    labs(title = title, x = paste("Scores for", component), y = meta_feature, fill = meta_feature) +
    theme(legend.position = "right")
  
  return(p)
}

#' Analyze Differences in Component Scores Across Groups
#'
#' This function performs statistical comparisons between groups for a given component in the scores matrix,
#' using parametric (ANOVA/t-test) or non-parametric (Kruskal-Wallis/Wilcoxon) tests.
#'
#' @param scores_matrix A numeric matrix where rows are samples and columns are components.
#' @param meta_df A data frame with rownames matching `scores_matrix`, containing metadata.
#' @param component Character. The component column (e.g., `"Component_1"`) to analyze.
#' @param meta_feature Character. The metadata column to group by (e.g., `"seqrun"`).
#' @param verbose Logical. Whether to print test results (default: TRUE).
#'
#' @return A list containing:
#' \item{summary_stats}{A summary table with mean, median, variance per group.}
#' \item{test_result}{Results of ANOVA/Kruskal-Wallis or t-test/Wilcoxon.}
#'
#' @importFrom stats aov kruskal.test t.test wilcox.test shapiro.test aggregate
#' @importFrom dplyr group_by summarise mutate
#' @export
#'
#' @examples
#' \dontrun{
#' analyze_component_scores(scores_matrix = results$scores, 
#'                          meta_df = MetaDF, 
#'                          component = "Component_1", 
#'                          meta_feature = "seqrun")
#' }
analyze_component_scores <- function(scores_matrix,
                                     meta_df,
                                     component,
                                     meta_feature, verbose = TRUE) {
  
  # Convert input to data frame and ensure rownames match
  scores_df <- as.data.frame(scores_matrix)
  scores_df$SampleID <- rownames(scores_df)
  
  # Check if component exists in the matrix
  if (!(component %in% colnames(scores_df))) {
    stop(paste("Component", component, "not found in scores matrix."))
  }
  
  # Merge scores with metadata
  merged_df <- merge(scores_df, meta_df, by = "row.names", all.x = TRUE)
  merged_df <- merged_df[, c("Row.names", component, meta_feature)]  # Keep relevant columns
  colnames(merged_df) <- c("SampleID", "ComponentScore", "MetaFeature")  # Standardized names
  
  # Summarize statistics by group
  summary_stats <- merged_df %>%
    dplyr::group_by(MetaFeature) %>%
    dplyr::summarise(
      mean = mean(ComponentScore, na.rm = TRUE),
      median = median(ComponentScore, na.rm = TRUE),
      variance = var(ComponentScore, na.rm = TRUE),
      n = n()
    )
  
  # Check normality
  shapiro_p <- shapiro.test(merged_df$ComponentScore)$p.value
  print("Shapiro P is: ")
  print(shapiro_p)
  is_normal <- shapiro_p > 0.05  # p > 0.05 suggests normal distribution
  
  # Choose appropriate statistical test
  unique_groups <- length(unique(merged_df$MetaFeature))
  
  if (unique_groups > 2) {
    if (is_normal) {
      test_result <- aov(ComponentScore ~ MetaFeature, data = merged_df)
      test_type <- "ANOVA"
    } else {
      test_result <- kruskal.test(ComponentScore ~ MetaFeature, data = merged_df)
      test_type <- "Kruskal-Wallis"
    }
  } else {
    if (is_normal) {
      test_result <- t.test(ComponentScore ~ MetaFeature, data = merged_df)
      test_type <- "t-test"
    } else {
      test_result <- wilcox.test(ComponentScore ~ MetaFeature, data = merged_df)
      test_type <- "Wilcoxon rank-sum"
    }
  }
  
  if (verbose) {
    message("Summary Statistics:\n")
    print(summary_stats)
    message("\nUsing ", test_type, " Test:")
    print(test_result)
  }
  
  return(list(summary_stats = summary_stats, test_result = test_result))
}


#' Screen All Components for Differences by Metadata Feature
#'
#' This function systematically analyzes all components in a score matrix, testing their differences across 
#' a specified metadata feature. Components are ranked by statistical significance.
#'
#' @param scores_matrix A numeric matrix where rows are samples and columns are components.
#' @param meta_df A data frame with rownames matching `scores_matrix`, containing metadata.
#' @param meta_feature Character. The metadata column to group by (e.g., `"seqrun"`).
#' @param top_n Integer. Number of top components to visualize (default: 5).
#' @param method Character. Statistical test to use (`"anova"` for normal data, `"kruskal"` for non-parametric, default: `"auto"`).
#'
#' @return A list with:
#' \item{ranked_results}{A dataframe ranking components by significance (p-values & effect sizes).}
#' \item{plots}{A list of ggplot objects for top components.}
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' screen_results <- screen_components_by_metadata(scores_matrix = results$scores, 
#'                                                 meta_df = MetaDF, 
#'                                                 meta_feature = "seqrun")
#' print(screen_results$ranked_results) # Ranked table of components
#' screen_results$plots[[1]]  # View first plot
#' }
screen_components_by_metadata <- function(scores_matrix, meta_df, meta_feature, 
                                          top_n = 5, method = "auto") {
  # Ensure required packages are installed
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggridges", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'ggridges' are required. Install them using install.packages('ggplot2', 'ggridges').")
  }
  
  # Convert input to data frame and ensure rownames match
  scores_df <- as.data.frame(scores_matrix)
  scores_df$SampleID <- rownames(scores_df)
  
  # Merge scores with metadata
  merged_df <- merge(scores_df, meta_df, by = "row.names", all.x = TRUE)
  colnames(merged_df)[1] <- "SampleID"  # Ensure the sample ID column is named
  
  # Store results
  results_list <- list()
  
  for (component in colnames(scores_matrix)) {
    # Extract relevant columns
    test_df <- merged_df[, c("SampleID", component, meta_feature)]
    colnames(test_df) <- c("SampleID", "ComponentScore", "MetaFeature")
    
    
    
    # Choose test method
    if (method == "auto") {
      # Check normality
      shapiro_p <- shapiro.test(test_df$ComponentScore)$p.value
      print("Shapiro P is: ")
      print(shapiro_p)
      
      is_normal <- shapiro_p > 0.05
      
      test_type <- if (is_normal) "ANOVA" else "Kruskal-Wallis"
    } else if (method == "anova") {
      test_type <- "ANOVA"
    } else {
      test_type <- "Kruskal-Wallis"
    }
    
    # Perform statistical test
    if (test_type == "ANOVA") {
      test_result <- aov(ComponentScore ~ MetaFeature, data = test_df)
      p_value <- summary(test_result)[[1]][["Pr(>F)"]][1]
    } else {
      test_result <- kruskal.test(ComponentScore ~ MetaFeature, data = test_df)
      p_value <- test_result$p.value
    }
    
    # Calculate effect size (Cohen's d for top 2 groups)
    top_groups <- names(sort(table(test_df$MetaFeature), decreasing = TRUE))[1:2]
    eff_size <- if (length(top_groups) >= 2) {
      effsize::cohen.d(
        test_df$ComponentScore[test_df$MetaFeature == top_groups[1]],
        test_df$ComponentScore[test_df$MetaFeature == top_groups[2]]
      )$estimate
    } else {
      NA
    }
    
    # Store results
    results_list[[component]] <- data.frame(
      Component = component,
      Test = test_type,
      P_Value = p_value,
      Effect_Size = eff_size
    )
  }
  
  # Combine results and rank by significance
  ranked_results <- do.call(rbind, results_list) %>%
    dplyr::arrange(P_Value) %>%
    dplyr::mutate(Significance = ifelse(P_Value < 0.05, "*", ""))
  
  # Generate plots for the top components
  top_components <- head(ranked_results$Component, top_n)
  # plots <- lapply(top_components, function(comp) {
  #   ggplot(merged_df, aes(x = get(comp), y = as.factor(get(meta_feature)), fill = get(meta_feature))) +
  #     ggridges::geom_density_ridges(alpha = 0.7, scale = 1.2) +
  #     theme_minimal() +
  #     labs(title = paste("Ridge Plot of", comp, "by", meta_feature),
  #          x = comp, y = meta_feature, fill = meta_feature) +
  #     theme(legend.position = "right")
  # })
  
  return(list(ranked_results = ranked_results#, plots = plots
  ))
}


#' Summarize Significant Components Across Metadata Features
#'
#' This function extracts significant components (p < 0.05) from multiple metadata feature tests
#' and creates a summary table showing which components are shared across features.
#'
#' @param screenLS A named list where each element contains the output of `screen_components_by_metadata()`
#'                 for a metadata feature.
#'
#' @return A list containing:
#' \item{summary_table}{A binary presence/absence matrix showing which components are significant for each feature.}
#' \item{upset_plot}{An UpSet plot visualizing shared and unique significant components.}
#'
#' @importFrom dplyr bind_rows distinct mutate across
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_minimal labs
#' @importFrom ComplexUpset upset
#'
#' @export
#'
#' @examples
#' \dontrun{
#' significance_summary <- summarize_significant_components(screenLS)
#' print(significance_summary$summary_table)
#' significance_summary$upset_plot
#' }
summarize_significant_components <- function(screenLS) {
  # Ensure required packages are installed
  if (!requireNamespace("tidyr", quietly = TRUE) || 
      !requireNamespace("dplyr", quietly = TRUE) || 
      !requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Packages 'tidyr', 'dplyr', 'ggplot2', and 'ComplexUpset' are required. Install them using install.packages().")
  }
  
  # Extract significant components (p < 0.05) from each metadata feature
  sig_components_list <- lapply(names(screenLS), function(feature) {
    sig_components <- screenLS[[feature]]$ranked_results %>%
      dplyr::filter(P_Value < 0.05) %>%
      dplyr::select(Component) %>%
      dplyr::mutate(MetaFeature = feature)
    return(sig_components)
  })
  
  # Combine into a single dataframe
  sig_components_df <- dplyr::bind_rows(sig_components_list)
  
  # Create a presence/absence matrix for components across metadata features
  summary_table <- sig_components_df %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = MetaFeature, values_from = MetaFeature,
                       values_fn = length, values_fill = 0) %>%
    dplyr::mutate(across(-Component, ~ ifelse(. > 0, 1, 0))) # Convert to binary matrix
  
  # Create UpSet plot
  upset_plot <- ComplexUpset::upset(summary_table, 
                                    names(summary_table)[-1], 
                                    name = "Metadata Features",
                                    width_ratio = 0.2) +
    ggplot2::labs(title = "Shared and Unique Significant Components Across Metadata Features")
  
  return(list(summary_table = summary_table, upset_plot = upset_plot))
}