#' GetKEGGPathways
#'
#' @param df data frame containing differential gene expression results
#' @param species species to map . Available option: "human" and "mouse".
#'
#' @return Data frame containing KEGG pathways
#' @export
GetKEGGPathways <- function(df = NULL, species = "human"){
  
  up <- df %>% filter(avg_log2FC > 0) %>% pull(gene)
  down <- df %>% filter(avg_log2FC < 0) %>% pull(gene)
  
  if(species == "human"){
    if (length(up) > 0)
      up <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, up, "ENTREZID", "SYMBOL")
    if (length(down) > 0)
      down <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, down, "ENTREZID", "SYMBOL")
    
    genes <- list(Up = up, Down = down)
    pathways <- limma::kegga(genes, species = "Hs")
    
  }else if(species == "mouse"){
    if (length(up) > 0)
      up <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, up, "ENTREZID", "SYMBOL")
    if (length(down) > 0)
      down <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, down, "ENTREZID", "SYMBOL")
    
    genes <- list(Up = up, Down = down)
    pathways <- limma::kegga(genes, species = "Mm")
  }
  
  # adjust p values
  pathways <- pathways %>% mutate(adj.P.Up = p.adjust(pathways[["P.Up"]], method = "bonferroni"))
  pathways <- pathways %>% mutate(adj.P.Down = p.adjust(pathways[["P.Down"]], method = "bonferroni"))
  
  # signed log p-values
  pathways <- pathways %>% mutate(logP.Up = -log10(adj.P.Up)) # -log to obtain positive values for UP
  pathways <- pathways %>% mutate(logP.Down = log10(adj.P.Down))# log to obtain negative values for UP
  
  return(pathways)
}

#' GetTopKEGG
#'
#' @param df Data frame containing KEGG pathways.
#' @param top.by Label used to select the top.
#' @param type Type of selected p-values. Available options are: "double", negatives and positive, "positive", and "negative".
#' @param ntop Number of top values.
#'
#' @return Data frame containing the top values.
#' @export
GetTopKEGG <- function(df = NULL, top.by = "Selected_PV", type = "double" , ntop = 20){
  
  if(type == "double"){
    idx <- sort(abs(df[[top.by]]), decreasing = TRUE, index.return = TRUE)$ix
  }else if(type == "positive"){
    idx <- sort(df[[top.by]], decreasing = TRUE, index.return = TRUE)$ix
  }else{
    idx <- sort(df[[top.by]], decreasing = FALSE, index.return = TRUE)$ix
  }
  
  idx <- idx[1:ntop]
  tmp <- df[idx,]
  
  return(tmp)
}

#' GetTopPathways
#'
#' @param dfs List of dataframes containing KEGG pathways
#' @param ntop Number of top pathways to select.
#'
#' @return List containing the selected pathways and their ranks (e.g. number of dfs containing the pathway).
#' @export
GetTopPathways <- function(dfs = NULL, ntop = 20){
  dfs <- lapply(dfs, GetTopKEGG, ntop = ntop)
  tmp <- Reduce(rbind,dfs)
  
  idx <- which(tmp$Selected_PV == 0)
  tmp <- tmp[-idx,]
  
  tmp <- tmp %>% group_by(Pathway) %>% count(Pathway)
  tmp <- SortDF(df = tmp, sort.by = "n")
  sel_pathways <- tmp[1:ntop,]
  
  return(sel_pathways)
}


#' SetDirectionKEGG
#'
#' @param df Data frame containing KEGG analysis results.
#'
#' @return Input data frame with two extra columns: Direction (e.g. Up and Down for up and down expressed pathways, respectively), and Selected_PV (e.g. -log10(pval) of the selected up or down p-value)
#' @export
SetDirectionKEGG <- function(df = NULL){
  # SETUP
  ## setup direction and selected p-values
  
  # first column are positive pv, second column are negative pv
  tmp <- cbind(abs(df$logP.Up), abs(df$logP.Down))
  idxMax <- apply(tmp, MARGIN = 1, which.max)
  direction <- rep(NA, nrow(tmp))
  direction[(idxMax == 1)] <- "up"
  direction[(idxMax == 2)] <- "down"
  
  selPV <- rep(NA, nrow(tmp))
  selPV[(idxMax == 1)] <- tmp[(idxMax == 1),1]
  selPV[(idxMax == 2)] <- -tmp[(idxMax == 2),2]
  
  df <- df %>% mutate(Direction = direction)
  df <- df %>% mutate(Selected_PV = selPV)
  
  return(df)
}


#' PlotKEGG
#'
#' @param df Data frame containing KEGG pathways
#' @param type Type of selected p-values. Available options are: "double", negatives and positive, "positive", and "negative".
#' @param color_by Label to color points.
#' @param point_size Size of points.
#' @param text_size Size of text in the plot.
#' @param sort.by Label to use to sort the pathways.
#' @param decreasing Sorting direction.
#' @param threshold P-value threshold.
#' @param vline_size Threshold line size.
#'
#' @return A ggplot object.
#' @export
PlotKEGG <- function(df = NULL, type = "double", color_by = "Direction", point_size = 3, text_size = 15, sort.by = "Selected_PV", decreasing = FALSE, threshold = NULL, vline_size = 1 ){
  
  df <- SortDF(df = df, sort.by = sort.by, decreasing = decreasing)
  levels <- unique(df$Pathway)
  df[["Pathway"]] <- factor(df[["Pathway"]], levels = levels)
  
  p <- ggplot(df, aes_string(x = "Selected_PV", y = "Pathway", color = color_by)) +
    geom_point(size = point_size) +
    theme_classic() +
    labs(x = expression(paste('Signed ', -log[10],'(adjusted p-value)')),
         y = "KEGG Pathways") +
    TextSize(size = text_size)
  
  if(!is.null(threshold)){
    idx <- which(abs(df$Selected_PV) < -log10(threshold))
    df_t <- df[idx,]
    
    p <- p +
      geom_point(data = df_t, aes_string(x = "Selected_PV", y = "Pathway"), color = "gray",size = point_size) + geom_vline(xintercept = c(-log10(threshold),log10(threshold)), linetype = "dashed", size = vline_size)
  }
  
  return(p)
}

