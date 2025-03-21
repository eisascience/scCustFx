

#' This function to run tsne on matrix
#'
#' @param myMat SDA score matrix where cols are the comps and rows are cells
#' @param save_path default NULL ignores saving, or set path to save the tsne rds obj
#' @param tsnepp tSNE preplexity
#' @param tSNE_n.iter N iter for tSNE
#' @param num_threads number of threads
#' @param check_duplicates check for dup rows
#' @param selectFeats select Features from Matrix to run the tSNE
#' @return A list of two vector text of enrichments, one for pos one for neg scored cells
#' @export
RuntSNE.Matrix <- function(myMat, save_path = NULL, tsnepp, tSNE_n.iter, 
                          num_threads = 8, check_duplicates = F, 
                          selectFeats = NULL, overwrite = F){
  
  if(is.null(selectFeats)) selectFeats = colnames(myMat)
  
  if(is.null(save_path)){
    tsne_res <- Rtsne::Rtsne(myMat[,selectFeats], 
                             verbose=T, pca=F, normalize = F,
                             perplexity = tsnepp,  
                             max_iter=tSNE_n.iter, 
                             num_threads = num_threads, 
                             check_duplicates = check_duplicates,
                             partial_pca = FALSE,
                             # is_distance = FALSE,
                             # Y_init = NULL,
                             pca_center = FALSE,
                             pca_scale = FALSE,
                             # stop_lying_iter = ifelse(is.null(Y_init), 250L, 0L),
                             # mom_switch_iter = ifelse(is.null(Y_init), 250L, 0L),
                             momentum = 0.5,
                             final_momentum = 0.8,
                             eta = 200,
                             exaggeration_factor = 12)
  } else {
    
    if(!file.exists(save_path) || overwrite == TRUE){    
      
      tsne_res <- Rtsne::Rtsne(myMat[,selectFeats], 
                               verbose=T, pca=F, normalize = F,
                               perplexity = tsnepp,  
                               max_iter=tSNE_n.iter, 
                               num_threads = num_threads, 
                               check_duplicates = check_duplicates,
                               partial_pca = FALSE,
                               # is_distance = FALSE,
                               # Y_init = NULL,
                               pca_center = FALSE,
                               pca_scale = FALSE,
                               # stop_lying_iter = ifelse(is.null(Y_init), 250L, 0L),
                               # mom_switch_iter = ifelse(is.null(Y_init), 250L, 0L),
                               momentum = 0.5,
                               final_momentum = 0.8,
                               eta = 200,
                               exaggeration_factor = 12)
      
      saveRDS(tsne_res, file=save_path)
      
      
      
    } else {
      tsne_res <- readRDS(save_path)
    }
  }
  
  
  return(tsne_res)
}



#' @title gg_tSNE_df
#'
#' @description Takes a data.frame with tSNE1 and tSNE2 precomputed and plots it.
#' @param DF A data.frame 
#' @param SavePath A data.frame 
#' @param Markers A char. vec. that that all must exist in the data.frame.
#' @param splitby a single character which data.frame must contain this value in the cols; the resulting plot is faceted by this.
#' @param verbose if T verbose output is printed 
#' @param Title A title else ''
#' @param gg_theme a ggplot theme or NULL
#' @return A ggplot 
#' @export
gg_tSNE_df <- function(DF, 
                       SavePath=NULL, 
                       Markers=NULL, 
                       splitby=NULL,
                       verbose = T,
                       size=0.2, alpha=0.6,
                       Title="",
                       width = 2000, height = 2300, pointsize = 20, ncol=4, 
                       limits = c(0, 1),
                       gg_theme=NULL, printPlot = T, vridis_col_opt = "A"){
  
  if(! ("tSNE1" %in% colnames(DF)) & "tSNE2" %in% colnames(DF)) stop("tSNE1 and/or tSNE2 not found in df")
  if(is.null(Markers)) stop("Markers is NULL")
  
  ExprPlots <- list()
  
  if(is.null(splitby)){
    
    for(xN in Markers){
      if(verbose) print(xN)
      
      
      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", xN)], id.vars=c("tSNE1", "tSNE2")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        ggthemes::scale_colour_gradient2_tableau('pi0', low = "blue", mid = "white", high = "red", midpoint = 0)
      scale_color_viridis(option=vridis_col_opt,limits = limits, oob = scales::squish) +
        ggtitle(paste0(Title, "\ntSNE - ", xN))
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  } else {
    
    if(! splitby %in% colnames(DF)) stop("splitby not in colnames of DF")
    
    DF$VarSplit = DF[,splitby]
    
    for(xN in Markers){
      # xN = Markers[1]
      if(verbose) print(xN)
      
      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", "VarSplit", xN)], id.vars=c("tSNE1", "tSNE2", "VarSplit")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        scale_color_viridis(option="A") +
        ggtitle(paste0(Title, "\ntSNE - ", xN, "\n per", splitby))+ facet_wrap(~VarSplit)
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  }
  
  
  if(!is.null(SavePath)) {
    
    png(SavePath, width = width, height = height, pointsize = pointsize)
    
    print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol))
    
    dev.off()
    
    
  } else { 
    if(printPlot) {
      return(print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol)))
    } else {
      return(ExprPlots)
    } 
    
  }
  
}


#' @title rotate_tsne
#' @description rotates a tsne object 
#' @param tsne, a tsne object with tsne$Y
#' @param angle, rotation angle
#' @return rotated tsne Y
#' @export
rotate_tsne <- function(tsne, angle){
  angle <- (-angle * pi) / (-180)
  rotm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2)
  tsne$Y %*% rotm
}