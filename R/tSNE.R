

#' This function to run tsne on matrix
#'
#' @param myMat SDA score matrix where cols are the comps and rows are cells
#' @param save_path path to save the tsne rds obj
#' @param tsnepp tSNE preplexity
#' @param tSNE_n.iter N iter for tSNE
#' @param num_threads number of threads
#' @param check_duplicates check for dup rows
#' @param selectComps select SDA compos to run the tSNE
#' @return A list of two vector text of enrichments, one for pos one for neg scored cells
#' @export
.run_tSNE_mat <- function(myMat, save_path, tsnepp, tSNE_n.iter, 
                          num_threads = 8, check_duplicates = F, selectComps = NULL, overwrite = F){
  
  if(is.null(selectComps)) selectComps = colnames(myMat)
  
  if(!file.exists(save_path) || overwrite == TRUE){    
    
    tsne_res <- Rtsne::Rtsne(myMat[,selectComps], 
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
  return(tsne_res)
}

