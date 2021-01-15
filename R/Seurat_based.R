
#' @title ReduceSerObj_HTO
#'
#' @description Remove unwanted HTOS from Seurat objects
#' @param so A Seurat Object 
#' @param removeHTOs Vector of names of HTOs to remove, if NULL, removeHTOs = c("Discordant", "Doublet", "ND")
#' @return A subset Seurat object.
#' @export
ReduceSerObj_HTO <- function(so, removeHTOs = NULL){
  if(is.null(removeHTOs)){
    print("removeHTOs is null, therefore, deep cleaning")
    removeHTOs = c("Discordant", "Doublet", "ND")
  }
  if (sum(removeHTOs %in% names(table(so$HTO)))==0){
    stop("removeHTOs defined not in seurat object ")
  } else {
    keepHTOs <- setdiff(names(table(so$HTO)), removeHTOs)
  }
  
  print("Keept HTOs: ")
  print(keepHTOs)
  return(subset(so, HTO %in% keepHTOs))
}



