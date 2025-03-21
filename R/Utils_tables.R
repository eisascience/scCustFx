#' Covert a table to a dataframe
#'

#'
#' @param tbl table
#' @param do.prop if T, does prop and need to set prop.marg
#' @param prop.marg either 1 or 2 passed to prop.table margin 
#' @return a dataframe version of the input table
#' 
#' @export
table2DF = function(tbl = NULL, do.prop = F, prop.marg = 1){
  
  
  if(is.null(tbl)) stop()
  
  if(do.prop){
    tmpDF = as.data.frame(as.matrix.data.frame(prop.table(tbl, margin = prop.marg)))
    colnames(tmpDF) = colnames(tbl)
    rownames(tmpDF) = rownames(tbl)
    
  } else {
    tmpDF = as.data.frame(as.matrix.data.frame(tbl))
    colnames(tmpDF) = colnames(tbl)
    rownames(tmpDF) = rownames(tbl)
  }
  
  tmpDF
}


#' @title discretize_table
#'
#' @description discretize a table; if tbl>0 then 1 else 0
#' @param tbl, numbers.
#' @return histo_numers
#' @export
discretize_table <- function(tbl=NULL){
  tbl[tbl>0] = 1
  tbl
} 