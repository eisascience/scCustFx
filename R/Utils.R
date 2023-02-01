#' @title is.even
#'
#' @description logical returns T if even
#' @param x, numbers
#' @export
is.even <- function(x) x %% 2 == 0

#' @title is.odd
#'
#' @description logical returns T if odd
#' @param x, numbers.
#' @return histo_numers
#' @export
is.odd <- function(x) x %% 2 != 0


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


#' @title transposedt
#'
#' @description Transpose a data.table
#' @return A transposed data.table
#' @keywords transpose, t
#' @param dt The datatable
#' @import data.table
#' @export
TransposeDT <- function(dt, varlabel="myVar") {
  dtrows = names(dt)
  dtcols = as.list(c(dt[,1]))
  dtt = transpose(dt)
  dtt[, eval(varlabel) := dtrows]
  setnames(dtt, old = names(dtt), new = c(dtcols[[1]], eval(varlabel)))
  dtt = dtt[-1,]
  setcolorder(dtt, c(eval(varlabel), names(dtt)[1:(ncol(dtt) - 1)]))
  return(dtt)
}



#' @title LogAdd
#'
#' @description Sum a vector of log valued
#' @return A vector.
#' @keywords log, sum
#' @export
LogAdd <- function(x) {
  mpi <- max(x)
  return(mpi + log(x = sum(exp(x = x - mpi))))
}


#' @title UniformSampleDF_FacPor
#'
#' @description uniformly sample a dataframe based on factor and porportion
#' @return index of the uniformly sampled samples
#' @keywords sample, uniform
#' @export
UniformSampleDF_FacPor <- function(x, ClassF, p){
  nr <- NROW(x)
  size <- (nr * p) %/% length(unique(ClassF))
  idx <- lapply(split(seq_len(nr), ClassF), function(.x) sample(.x, size))
  unlist(idx)
}


#' @title RearrangeLabelsByFrequency
#'
#' @description reassign the labels of the elements in the vector based on the frequency of their original labels. 
#' @return a vector relabeled frequency wise
#' @export
RearrangeLabelsByFrequency <- function(vec) {
  tab <- table(vec)
  ordered_tab <- sort(tab, decreasing = TRUE)
  relabeled_vec <- vec
  label <- 0
  for (val in names(ordered_tab)) {
    relabeled_vec[vec == val] <- label
    label <- label + 1
  }
  return(relabeled_vec)
}