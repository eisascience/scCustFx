
#' Calculate gene overlap between two data frames
#'
#' This function takes in two data frames and calculates the gene overlap within each data frame as well as the overlap between the two data frames.
#'
#' @param df1 A data frame containing genes of interest.
#' @param df2 A data frame containing genes of interest.
#'
#' @return A list containing three matrices:
#' \describe{
#' \item{\code{df1_ol}}{A matrix of gene overlap within \code{df1}.}
#' \item{\code{df2_ol}}{A matrix of gene overlap within \code{df2}.}
#' \item{\code{df12_ol}}{A matrix of gene overlap between \code{df1} and \code{df2}.}
#' }
#'
#'
#' @export
string_overlap <- function(df1, df2) {
  
  # Initialize empty matrices for gene overlap
  overlaps1 <- matrix(0, ncol(df1), ncol(df1))
  overlaps2 <- matrix(0, ncol(df2), ncol(df2))
  overlaps12 <- matrix(0, ncol(df1), ncol(df2))
  
  # Calculate gene overlap within df1
  for (i in 1:(ncol(df1) - 1)) {
    for (j in (i+1):ncol(df1)) {
      overlaps1[i,j] <- sum(df1[,i] == df1[,j])
    }
  }
  
  # Calculate gene overlap within df2
  for (i in 1:(ncol(df2) - 1)) {
    for (j in (i+1):ncol(df2)) {
      overlaps2[i,j] <- sum(df2[,i] == df2[,j])
    }
  }
  
  # Calculate gene overlap between df1 and df2
  for (i in 1:(ncol(df1) - 1)) {
    for (j in 1:ncol(df2)) {
      overlaps12[i,j] <- sum(df1[,i] == df2[,j])
    }
  }
  
  SumMat = overlaps1 + overlaps2 + overlaps12
 
  list(df1_ol=overlaps1, df2_ol=overlaps2, , df12_ol=overlaps12)
}



#' @title similarity overlap dataframe
#'
#' @description returns a jaccard similarity matrix comparing similarity overlap of character items
#' @param df,input dataframe where each column contains N genes/strings
#' @export
similarity_jaccard_DF <- function(df, plot=T) {

  # ncol(df)
  
  
  dist <- apply(expand.grid(1:ncol(df), 1:ncol(df)), 1, function(x) {
    # print(x)
    bayesbio::jaccardSets(df[,x[1]], df[,x[2]])
  })
  
  
  dist_df <- cbind(( expand.grid(1:ncol(df), 1:ncol(df)) ), dist)
  
  head(dist_df)
  dist_df$Var1 = naturalsort::naturalfactor(paste0("v", dist_df$Var1))
  dist_df$Var2 = naturalsort::naturalfactor(paste0("v", dist_df$Var2))
  
  jaccard_mat = reshape2::dcast(dist_df, Var1 ~ Var2, value.var = "dist")[,-1]
  
  
  # jaccard_mat[1:5, 1:5]
  # dim(jaccard_mat)
  
  colnames(jaccard_mat) = colnames(df)
  rownames(jaccard_mat) = colnames(df)
  
 if(plot) pheatmap::pheatmap(jaccard_mat[,])
  
  return(jaccard_mat)
}




#' @title similarity_overlap list
#'
#' @description returns a jaccard similarity matrix comparing similarity overlap of character items
#' @param VecLs, a list of character vectors
#' @export
similarity_jaccard_LS = function(VecLs, plot=F){
  similarity_matrix <- matrix(0, nrow = length(VecLs), ncol = length(VecLs))
  for (i in 1:(length(VecLs)-1)) {
    for (j in (i+1):length(VecLs)) {
      set1 <- VecLs[[i]]
      set2 <- VecLs[[j]]
      intersection <- length(intersect(set1, set2))
      union <- length(union(set1, set2))
      jaccard_similarity <- intersection / union
      similarity_matrix[i, j] <- jaccard_similarity
      similarity_matrix[j, i] <- jaccard_similarity
    }
  }
  if(!is.null(names(VecLs))){
    colnames(similarity_matrix) = names(VecLs)
    rownames(similarity_matrix) = names(VecLs)
  }
  
  if(plot) pheatmap::pheatmap(jaccard_mat[,])
  
  return(similarity_matrix)
}

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

#' @title odds
#'
#' @description returns odd values only
#' @param x, numbers.
#' @return odd values numerical vector
#' @export
odds <- function(x){
  x[is.odd(x)]
}

#' @title evens
#'
#' @description returns even values only
#' @param x, numbers.
#' @return even values numerical vector
#' @export
evens <- function(x){
  x[is.even(x)]
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


#' Counts the overlap of items in rows per column feature
#'
#' This function takes in a list of ggplots and makes a GIF
#'
#' @param df a dataframe
#' @return rturns a count matrix
#' @export
count_overlap <- function(df) {
  overlap <- unique(unlist(df))
  
  count_matrix <- matrix(0, nrow = length(overlap), ncol = ncol(df), dimnames = list(overlap, colnames(df)))
  
  for (i in 1:ncol(df)) {
    overlapping_genes <- intersect(overlap, df[, i])
    count_matrix[overlapping_genes, i] <- rep(1, length(overlapping_genes))
  }
  return(count_matrix)
}