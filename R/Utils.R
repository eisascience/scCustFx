


#' Compute the percentage of intersection between two character vectors
#'
#' This function calculates the percentage of intersection between two character vectors.
#' It converts the vectors to sets, finds the intersection, and computes the percentage
#' based on the size of the intersection relative to the first vector.
#'
#' @param vector1 A character vector.
#' @param vector2 A character vector.
#' @return The percentage of intersection between the two character vectors.
#' @examples
#' vector1 <- c("A", "B", "C", "D", "E")
#' vector2 <- c("C", "D", "E", "F", "G")
#' intersect_perc(vector1, vector2)
#' @export
intersect_perc <- function(vector1, vector2) {
  intersection <- intersect(unique(vector1), unique(vector2))
  percent_intersect <- length(intersection) / length(unique(vector1)) * 100
  return(percent_intersect)
}



#' Quantile Breaks
#'
#' @param xs A numeric vector of values.
#' @param n An integer specifying the number of quantile breaks to be calculated.
#'
#' @return A numeric vector representing the quantile breaks.
#'
#' @examples
#' xs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' quantile_breaks(xs, n = 5)
#' 
#' @export
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}


#' Check if a value falls within a distribution range defined by a set of numbers
#'
#' This function checks if a value falls within the distribution range defined by a set
#' of numbers, represented by the minimum and maximum values of the set.
#'
#' @param value The value to be checked
#' @param my_set A set of numbers defining the distribution range
#'
#' @return A logical value, TRUE if the value falls within the distribution range, FALSE otherwise
#'
#' @examples
#' my_set <- c(10, 20, 30, 40, 50)
#' is_in_range(25, my_set)
#'
#' @export
is_in_range <- function(value, my_set) {
  min_val <- min(my_set)
  max_val <- max(my_set)
  return(value >= min_val && value <= max_val)
}



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


#' Find the closest multiple of 10 to a given number.
#'
#' This function takes a numeric input \code{x} and returns the closest multiple of 10 to \code{x}.
#'
#' @param x a numeric input.
#' @return the closest multiple of 10 to \code{x}.
#' @examples
#' closest_1_0s(12345)
#' closest_1_0s(6.75)
#' @export
closest_1_0s <- function(x) {
  n <- floor(log10(x))+1 #number of digits
  closest_value <- round(x/10^(n-1))*10^(n-1)
  return(closest_value)
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





#' @title frac_to_numeric
#' @description Converts fraction to numeric
#' @return A vector
#' @keywords numeric
#' @export
frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))






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


#' @title WichIn1not2
#'
#' @description compares unique association between two labels from A dataftrame usually from GO annotation with a column named cluster
#' @return unique genes to the comparison
#' @keywords cluster, gene, 
WichIn1not2 <- function(Clus1N = c(1), DataT = "", Clus2N = c(2)){
  Gs1 <- subset(DataT, cluster %in% Clus1N)$gene 
  Gs2 <- subset(DataT, cluster %in% Clus2N)$gene
  Gs1[which(!Gs1 %in% Gs2)]
  
}

