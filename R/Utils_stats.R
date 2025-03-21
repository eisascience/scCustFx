
#' @title find_peaks
#' @description Returns the maxima of points. To get minima, -1*x
#' @param x, A vector of numbers, if small or not very smooth, use a smoothing function. Like density(x).
#' @param m, An integer that acts as a loose hand for resolution.
#' @return vector of peaks positions.
#' @export
find_peaks <- function (x, m = 4){
  #https://github.com/stas-g/findPeaks
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
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