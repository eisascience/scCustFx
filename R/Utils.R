
#' Create an UpSet plot from a list of sets
#'
#' This function generates an UpSet plot from a list of sets and saves it as an image.
#'
#' @param set_list A list of sets to visualize.
#' @param output_filename The filename for the output image.
#' @param image_format The format in which to save the image (default: "svg").
#' @param do.sqr Should the y-axis of the main plot and overlap plot be square-root transformed? (default: TRUE)
#'
#' @return A ggplot2 UpSet plot.
#'
#' @import ggplot2
#' @import scales
#' @import gtable
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import grid
#' @import cowplot
#' @import RColorBrewer
#' @import ComplexHeatmap
#' @examples
#' # Example usage:
#' create_upset_plot(list("A" = c("1", "2", "3"), "B" = c("2", "3", "4")), "upset_plot.svg")
#'
#'
#' @export
create_upset_plot <- function(set_list, 
                              output_filename = NULL, 
                              image_format = "png", do.sqr =T,
                              plot.nrow = 2,
                              plot.ncol = 2,
                              plot.heights = c(1, 2),
                              plot.widths = c(1, 0.8),
                              save.height = 3.5, save.width = 3) {
  library(ComplexHeatmap)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(tidyr)
  
  # Create a matrix from the list of sets
  comb_mat <- make_comb_mat(set_list)
  
  # Get set names
  my_names <- set_name(comb_mat)
  
  # Total set size
  my_set_sizes <- set_size(comb_mat) %>%
    as.data.frame() %>%
    rename(sizes = ".") %>%
    mutate(Set = row.names(.))
  
  total_size <- sum(my_set_sizes$sizes)
  my_set_sizes$percentage <- paste0(my_set_sizes$sizes, " (", scales::percent(my_set_sizes$sizes / total_size), ")")
  
  
  p1 <-  my_set_sizes %>%
    mutate(Set = reorder(Set, sizes)) %>%
    ggplot(aes(x = Set, y = sizes)) +
    geom_bar(stat = "identity", aes(fill = Set), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = percentage), size = 5, angle = 90, hjust = 0, y = 1) +
    scale_fill_manual(values = brewer.pal(4, "Set2"), limits = my_names) +
    labs(x = NULL, y = "Set size", fill = NULL) +
    theme_classic() +
    theme(legend.position = "right", text = element_text(size = 14),
          axis.ticks.y = element_blank(), axis.text = element_blank())
  
  if(do.sqr) p1 = p1 + scale_y_continuous(trans = "sqrt") +
    labs(x = NULL, y = "sqrt(Set size)", fill = NULL)
  
  # Legend
  get_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  
  p2 <- get_legend(p1)
  
  # # Overlap sizes
  my_overlap_sizes <- comb_size(comb_mat) %>%
    as.data.frame() %>%
    rename(overlap_sizes = ".") %>%
    mutate(category = row.names(.))
  
  # p3 <- my_overlap_sizes %>%
  #   mutate(category = reorder(category, -overlap_sizes)) %>%
  #   ggplot(aes(x = category, y = overlap_sizes)) +
  #   geom_bar(stat = "identity", fill = "grey80", color = NA, alpha = 0.8, width = 0.7) +
  #   geom_text(aes(label = overlap_sizes, y = 0), size = 5, hjust = 0, vjust = 0.5, angle = 90) +
  #   labs(y = "Intersect size", x = NULL) +
  #   theme_classic() +
  #   theme(text = element_text(size = 14, color = "black"),
  #         axis.text = element_blank(),
  #         axis.ticks.x = element_blank(),
  #         axis.title.x = element_text(hjust = 0))
  
  # Calculate the percentages
  my_overlap_sizes <- my_overlap_sizes %>%
    mutate(overlap_sizes_percentage = (overlap_sizes / sum(overlap_sizes)) * 100)
  
  # Create the ggplot barplot with percentages
  p3 <- my_overlap_sizes %>%
    mutate(category = reorder(category, -overlap_sizes_percentage)) %>%
    ggplot(aes(x = category, y = overlap_sizes_percentage)) +
    geom_bar(stat = "identity", fill = "grey80", color = NA, alpha = 0.8, width = 0.7) +
    # geom_text(aes(label = paste(round(overlap_sizes_percentage, 1), "%")), size = 5, hjust = 0, vjust = 0.5, angle = 90) +
    geom_text(aes(label = paste(round(overlap_sizes_percentage, 1), "%")), 
              size = 5, angle = 90, hjust = 0, y = 1) +
    labs(y = "Intersect size (%)", x = NULL) +
    theme_classic() +
    theme(text = element_text(size = 14, color = "black"),
          axis.text = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(hjust = 0))
  
  if(do.sqr) p3 = p3 + scale_y_continuous(trans = "sqrt") +
    labs(x = NULL, y = "sqrt(Intersect size)", fill = NULL)
  
  p3 = p3 + scale_x_discrete(limits = rev(levels(p3$data$category)))
  
  # Overlap matrix
  my_overlap_matrix <- str_split(string = my_overlap_sizes$category, pattern = "", simplify = T) %>%
    as.data.frame()
  
  colnames(my_overlap_matrix) <- my_names
  
  my_overlap_matrix_tidy <- my_overlap_matrix %>%
    cbind(category = my_overlap_sizes$category) %>%
    pivot_longer(cols = !category, names_to = "Set", values_to = "value") %>%
    full_join(my_overlap_sizes, by = "category") %>%
    full_join(my_set_sizes, by = "Set")
  
  p4 <- my_overlap_matrix_tidy %>%
    mutate(category = reorder(category, -overlap_sizes)) %>%
    mutate(Set = reorder(Set, sizes)) %>%
    ggplot(aes(x = Set, y = category)) +
    geom_tile(aes(fill = Set, alpha = value), color = "grey30", linewidth = 1) +
    scale_fill_manual(values = brewer.pal(4, "Set2"), limits = my_names) +
    scale_alpha_manual(values = c(0.8, 0), limits = c("1", "0")) +
    labs(x = "Sets", y = "Overlap") +
    theme_minimal() +
    theme(legend.position = "none", text = element_text(color = "black", size = 14),
          panel.grid = element_blank(), axis.text = element_blank())
  
  # Put them together
  upset_plot <- wrap_plots(p1, p2, p4, p3,
                           nrow = plot.nrow,
                           ncol = plot.ncol,
                           heights = plot.heights,
                           widths = plot.widths,
                           guides = "collect") &
    theme(legend.position = "none")
  
  
  print(upset_plot)
  
  # Save the plot as an image
  if(!is.null(output_filename)) ggsave(output_filename, 
                                       plot = upset_plot, 
                                       height = save.height, width = save.width, bg = "white", 
                                       device = image_format)
}


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


#' @title frac_to_numeric
#' @description Converts fraction to numeric
#' @return A vector
#' @keywords numeric
#' @export
frac_to_numeric <- function(x) sapply(x, function(x) eval(parse(text=x)))



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

#' @title rotate_df
#' @description rotates a 2D df like tsne UMAP etc any scatter 
#' @param tsne, a 2D df scatter sctructure
#' @param angle, rotation angle
#' @return rotated tsne Y
#' @export
rotate_df <- function(df, angle){
  angle <- (-angle * pi) / (-180)
  rotm <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2)
  df %*% rotm
}
