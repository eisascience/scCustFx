


#' Plot k-Nearest Neighbors Graph
#'
#' This function takes a dataframe that contains x and y coordinates and a label
#' that defines the classes and makes a plot of the k-nearest neighbors graph using 
#' Euclidean distance and the cover tree algorithm.
#'
#' @param df A dataframe with x and y coordinates and a label column defining the classes.
#' @param k An integer that defines the number of nearest neighbors to consider. Default is 1.
#' @param labels An optional character vector specifying the label column in the dataframe. 
#'               If not specified, it is assumed that the label column is named "label".
#' 
#' @return A plot of the k-nearest neighbors graph using ggplot2.
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggplot2 geom_point geom_segment scale_color_discrete theme_classic
#' @importFrom FNN get.knn
#'
#' @examples
#' # create sample data
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' label <- sample(c("A", "B", "C"), 100, replace = TRUE)
#' df <- data.frame(x = x, y = y, label = label)
#' 
#' # plot with k = 2
#' plot_knn(df, k = 2)
#'
#' @export
plot_knn <- function(df, k = 1, labels=NULL) {
  
  
  library(igraph)
  library(ggplot2)
  library(FNN)
  
  
  # compute the k-nearest neighbors using Euclidean distance and cover tree algorithm
  knn <- FNN::get.knn(data = df[, c("x", "y")], k = k, algorithm = "cover_tree")
  
  # get the labels and coordinates
  if(is.null(labels)) labels <- df$label
  coords <- data.frame(x = df$x, y = df$y)
  
  # extract the nearest neighbor indices
  nn_index <- knn$nn.index
  
  # create a data frame with edges between nearest neighbors
  edges <- data.frame(from = nn_index[, 1], to = nn_index[, 2])
  
  # plot the graph with ggplot2
  ggplot(coords, aes(x = x, y = y)) +
    geom_point(aes(color = factor(labels))) +
    geom_segment(data = edges, aes(x = coords[from, "x"], y = coords[from, "y"], 
                                   xend = coords[to, "x"], yend = coords[to, "y"])) +
    theme_classic() +
    scale_color_discrete(name = "Label")
}





#' Make a GIF from plot list
#'
#' This function takes in a list of ggplots and makes a GIF
#'
#' @param count_matrix a count matrix 
#' @param method dist or cor
#' @return plots a network graph
#' @export
plot_gene_network <- function(count_matrix, method = "dist") {
  library(igraph)
  
  
  
  # sort( rowSums(count_matrix), decreasing = T)
  n <- nrow(count_matrix)
  m <- ncol(count_matrix)
  
  if(method=="cor"){
    ajm = round(asinh(cor(count_matrix)^4 * 10), 2)
    
    graph <- graph_from_adjacency_matrix(ajm, mode = "undirected", weighted = TRUE)
    
    
  }
  
  if(method=="dist"){
    
    
    # create a similarity matrix from the count matrix based on Euclidean distance
    distance_matrix <- dist(t(count_matrix))
    similarity_matrix <- 1 / (1 + as.matrix(distance_matrix)^2)
    
    ajm = round(asinh((similarity_matrix) * 10) , 2)
    
    # create an igraph object from the similarity matrix
    graph <- graph_from_adjacency_matrix(, mode = "undirected", weighted = TRUE)
    
  }
  # plot(graph)
  
  # calculate degree centrality and betweenness centrality
  degree <- degree(graph)
  betweenness <- betweenness(graph)
  
  # set node attributes based on centrality measures
  V(graph)$label <- V(graph)$name
  V(graph)$size <- degree 
  V(graph)$color <- "lightblue"
  V(graph)$frame.color <- "white"
  V(graph)$label.color <- "black"
  V(graph)$label.cex <- 0.8
  V(graph)$label.family <- "Helvetica"
  V(graph)$label.dist <- 0
  V(graph)$label.degree <- 0
  
  # set edge attributes based on centrality measures
  E(graph)$color <- "gray"
  # E(graph)$width <- betweenness# * 2
  
  # set edge label attribute with edge weight
  E(graph)$label <- E(graph)$weight
  E(graph)$label.cex <- 0.8
  
  # plot the graph with layout and style adjustments
  plot(graph,
       layout = layout_with_fr,
       vertex.label.color = V(graph)$label.color,
       vertex.label.family = V(graph)$label.family,
       vertex.label.cex = V(graph)$label.cex,
       vertex.label.dist = V(graph)$label.dist,
       vertex.label.degree = V(graph)$label.degree,
       vertex.color = V(graph)$color,
       vertex.frame.color = V(graph)$frame.color,
       vertex.size = V(graph)$size,
       edge.color = E(graph)$color,
       edge.width = E(graph)$width,
       edge.label = E(graph)$label,
       edge.label.cex = E(graph)$label.cex
  )
  
  
}






