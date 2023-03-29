

#' Plot a word cloud from a named vector of word counts
#'
#' @param word_count A named vector where the names are words and the values are counts
#' @param max.words The maximum number of words to display in the word cloud (default: 100)
#' @param colors The color palette to use for the word cloud (default: brewer.pal(8, "Dark2"))
#' @param min.freq The minimum frequency of a word to be included in the word cloud (default: 1)
#'
#' @return A plot of the word cloud
#'
#' @importFrom wordcloud wordcloud
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_wordcloud <- function(word_count, max.words = 100, colors = brewer.pal(8, "Dark2"), min.freq = 1) {
  library(wordcloud)
  # Convert named vector to data frame
  word_df <- data.frame(word = names(word_count), freq = as.integer(word_count))
  
  # Plot word cloud
  wordcloud(words = word_df$word, freq = word_df$freq, min.freq = min.freq,
            max.words = max.words, random.order = FALSE, rot.per = 0.35,
            colors = brewer.pal(8, "Dark2"))
}






#' Make a GIF from plot list
#'
#' This function takes in a list of ggplots and makes a GIF
#'
#' @param plot_list list of ggplots 
#' @param filename path and filename of gif
#' @return saves filename gif
#' @export
create_gif_LS <- function(plot_list, filename = "animation.gif", delay = 100, res = 72) {
  library(animation)
  n_frames <- length(plot_list)
  ani.options(interval = delay/1000, nmax = n_frames, 
              ani.width = 600, ani.height = 400, 
              outdir = getwd(), 
              title = "Animation", description = "My animation")
  
  saveGIF({
    for (i in 1:n_frames) {
      print(plot_list[[i]])
    }
  }, movie.name = filename, res = res, loop = TRUE)
}



#' Add quadrant labels to a ggplot representing a scatterplot
#'
#' This function takes in a ggplot that is representing a scatterplot, a horizontal y coordinate, and a vertical x coordinate, and adds quadrant boundaries and text labels to the plot to indicate how many points are in each section of the scatterplot.
#'
#' @param plot ggplot object representing a scatterplot
#' @param x horizontal y coordinate to separate the quadrants
#' @param y vertical x coordinate to separate the quadrants
#'
#' @return ggplot object with quadrant boundaries and text labels added
#' @export
add_quadrant_labels <- function(plot, x, y) {
  # Get data from plot
  data <- layer_data(plot, 1)
  
  # Calculate quadrant boundaries
  x_median <- x
  y_median <- y
  
  # Count number of points in each quadrant
  q1_count <- sum(data$x > x_median & data$y > y_median)
  q2_count <- sum(data$x <= x_median & data$y > y_median)
  q3_count <- sum(data$x <= x_median & data$y <= y_median)
  q4_count <- sum(data$x > x_median & data$y <= y_median)
  
  # Calculate percentages
  total_count <- nrow(data)
  q1_percent <- round(q1_count / total_count * 100)
  q2_percent <- round(q2_count / total_count * 100)
  q3_percent <- round(q3_count / total_count * 100)
  q4_percent <- round(q4_count / total_count * 100)
  
  # Get min and max values of x and y coordinates
  maxX <- max(data$x)
  minX <- min(data$x)
  maxY <- max(data$y)
  minY <- min(data$y)
  
  # Calculate position of text labels
  q1_x <- mean(c(x_median, maxX))
  q1_y <- mean(c(y_median, maxY))
  
  q2_x <- mean(c(minX, x_median))
  q2_y <- mean(c(y_median, maxY))
  
  q3_x <- mean(c(minX, x_median))
  q3_y <- mean(c(minY, y_median))
  
  q4_x <- mean(c(x_median, maxX))
  q4_y <- mean(c(minY, y_median))
  
  # Add lines for quadrant boundaries
  plot <- plot +
    geom_vline(xintercept = x_median, linetype = "dashed") +
    geom_hline(yintercept = y_median, linetype = "dashed")
  
  # Add text labels for each quadrant
  plot <- plot +
    annotate("text", x = q1_x, y = q1_y, 
             label = paste0(q1_count, " (", q1_percent, "%)"), size = 4) +
    annotate("text", x = q2_x, y = q2_y, 
             label = paste0(q2_count, " (", q2_percent, "%)"), size = 4) +
    annotate("text", x = q3_x, y = q3_y, 
             label = paste0(q3_count, " (", q3_percent, "%)"), size = 4) +
    annotate("text", x = q4_x, y = q4_y, 
             label = paste0(q4_count, " (", q4_percent, "%)"), size = 4)
  
  return(plot)
}




#' barplot vectors side by side 
#'
#' barplot vectors side by side
#' 
#' @param vec1 a named vector 
#' @param vec2 a secibd named vector 
#' @param title title
#' @param QuantThr color pink via quantiles
#' @return If returnDF is `TRUE`, the plotted data frame is returned.
#' @export
barplot_vectors <- function(vec1, vec2, title="", QuantThr = .66) {
  # Convert vectors to data frames
  # CellMembrane:::.UpdateGeneModel()
  names(vec1) =  CellMembrane::RenameUsingCD(CellMembrane:::.UpdateGeneModel(names(vec1)))
  names(vec2) =  CellMembrane::RenameUsingCD(CellMembrane:::.UpdateGeneModel(names(vec2)))
  
  df1 <- data.frame(names = names(vec1), value = vec1, stringsAsFactors = FALSE)
  df2 <- data.frame(names = names(vec2), value = vec2, stringsAsFactors = FALSE)
  df1$col = "gray"
  df2$col = "grey"
  df1 = df1[order(abs(df1$value), decreasing = T),]
  df2 = df2[order(abs(df2$value), decreasing = T),]
  
  
  df1[abs(df1$value) > quantile(abs(df1$value), .66),]$col = "pink"
  df2[abs(df2$value) > quantile(abs(df2$value), .65),]$col = "pink"
  
  # if(nrow(df2)>N_highlight){
  #   df2[1:N_highlight, ]$col = "pink"
  # }
  # 
  # if(nrow(df1)>N_highlight){
  #   df1[1:N_highlight, ]$col = "pink"
  # }
  
  df1$col <- factor(df1$col, levels = c("pink", "gray"))
  df2$col <- factor(df2$col, levels = c("pink", "gray"))
  
  
  
  plot1 = ggplot(df1, aes(x = reorder(names, value), y = value, fill=col)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    geom_text(aes(label = names), position = position_stack(vjust = 0.5), color = "black") +
    coord_flip() +
    labs(x = NULL, y = NULL, title = title)+
    theme_classic() +
    theme(legend.position = "none",
          # axis.line = element_blank(),
          # axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_blank())+
    scale_fill_manual(values = c("pink", "gray"))
  
  
  
  
  
  plot2 = ggplot(df2, aes(x = reorder(names, -value), y = value, fill=col)) + 
    geom_bar(stat = "identity") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
    geom_text(aes(label = names), position = position_stack(vjust = 0.5), color = "black") +
    coord_flip() +
    labs(x = NULL, y = NULL, title = title)+
    theme_classic() +
    theme(legend.position = "none",
          # axis.line = element_blank(),
          axis.line.y = element_blank(),
          # axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_blank())+
    scale_fill_manual(values = c("pink", "gray"))
  
  
  
  wrap_plots(plot2, plot1, ncol = 2)
  
  # cowplot::plot_grid(plot1, plot2, title = title, ncol = 2)
  
  # Combine plots using patchwork
  
  
}



#' Plot 2D Dimensionality Reduction
#'
#' Plot 2D dimensional reduction of a data frame.
#' 
#' @param plotDF A data frame of numeric values to plot.
#' @param centers Number of centers to use for k-means clustering.
#' @param NfeatPerClus Number of features to include in the plot per cluster.
#' @param returnDF Logical indicating if the plotted data frame should be returned.
#' @param title Character string for the plot title.
#' @return If returnDF is `TRUE`, the plotted data frame is returned.
#' @export
plot2Dredux <- function(plotDF, centers = 8, NfeatPerClus = 10, returnDF = F, title = "", base_size = 12){
  
  colnames(plotDF) = c("Dim1", "Dim2")
  plotDF$labs1 = rownames(plotDF)
  
  plotDF$km = factor(kmeans(plotDF[,1:2] , centers = centers)$cluster)
  
  # group by km
  keepLabs <- plotDF %>% group_by(km) %>% 
    sample_n(NfeatPerClus, replace = F) %>% 
    pull(labs1)
  
  plotDF$labs2 = ""
  
  plotDF[keepLabs, ]$labs2 = keepLabs
  
  
  
  
  # ggplot(plotDF, aes(tSNE1, tSNE2, col=sums)) +
  #   geom_point(size = 1) + theme_bw() + 
  #   scale_color_distiller(palette = "Spectral") +
  #   theme(legend.position = "bottom") +
  #   ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ")
  
  gg = ggplot(plotDF, aes(Dim1, Dim2, col=km)) +
    geom_point(size = 1) + ggplot2::theme_bw(base_size = base_size) +
    # scale_color_distiller(palette = "Spectral") +
    theme(legend.position = "bottom") +
    ggrepel::geom_text_repel(aes(label=labs2), size=3, max.overlaps =150, segment.color="grey50", nudge_y = 0.25) +
    ggtitle(title) + 
    scale_color_manual(values = scCustFx:::ColorTheme()$col_vector)
  
  
  # ggplot(subset(plotDF, km==6), aes(tSNE1, tSNE2, col=km)) +
  #   geom_point(size = 1) + theme_bw() +
  #   # scale_color_distiller(palette = "Spectral") +
  #   theme(legend.position = "bottom") +
  #   ggrepel::geom_text_repel(aes(label=labs1), size=3, max.overlaps =150, segment.color="grey50", nudge_y = 0.25) +
  #   ggtitle("tSNE SDA qc Components\n Sum absolute-cell-scores normalized by its mean \n ") + 
  #   scale_color_manual(values = scCustFx:::ColorTheme()$col_vector)
  # 
  if(returnDF) return(plotDF) else return(gg)
  
  
}






