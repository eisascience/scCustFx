#' Create a QQ Plot
#'
#' This function generates a quantile-quantile (QQ) plot using ggplot2.
#'
#' @param x Not used (for flexibility).
#' @param y The sample data for which you want to create the QQ plot.
#' @param x_title The title for the x-axis.
#' @param y_title The title for the y-axis.
#' @param title The title for the plot.
#'
#' @return A ggplot object representing the QQ plot.
#'
#' @examples
#' set.seed(123)
#' sample_data <- rnorm(100)
#' qq_plot <- create_qq_plot(NULL, sample_data, "Custom X-axis Title", "Custom Y-axis Title", "Custom QQ Plot Title")
#' print(qq_plot)
#'
#' @export
create_qq_plot <- function(x, y, x_title = "Theoretical Quantiles", 
                           y_title = "Sample Quantiles", title = "QQ Plot", theme = theme_bw(base_size = 14)) {
  
  gg_qqplot <- ggplot(data = data.frame(x = x, y = y),
                      aes(x = x, y = y)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = title,
         x = x_title,
         y = y_title) + theme
  
  return(gg_qqplot)
}


#' Save a patchwork of ggplot figures to a file
#'
#' This function takes a list of ggplot figures and saves them as a patchwork to a file specified by a base path and filename, split into multiple files if the list is longer than a given split N.
#'
#' @param ggLS A list of ggplot figures.
#' @param basepath The base path for the output file.
#' @param filename The name of the output file, not including any extension.
#' @param splitN The number of ggplots to include in each patchwork file (default is 20).
#' @param width The width of the output file in inches (default is 12).
#' @param height The height of the output file in inches (default is 7).
#' @param units The units of the output file (default is "in").
#' @param dpi The resolution of the output file in dots per inch (default is 150).
#'
#' @importFrom patchwork wrap_plots
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' ## Save a list of ggplot figures to a file
#' my_ggplots <- list(ggplot(data = mtcars, aes(x = mpg, y = wt)) + geom_point(),
#'                    ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width)) + geom_point())
#' save_ggplot_list(ggLS = my_ggplots, basepath = "/path/to/file", filename = "my_ggplots")
#'
#' @export
save_ggplot_list = function(ggLS, basepath, filename, splitN = 20, width = 12, height = 7, units = "in", dpi = 150) {
  n_plots <- length(ggLS)
  n_patches <- ceiling(n_plots / splitN)
  
  for (i in 1:n_patches) {
    start_idx <- (i - 1) * splitN + 1
    end_idx <- min(i * splitN, n_plots)
    patch <- patchwork::wrap_plots(ggLS[start_idx:end_idx], ncol = 5)
    patch_filename <- file.path(basepath, sprintf("%s_%02d.png", filename, i))
    ggsave(filename = patch_filename, plot = patch, width = width, height = height, units = units, dpi = dpi)
  }
}





#' @title PlotHistDensOptima
#' @description This function plots histogram+density and finds all optima from a df usually from FetchData() of Seurat. 
#' @param dftemp, a dataframe with column var1 to be plotted. 
#' @param col_vector, color vector.
#' @param nn, a color from color_vector.
#' @return plot with ggplot tabulated bar plots. 
#' @export
PlotHistDensOptima <- function(dftemp, Print = F, col_vector = col_vector, nn = NULL, Title = ""){
  dens <- density(dftemp$var1, n = 1000)
  peaks <- dens$x[find_peaks(dens$y, m=20)]
  dips <- dens$x[find_peaks(-1*dens$y, m=20)]
  
  if(is.null(nn)) nn=1
  
  
  p1 <- ggplot(dftemp, aes(x=var1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth = 0.05) +
    geom_density(alpha=.2, fill=col_vector[nn]) +
    geom_vline(aes(xintercept=mean(var1)), color="blue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=mean(var1) - 2*sd(var1)), color="dodgerblue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=mean(var1) + 2*sd(var1)), color="dodgerblue", linetype="dashed", size=.5) +
    geom_vline(aes(xintercept=peaks[1] ), color="plum", size=.8) +
    geom_vline(aes(xintercept=dips[1] ), color="red", size=1)  +
    theme_bw() + ggtitle(Title)
  
  
  if((length(peaks) >= 2)) p1 <- p1 + geom_vline(aes(xintercept=peaks[2] ), color="plum", size=.8)
  if((length(peaks) >= 3)) p1 <- p1 + geom_vline(aes(xintercept=peaks[3] ), color="plum", size=.8)
  if((length(peaks) >= 4)) p1 <- p1 + geom_vline(aes(xintercept=peaks[4] ), color="plum", size=.8)
  
  if((length(dips) >= 2)) p1 <- p1 + geom_vline(aes(xintercept=dips[2] ), color="brown", size=.8)
  if((length(dips) >= 3)) p1 <- p1 + geom_vline(aes(xintercept=dips[3] ), color="brown", size=.8)
  if((length(dips) >= 4)) p1 <- p1 + geom_vline(aes(xintercept=dips[4] ), color="brown", size=.8)
  
  if(Print) print(p1) else return(p1)
}



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

