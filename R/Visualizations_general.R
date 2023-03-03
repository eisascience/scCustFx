

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
  
  plotDF$km = factor(kmeans(plotDF[,1:2] , centers = 8)$cluster)
  
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




#' Plot a pie and donut chart
#' 
#' @param data Data to plot.
#' @param mapping Mapping for colors.
#' @param start Starting angle for pie chart.
#' @param addPieLabel Add labels to pie chart.
#' @param addDonutLabel Add labels to donut chart.
#' @param showRatioDonut Show ratios in donut chart.
#' @param showRatioPie Show ratios in pie chart.
#' @param ratioByGroup Show ratio by group.
#' @param showRatioThreshold Threshold for showing ratios.
#' @param labelposition Label position.
#' @param labelpositionThreshold Threshold for label position.
#' @param r0 Inner radius of donut chart.
#' @param explode Explode pie chart.
#' @param selected Selected data to highlight.
#' @param explodePos Explode position.
#' @param color Color for pie chart.
#' @param pieAlpha Alpha for pie chart.
#' @param donutAlpha Alpha for donut chart.
#' @param maxx Maximum value for x axis.
#' @param showPieName Show name for pie chart.
#' @param showDonutName Show name for donut chart.
#' @param title Plot title.
#' @param pieLabelSize Font size for pie chart labels.
#' @param donutLabelSize Font size for donut chart labels.
#' @param titlesize Font size for plot title.
#' @param explodePie Explode pie chart.
#' @param explodeDonut Explode donut chart.
#' @param use.label Use label.
#' @param use.labels Use labels.
#' @param family Font family.
#' @return A pie and donut chart.
#' @export
PieDonut2 <- function(data, mapping, start = getOption("PieDonut.start", 0),
                      addPieLabel = TRUE, addDonutLabel = TRUE, showRatioDonut = TRUE, 
                      showRatioPie = TRUE, ratioByGroup = TRUE, 
                      showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.02), 
                      labelposition = getOption("PieDonut.labelposition", 2), 
                      labelpositionThreshold = 0.1, 
                      r0 = getOption("PieDonut.r0",                                                                                                                                                                                                                   1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                      color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                      showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                      pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                      explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                      family = getOption("PieDonut.family", "")){
  require(moonBook)
  require(ggforce)
  require(grid)
  
  (cols = colnames(data))
  if (use.labels) 
    data = addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = getMapping(mapping, "x"))
  (donuts = getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    # df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
    #                                                          "\n(", scales::percent(df$ratio, accuracy=.1), ")"), as.character(df$label))
    
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                             " (", scales::percent(df$ratio, accuracy=.1), ")"), as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  # mainCol = gg_color_hue(nrow(df))
  mainCol = colorRampPalette(RColorBrewer::brewer.pal(11, 'Spectral'))(nrow(df))
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(round(df3$Freq/df3$group, 3), accuracy=.1)
    }
    else {
      df3$ratio <- scales::percent(round(df3$ratio1, 3), accuracy=.1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$ratio, 
                           ")")
      else {
        # df3$label = paste0(df3$label, "\n(", df3$ratio, 
        #                    ")")
        df3$label = paste0(df3$label, " (", df3$ratio, 
                           ")")
      }
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  ####################
  # df3$angle = 30 #seq(0, 360, nrow(df3))
  # myAng <-seq(-20,-340, length.out = nrow(df3))
  
  p <- ggplot() + theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", r0 = as.character(r0), 
                                    r = as.character(r1), start = "start1", end = "end1", 
                                    fill = pies), alpha = pieAlpha, color = color, data = df) + 
    transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = "none")
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", vjust = "vjust"), 
                size = pieLabelSize, data = df, family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", y = "segy", 
                                       xend = "segxend", yend = "segyend"), data = df[df$ratio < 
                                                                                        labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                                                                                                                                          y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                          vjust = "vjust"), size = pieLabelSize, data = df[df$ratio < 
                                                                                                                                                                                             labelpositionThreshold, ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                                                                                                                                                                                                y = "labely", label = "label"), size = pieLabelSize, 
                                                                                                                                                                                                                                                     data = df[df$ratio >= labelpositionThreshold, ], 
                                                                                                                                                                                                                                                     family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no", 
                                        explode = "focus"), alpha = donutAlpha, color = color, 
                             data = df3)
    }
    else {
      p3 <- p + geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                        r0 = as.character(r1), r = as.character(r2), 
                                        start = "start1", end = "end1", fill = "no"), 
                             alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = "none")
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3) + 
        geom_text(aes_string(x = "segxend", y = "segyend", 
                             label = "label", hjust = "hjust", vjust = "vjust"), 
                  size = donutLabelSize, data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", y = "labely", 
                                      label = "label"), size = donutLabelSize, data = df3, 
                           family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", y = "segy", 
                                         xend = "segxend", yend = "segyend"), data = df3[df3$ratio1 < 
                                                                                           labelpositionThreshold, ]) + geom_text(aes_string(x = "segxend", 
                                                                                                                                             y = "segyend", label = "label", hjust = "hjust", 
                                                                                                                                             vjust = "vjust"), size = donutLabelSize, data = df3[df3$ratio1 < 
                                                                                                                                                                                                   labelpositionThreshold, ], family = family) + 
        geom_text(aes_string(x = "labelx", y = "labely", 
                             label = "label"), size = donutLabelSize, data = df3[df3$ratio1 >= 
                                                                                   labelpositionThreshold, ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, label = title, 
                          size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, y = r3, 
                          label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid.newpage()
    print(p1, vp = viewport(height = 1, width = 1))
    print(p3, vp = viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}




#' @title heatmap_table
#'
#' @description Takes a table and plots it. Either as a table or a pheatmap heatmap with clustering
#' @param TBL A table (base R table()) 
#' @param clust.meth default is "ward.D2", see pheatmap for more options
#' @param clust.dist default is "euclidean", see pheatmap for more options
#' @param Title A title else ''
#' @param plotTable default is F, if T, it will plot an actual table
#' @return A pheatmap heatmap or a plotted table
#' @export
heatmap_table <- function(TBL, clust.meth="ward.D2", clust.dist="euclidean", Title="", plotTable=F){
  if(plotTable){
    grid.table(TBL)
  } else {
    pheatmap::pheatmap(as.matrix(TBL),
                       cluster_rows = T, cluster_cols = T,
                       main = paste0(Title, "\n", clust.dist, " dist w. ", clust.meth, "clustering"),
                       color = viridis(20), cellwidth = 10, show_colnames = T,
                       clustering_distance_rows = clust.dist,
                       clustering_distance_cols = clust.dist, #euclidean
                       clustering_method = clust.meth,
                       # annotation_col = ColAnnRow,
                       # annotation_colors = my_colour
                       #, labels_row = IDgenes
    )
  }
}

#' @title gg_barplot_table, aka OOSAP::PlotMyTable
#'
#' @description Takes a table and plots a gg barplot 
#' @param MyTable A 2way table (base R table()) 
#' @param Title A title else ''
#' @param legend.position legend position defautl "bottom"
#' @param PlotCombo if T plots % and count plots together
#' @param xlab x axis label
#' @param xtext_angle x axis label
#' @param hjust if angle is 90 set to 1 
#' @param vjust if angle is 90 set to 0.5
#' @param returnlist logical
#' @param base_size numerical
#' @return A ggplot 
#' @export
gg_barplot_2waytable = function(MyTable, Title="", legend.position="bottom", PlotCombo = F, xlab="", 
                                xtext_angle = 45, hjust = 0.9,vjust = 0.8,base_size = 10,
                                returnlist = T){
  tempMeltDF <- reshape2::melt(MyTable)
  
  
  p1 <- (ggplot(tempMeltDF) +
           geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7) +
           theme_bw(base_size = base_size)  + scale_fill_manual(values=col_vector) +
           theme(legend.position=legend.position,
                 legend.direction="horizontal",
                 axis.text.x = ggplot2::element_text(angle = xtext_angle, vjust = vjust, hjust=hjust),
                 legend.title = element_blank()) +
           ggtitle(paste0(Title, "\n Total Contribution")) + ylab("Total No.")) + xlab(xlab)
  
  
  p2<- (ggplot(tempMeltDF) +
          geom_bar(aes(x=Var2, y=value, fill=factor(Var1)), stat="identity", width = 0.7, position="fill") +
        theme_bw(base_size = base_size)  + scale_fill_manual(values=col_vector) +
          theme(legend.position=legend.position,
                legend.direction="horizontal",
                axis.text.x = ggplot2::element_text(angle = xtext_angle, vjust = vjust, hjust=hjust),
                legend.title = element_blank())+
          scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                             labels = scales::percent(c(0, 0.25, 0.5, 0.75, 1))) +
          ggtitle(paste0(Title, "\n Relative % Contribution")) + ylab("(%)")) + xlab(xlab)
  


  
  if(PlotCombo) {
    print(cowplot::plot_grid(p1, p2, ncol = 1)) 
    p1 + theme(legend.position="none") / p2
    } else {
      if(returnlist) {
        list(p1 = p1, p2=p2)
      } else {
        p1
        p2
      }
    
  }
  
}



#' @title gg_barplot_table
#'
#' @description Takes a table and plots a gg barplot 
#' @param TBL A table (base R table()) 
#' @param Title A title else ''
#' @param gg_theme a ggplot theme or NULL
#' @return A ggplot 
#' @export
gg_barplot_table <- function(TBL, Title = "", gg_theme=NULL){
  
  tbl.m = reshape2::melt(TBL)
  if("Var1" %in% colnames(tbl.m)) {
    if(class(tbl.m$Var1)=="integer") tbl.m$Var1 = as.character(tbl.m$Var1)
    ggp = ggplot(tbl.m, aes(x=Var1, y=value)) 
  }
  
  if("variable" %in% colnames(tbl.m)) {
    if(class(tbl.m$variable)=="integer") tbl.m$variable = as.character(tbl.m$variable)
    ggp = ggplot(tbl.m, aes(x=variable, y=value)) 
  }
  
  if(!is.null(gg_theme)){
    ggp = ggp + gg_theme()
  } else {
    ggp = ggp + ggplot2::theme_bw()
  }
  
  
  ggp   +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggplot2::ggtitle(Title)
}



#' @title gg_tSNE_df
#'
#' @description Takes a data.frame with tSNE1 and tSNE2 precomputed and plots it.
#' @param DF A data.frame 
#' @param SavePath A data.frame 
#' @param Markers A char. vec. that that all must exist in the data.frame.
#' @param splitby a single character which data.frame must contain this value in the cols; the resulting plot is faceted by this.
#' @param verbose if T verbose output is printed 
#' @param Title A title else ''
#' @param gg_theme a ggplot theme or NULL
#' @return A ggplot 
#' @export
gg_tSNE_df <- function(DF, 
                       SavePath=NULL, 
                       Markers=NULL, 
                       splitby=NULL,
                       verbose = T,
                       size=0.2, alpha=0.6,
                       Title="",
                       width = 2000, height = 2300, pointsize = 20, ncol=4, 
                       limits = c(0, 1),
                       gg_theme=NULL, printPlot = T, vridis_col_opt = "A"){
  
  if(! ("tSNE1" %in% colnames(DF)) & "tSNE2" %in% colnames(DF)) stop("tSNE1 and/or tSNE2 not found in df")
  if(is.null(Markers)) stop("Markers is NULL")
  
  ExprPlots <- list()
  
  if(is.null(splitby)){
    
    for(xN in Markers){
      if(verbose) print(xN)
      

      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", xN)], id.vars=c("tSNE1", "tSNE2")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        ggthemes::scale_colour_gradient2_tableau('pi0', low = "blue", mid = "white", high = "red", midpoint = 0)
        scale_color_viridis(option=vridis_col_opt,limits = limits, oob = scales::squish) +
        ggtitle(paste0(Title, "\ntSNE - ", xN))
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  } else {
    
    if(! splitby %in% colnames(DF)) stop("splitby not in colnames of DF")
    
    DF$VarSplit = DF[,splitby]
    
    for(xN in Markers){
      # xN = Markers[1]
      if(verbose) print(xN)

      
      ExprPlots[[xN]] <- ggplot(reshape2::melt(DF[, c("tSNE1", "tSNE2", "VarSplit", xN)], id.vars=c("tSNE1", "tSNE2", "VarSplit")),
                                aes(tSNE1, tSNE2, color=asinh(value))) +
        geom_point(size=size, alpha=alpha)+theme_bw() +
        theme(legend.position = "bottom") +
        scale_color_viridis(option="A") +
        ggtitle(paste0(Title, "\ntSNE - ", xN, "\n per", splitby))+ facet_wrap(~VarSplit)
      
      if(!is.null(gg_theme)) ExprPlots[[xN]] <- ExprPlots[[xN]] + gg_theme
    }
    
  }
  
  
  if(!is.null(SavePath)) {
    
    png(SavePath, width = width, height = height, pointsize = pointsize)
    
    print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol))
    
    dev.off()
    
    
  } else { 
    if(printPlot) {
      return(print(cowplot::plot_grid(plotlist = ExprPlots, ncol=ncol)))
    } else {
      return(ExprPlots)
    } 
    
    }
  
}