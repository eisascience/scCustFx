

#' Theme Publication
#'
#' This function defines a custom theme for publication-quality plots.
#'
#' @param base_size The base font size for the plot elements.
#' @param base_family The base font family for the plot elements.
#' @export
theme_Publication <- function(base_size=14, base_family="") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ))
  
}


#' Scale Colour Publication
#'
#' This function sets the color scale to a predefined set of colors commonly used for publication purposes.
#'
#' @param ... Additional arguments to be passed to the discrete_scale function.
#'
#' @export
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#' Scale Fill Publication
#'
#' This function sets the fill color scale to a predefined set of colors commonly used for publication purposes.
#'
#' @param ... Additional arguments to be passed to the discrete_scale function.
#'
#' @export
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}




#' @title theme_simple remove legend, and axis text
#'
#' @description ggplot theme, use as + theme_simple()
#' @return A ggplot theme
#' @export
theme_simple <- ggplot2::theme(legend.position = "none",
                           axis.title.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.ticks = ggplot2::element_blank())


#' @title theme_simple remove legend, axis, and title with base size 20
#'
#' @description ggplot theme, use as + theme_simple()
#' @return A ggplot theme
#' @export
theme_simple2 <- ggplot2::theme_classic(base_size = 20) +
                 ggplot2::theme(legend.position = "none",
                                axis.line = ggplot2::element_blank(),
                                axis.text.x = ggplot2::element_blank(),
                                axis.text.y = ggplot2::element_blank(),
                                axis.ticks = ggplot2::element_blank(),
                                axis.title = ggplot2::element_blank(),
                                plot.title = ggplot2::element_blank())

