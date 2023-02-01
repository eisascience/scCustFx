

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

