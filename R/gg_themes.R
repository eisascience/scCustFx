

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
theme_simple2 <- theme_classic(base_size = 20) +
                  theme(legend.position = "none",
                        axis.line = element_blank(),
                        axis.text.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        plot.title = element_blank())

