

#' @title theme_simple remove duplicated extra axes
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

