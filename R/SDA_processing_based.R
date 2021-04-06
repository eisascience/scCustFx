#' @title SDA.AddCompStats
#'
#' @description SDA.AddCompStats
#' @param SDAres An SDA results list/object 
#' @return An SDA results list/object with Comp statistics
#' @import data.table
#' @export
SDA.AddCompStats <- function(SDAres, sdThrsh = 0.04, maxscoreThrsh=20, 
                             maxloadThrsh = 1, redoCalc=T){
  if (redoCalc)   SDAres$component_statistics <- NULL
  if (is.null(SDAres$component_statistics) ) {
    
    SDAres$component_statistics <- as.data.frame(data.table(
      Component = 1:SDAres$n$components, 
      Component_name = dimnames(SDAres$scores)[[2]], 
      max_score = apply(abs(SDAres$scores),  2, max),
      max_loading = apply(abs(SDAres$loadings[[1]]), 1, max),
      mean_score = apply(abs(SDAres$scores),  2, mean),
      mean_loading = apply(abs(SDAres$loadings[[1]]), 1, mean),
      sd_score = apply(abs(SDAres$scores),  2, sd),
      sd_loading = apply(abs(SDAres$loadings[[1]]), 1, sd),
      ssqrd_score = apply(SDAres$scores,  2, function(x) sum(x^2)),
      ssqrd_loading = apply(SDAres$loadings[[1]], 1, function(x) sum(x^2))
    )[order(-Component)])
    
  }
  
  SDAres$component_statistics$Component_namev2    <- SDAres$component_statistics$Component_name
  SDAres$component_statistics$Component_name_plot <- SDAres$component_statistics$Component_name
  
  SDAres$component_statistics[which(SDAres$component_statistics$sd_loading<sdThrsh),]$Component_namev2         <- rep("", length(which(SDAres$component_statistics$sd_loading<sdThrsh)))
  
  SDAres$component_statistics[which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh),]$Component_name_plot <- rep("", length(which(SDAres$component_statistics$max_score<maxscoreThrsh & SDAres$component_statistics$max_loading<maxloadThrsh)))
  
  SDAres$component_statistics <- data.table(SDAres$component_statistics)
  return(SDAres)
  
}