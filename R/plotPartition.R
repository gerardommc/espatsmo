#' @title Plot a spatPart
#' @description
#' Plot a spatial partition of a point process model and its associated covariates
#' @param x A spatPart object generated with the spatialPartition function, or a list containing two spatial splitswith class coerced to spatPart
#' @param partition A single numeric value off 0 or 1, or a vector containing both indices, used to select the data partittion to plot 
#' @return A r graphics plot
#' @examples 
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#'
#' part <- spatialPartition(covariates = r,
#'                          points = p, 
#'                          no.blocks = 100,
#'                          part.criteria = "random",
#'                          mask.criteria = "random",
#'                          seed = 432)
#'                          
#' plot(part, partition = c(1, 0))
#' @export
#' @method plot spatPart

plot.spatPart <- function(x, partition = 1){
  
  if(length(partition) == 1){
      if(partition == 0){
        terra::plot(x$part0$covariates[[1]], main = "Partition 0")
        points(x$part0$points, col = "red", pch = 20)
      }
      
      if(partition == 1){
        terra::plot(x$part1$covariates[[1]], main = "Partition 1")
        points(x$part1$points, col = "red", pch = 20)
      }
  }
  
  if(length(partition) > 1){
    par(mfrow = c(1,2))
    
    terra::plot(x$part0$covariates[[1]], main = "Partition 0")
    points(x$part0$points, col = "red", pch = 20)
    
    terra::plot(x$part1$covariates[[1]], main = "Partition 1")
    points(x$part1$points, col = "red", pch = 20)
  }
}
