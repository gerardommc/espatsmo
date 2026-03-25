#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmSingle, if the argument as.ppmSingle = T, otherwise 
#' summary will call spatstat.model::summary.ppm
#' @param object A model fitted with ppmSingle
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |>  scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast("")
#' 
#' model <- ppmSingleFit(points= p, 
#'                      covariates = r, 
#'                      formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                      bias.data = bias, #Data frame with sampling localities or raster layer
#'                      bias.correction = "weights",
#'                      as.ppmSingle = F)
#' 
#' summary(model)
#' @export
#' @method summary ppmSingle

summary.ppmSingle <- function(object = NULL){
  spatstat.model::summary.ppm(object = object$model)
}


