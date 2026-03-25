#' @title Compute and print the AIC of single or batch point process models
#' @description
#' Print the Akaike information criterion of a point process models fitted with spatstat's ppm, but via the inteface of the presente package.
#' @param object A ppmSingle, fitted with the corresponding function
#' @return The function returns the Akaike Information Criterion, one for a ppmSigle object, or the values of all the returned mdels with ppmBatchFit
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' model <- ppmSingleFit(points= p, 
#'                      covariates = r, 
#'                      formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                      bias.data = bias, #Data frame with sampling localities or raster layer
#'                      bias.correction = "weights",
#'                      as.ppmSingle = F)
#' 
#' AIC(model)
#' @export
#' @method AIC ppmSingle


AIC.ppmSingle <- function(object){
  spatstat.model::AIC.ppm(object)
}


