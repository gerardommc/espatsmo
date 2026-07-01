#' @title Compute and print the BIC of single or batch point process models
#' @description
#' Print the Bayesian information criterion of a point process models fitted with spatstat's ppm, but via the inteface of the presente package.
#' @param object A ppmSingle, fitted with the corresponding function
#' @return The function returns the Bayesian Information Criterion, one for a ppmSigle object, or the values of all the returned mdels with ppmBatchFit
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#' 
#' #Model with random sampling
#' 
#' model <- ppmSingleFit(points = pr,
#'                       covariates = r,
#'                       formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)",
#'                       as.ppmSingle = FALSE)
#' 
#' BIC(model)
#' @export
#' @method BIC ppmSingle


BIC.ppmSingle <- function(object){
  stats::BIC(object$model)
}


