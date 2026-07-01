#' @title Compute and print the log-Likelihood of single or batch point process models
#' @description
#' Print the log-Likelihood of a point process models fitted with spatstat's ppm, but via the espatsmo interface.
#' @param object A ppmSingle, fitted with the corresponding function
#' @return The function returns the log-Likelihood for a ppmSigle object.
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
#' logLik(model)
#' @export
#' @method logLik ppmSingle


logLik.ppmSingle <- function(object){
  stats::logLik(object$model)
}


