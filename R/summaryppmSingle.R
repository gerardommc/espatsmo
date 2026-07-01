#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmSingle, if the argument as.ppmSingle = TRUE, otherwise 
#' summary will call spatstat.model::summary.ppm
#' @param object A model fitted with ppmSingle
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
#' summary(model)
#' @export
#' @method summary ppmSingle

summary.ppmSingle <- function(object = NULL){
  spatstat.model::summary.ppm(object = object$model)
}


