#' @title Compute and print the log-Likelihood of a batch of point process models
#' @description
#' Print the log-Likelihood of one or multiple point proces models fitted with spatstat's ppm, but via the inteface of the presente package.@
#' @param object A ppmBatch object produced by ppmBatchFit
#' @return The function returns the log-Likelihood, one for a ppmSigle object, or the values of all the returned mdels with ppmBatchFit
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#'
#' resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getPolyFormulas(respDF = resp, 
#'                          compatMat = compat)
#' 
#' models <- ppmBatchFit(points = pr,
#'                       covariates = r,
#'                       formulas = forms,
#'                       parallel = FALSE,
#'                       top.models = 3)
#' 
#' logLik(models)
#' @export
#' @method logLik ppmBatch

logLik.ppmBatch <- function(object){
  sapply(object$models, stats::logLik) |> print(quote = FALSE)
}


