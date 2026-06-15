#' @title Compute and print the BIC of a batch of point process models
#' @description
#' Print the Bayesian information criterion of one or multiple point proces models fitted with spatstat's ppm, but via the inteface of the presente package.@
#' @param object A ppmBatch object produced by ppmBatchFit
#' @return The function returns the Bayesian Information Criterion, one for a ppmSigle object, or the values of all the returned mdels with ppmBatchFit
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' resp <- system.file("extdata", "Exponents.csv", package = "espatsmo") |>  read.csv()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getPolyFormulas(respDF = resp, 
#'                          compatMat = compat)
#' 
#' models <- ppmBatchFit(points = p,
#'                       covariates = r,
#'                       formulas = forms,
#'                       bias.data = bias,
#'                       bias.correction = "weights",
#'                       parallel = FALSE,
#'                       top.models = 3)
#' 
#' BIC(models)
#' @export
#' @method BIC ppmBatch

BIC.ppmBatch <- function(object){
  sapply(object$models, stats::BIC) |> print(quote = FALSE)
}
