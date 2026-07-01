#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmBatch, the summary is for the id models contained in the list of fitted models
#' @param object A batch of point process models fitted with ppmBatchFit
#' @param id The id of the model(s) to print the summary.
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
#' summary(models, id = 2)
#' @export
#' @method summary ppmBatch


summary.ppmBatch <- function(object = NULL, id = NULL){
  spatstat.model::summary.ppm(object = object$models[[id]])
}


