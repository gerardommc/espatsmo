#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmBatch, the summary is for the id models contained in the list of fitted models
#' @examples
#' \dontrun{
#' r <- terra::rast("inst/extdata/ChelsaBio.tif")
#' 
#' p <- read.csv("inst/extdata/points.csv")
#' 
#' bias <- terra::rast("inst/extdata/Target-group.tif")
#' 
#' resp <- read.csv("inst/extdata/Exponents.csv")
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getPolyFormulas(respDF = resp, 
#'                          compatMat = compat)
#' 
#' models <- ppmBatchFit(points = ,
#'                       covariates = r,
#'                       formulas = forms,
#'                       bias.data = bias,
#'                       bias.correction = "weights",
#'                       parallel = F,
#'                       top.models = 3)
#' 
#' summary(models, id = 2)
#' }
#' @export
#' @method summary ppmBatch


summary.ppmBatch <- function(model = NULL, id = NULL){
  spatstat.model::summary.ppm(object = model$models[[id]])
}


