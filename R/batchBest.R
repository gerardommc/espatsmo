#' @title Extract the model with the lowest AIC from a ppmBatch object
#' @description The function requires a batch of point process models
#' fitted with ppmBatchFit, and whether the returned object is a plain
#' spatstat ppm object or a ppmSingle objectm via the as.ppmSingle argument
#' Independently, of the latter being set as TRUE, the ppm object
#' is contained within the model slot of the returned ppmSingle object
#' @param batch A ppmBatch object from which to extract the "best" model, on the basis of the AIC.
#' @param as.ppmSingle Logical, whether the returned object is a ppSingle, or spatstat.model::ppm to be handled with spatstat.
#' @return A ppmSingle or ppm object.
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
#' models <- ppmBatchFit(points = ,
#'                       covariates = r,
#'                       formulas = forms,
#'                       bias.data = bias,
#'                       bias.correction = "weights",
#'                       parallel = F,
#'                       top.models = 3)
#' 
#' best <- batchBest(batch = models,
#'                   as.ppmSingle = F)
#' 
#' spatstat.model::summary(best)
#' }
#' @export

batchBest <- function(batch = NULL, as.ppmSingle = T){
  if(is.null(batch)){
    stop("Please provide a batch of models fitted with ppmBatchFit")
  }
  
  m <- batch$models[[1]]
  
  if(as.ppmSingle){
    ret.list <- list(model = m,
                     call = list(
                       bias.data = batch$call$bias.data,
                       bias.correction = batch$call$bias.correction,
                       weight.bias.conf = batch$call$weight.bias.conf))
    
    class(ret.list) <- "ppmSingle"
    
    return(ret.list)} else {
      
      return(m)
    }
}