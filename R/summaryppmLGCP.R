#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmLGCP, the summary is for either fixed or random effects
#' @param object A model fitted with ppmLGCP
#' @param effects A character string specifying whether to print `fixed` or `random` effects or both
#' @examples
#' \dontrun{
#' r <- terra::rast("inst/extdata/ChelsaBio.tif")
#' 
#' p <- read.csv("inst/extdata/points.csv")
#' 
#' bias <- terra::rast("inst/extdata/Target-group.tif")
#' 
#' model <- ppmLGCP(points= p, 
#'                  covariates = r, 
#'                  formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                  bias.data = bias,
#'                  bias.correction = "weights",
#'                  as.ppmSingle = F)
#' 
#' summary(model)
#' }
#' @export
#' @method summary ppmLGCP


summary.ppmLGCP <- function(object = NULL, effects = NULL){
  if(is.null(effects)){
    pl <- list(`Fixed effects` = object$model$summary.fixed,
               `Random effects` = object$model$summary.hyperpar)
    print(pl, quote = F)
  }
  
  if(!is.null(effects)){
      if(effects == "fixed"){
          print(object$model$summary.fixed, quote = F)
      }
      
      if(effects == "random"){
          print(object$model$summary.hyperpar, quote = F)
      }
    }
}



