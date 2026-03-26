#' @title Generate summaries from model objects
#' @description These functions require a range of inputs depending on the class of the supplied
#' objects. For objects of class ppmLGCP, the summary is for either fixed or random effects
#' @param object A model fitted with ppmLGCP
#' @param effects A character string specifying whether to print `fixed` or `random` effects or both
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#'
#' model <- ppmLGCP(points= p, 
#'                  covariates = r, 
#'                  formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                  bias.data = bias,
#'                  bias.correction = "weights")
#' 
#' summary(model)
#' @export
#' @method summary ppmLGCP


summary.ppmLGCP <- function(object = NULL, effects = NULL){
  if(is.null(effects)){
    pl <- list(`Fixed effects` = object$model$summary.fixed,
               `Random effects` = object$model$summary.hyperpar)
    print(pl, quote = FALSE)
  }
  
  if(!is.null(effects)){
      if(effects == "fixed"){
          print(object$model$summary.fixed, quote = FALSE)
      }
      
      if(effects == "random"){
          print(object$model$summary.hyperpar, quote = FALSE)
      }
    }
}



