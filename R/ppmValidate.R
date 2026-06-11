#' @title Generate validation statistics
#' @description
#' Quantify the predictive capacity of a point process model by block partitioning 
#' @param model A model object of class ppm, ppmSingle or ppmBatch
#' @return A summary of the validation statistics generated
#' @export
#' @method ppmValidate

ppmValidate <- function(model, ...){
    UseMethod(generic = "ppmValidate", model)
}
