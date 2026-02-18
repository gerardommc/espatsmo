#' @title Compute and print the AIC of single or batch point process models
#' @description
#' Print the Akaike information criterion of one or multiple point proces models fitted with spatstat's ppm, but via the inteface of the presente package.@
#' @param model A ppmBatch or ppmSingle, fitted with the corresponding functions
#' @return The function returns the Akaike Information Criterion, one for a ppmSigle object, or the values of all the returned mdels with ppmBatchFit

AIC <- function(model){
  UseMethod("AIC", model)
}

AIC.ppmBatch <- function(model){
  sapply(model$models, spatstat.model::AIC.ppm) |> print(quote = F)
}

AIC.ppmSingle <- function(model){
  spatstat.model::AIC.ppm(model)
}


