#' @title Estract the latent field 
#' @description
#' Thin function takes a ppmLGCP object and extracts the latent field as a SpatRaster
#' @param model A model fitted with ppmLGCP
#' @return A terra package SpatRaster object.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' s <-  system.file("extdata", "ClusterRandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' p <- system.file("extdata", "pointsCluster.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#'  
#' model <- ppmLGCP(points= pr, 
#'                  covariates = r, 
#'                  formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                  dist.ar = FALSE,
#'                  weight.units = "km",
#'                  coordinates = "m")
#' 
#' latent <- latentField(model)
#' @export

latentField <- function(model = NULL){
  
  if(!inherits(model, "ppmLGCP")){
      stop("Cannot extract the latent field from the provided object, please provide a ppmLGCP model.")
    }
    
  latent <- model$predictions$mean/model$base.predictions$preds.mean

  return(latent)
}