#' @title Automatically generate GAM formulas using a covariate compatibility matrix and the desired type of smoothers for each variable
#' @description This function takes two arguments, a data.fame with two columns with the specific names "Variable" and "Smoother",
#' and a covariate compatibility matrix, usually generated with the function findCompatibles. The variable names provided in the respDF argument
#' MUST match those in the covariate compatibility matrix. Both arguments are used to generate a "character" type of vector, where each element
#' is a string convertible to a formula object passed on to a gam or ppm with use.gam = T. 
#' @param respDF A data.frame containing two columns with names "Variable" and "Smoother", where foreach variable name in the compatibility matrix returned by findCompatibles, a smoother type used by gam, is specified.If no smoother is desired for a variable, this value should be set to 1.
#' @param compatMat A covariate compatibility matrix generated with findCompatibles.
#' @return A character vector containing the model formulas for each row of variable cobinations in compatMat.
#' @examples
#' \dontrun{
#' r <- terra::rast("inst/extdata/ChelsaBio.tif")
#' resp <- read.csv("inst/extdata/Smoothers.csv")
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getGAMFormulas(respDF = resp, 
#'                          compatMat = compat)
#' }
#' @export

getGAMFormulas <- function(respDF, compatMat){
  
  formulas <-  foreach::foreach(i = 1:nrow(compatMat), .combine = c) %do% {
    v <- compatMat[i, ]
    rd <- respDF[which(respDF$Variable %in% v), ]
    
    smoothType <- rd$Smoother
    
    responseTypes <- foreach::foreach(k = seq_along(smoothType), .combine = c) %do% {
      ifelse(smoothType[k] == 1, v[k], paste0(smoothType[k],"(", v[k], ")"))
    }
    
    f1 <- paste0("~ ", paste("", responseTypes, collapse =  " + "))

    return(f1)
  }

  return(formulas) 
}


