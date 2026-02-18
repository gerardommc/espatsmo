#' @title Automatically generate polynomial formulas using a covariate compatibility matrix and the exponents for each variable
#' @description This function takes two arguments, a data.fame with two columns with the specific names "Variable" and "Power",
#' and a covariate compatibility matrix, usually generated with the function findCompatibles. The variable names provided in the respDF argument
#' MUST match those in the covariate compatibility matrix. Both arguments are used to generate a "character" type of vector, where each element
#' is a string convertible to a formula object. the resulting formulas contain a polynomial of the order specified in respDF funciton of 
#' each subset of variables in compatMat
#' @param respDF A data.frame containing two columns with names "Variable" and "Power", where foreach variable name in the compatibility matrix returned by findCompatibles, there should be a proposed positive integer exponent.
#' @param compatMat A covariate compatibility matrix generated with findCompatibles.
#' @return A character vector containing the model formulas for each row of variable cobinations in compatMat, and here each term in the formula has an exponent 1:Power.


getPolyFormulas <- function(respDF = NULL, 
                            compatMat = NULL){
                              
  require(foreach)

  
  if(is.null(respDF)){
    stop("Cannot produce formulas, please provide a valid description of variable responses")
  }

  if(is.null(compatMat)){
    stop("Cannot produce formulas, please provide a valid compatibility matrix")
  }
    
  formulas <-  foreach::foreach(i = 1:nrow(compatMat), .combine = c) %do% {
    v <- compatMat[i, ]
    rd <- respDF[which(respDF$Variable %in% v), ]
    
    exponents <- foreach::foreach(j = 1:nrow(rd)) %do% {
      p <- 1:rd$Power[j]
    }
    
    var.exponents <- foreach::foreach(k = seq_along(exponents), .combine = c) %do% {
      ifelse(exponents[[k]] == 1, v[k], paste0("I(", v[k], "^", exponents[[k]], ")"))
    }
    
    f1 <- paste0("~",var.exponents[1])
    for(ii in 2:length(var.exponents)){
      f1 <- paste(f1, var.exponents[ii], sep = " + ")
    }
    return(f1)
  }
  
  return(formulas) 
}


