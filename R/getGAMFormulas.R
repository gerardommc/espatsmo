#' @title Automatically generate GAM formulas using a covariate compatibility matrix and the desired type of smoothers for each variable
#' @description This function takes two arguments, a data.fame with two columns with the specific names "Variable" and "Smoother",
#' and a covariate compatibility matrix, usually generated with the function findCompatibles. The variable names provided in the respDF argument
#' MUST match those in the covariate compatibility matrix. Both arguments are used to generate a "character" type of vector, where each element
#' is a string convertible to a formula object passed on to a gam or ppm with use.gam = T. 
#' @param respDF A data.frame containing two columns with names "Variable" and "Smoother", where foreach variable name in the compatibility matrix returned by findCompatibles, a smoother type used by gam, is specified.If no smoother is desired for a variable, this value should be set to 1.
#' @param compatMat A covariate compatibility matrix generated with findCompatibles.
#' @return A character vector containing the model formulas for each row of variable cobinations in compatMat.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast()
#' 
#' resp <- system.file("extdata", "Smoothers.csv", package = "espatsmo") |>  read.csv()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getGAMFormulas(respDF = resp, 
#'                          compatMat = compat)
#' @export

getGAMFormulas <- function(respDF, compatMat){

  `%do%` <- foreach::`%do%`
  
  formulas <- foreach::foreach(i = 1:nrow(compatMat), .combine = c) %do% {
    v <- compatMat[i, ]
    rd <- respDF[which(respDF$Variable %in% v), ]
    smoothType <- rd$Smoother
    
    conf <- list(k = rd$k, bs = rd$bs, m = rd$m,
                 d = rd$d, by = rd$by, fx = rd$fx,
                 np = rd$np, xt = rd$xt,id = rd$id,
                 sp = rd$sp, mc = rd$mc,pc = rd$pd)
    
    if(is.null(conf$k)){ conf$k <- "NA" }  
    if(is.null(conf$bs)){ conf$bs <- "cr" }
    if(is.null(conf$m)){ conf$m <- "NA" } 
    if(is.null(conf$d)){ conf$d <- "NA" }  
    if(is.null(conf$by)){ conf$by <- "NA" }
    if(is.null(conf$fx)){ conf$fx <- "FALSE" }
    if(is.null(conf$np)){ conf$np <- "TRUE" }
    if(is.null(conf$xt)){ conf$xt <- "NULL" }
    if(is.null(conf$id)){ conf$id <- "NULL" }
    if(is.null(conf$sp)){ conf$sp <- "NULL" }
    if(is.null(conf$mc)){ conf$mc <- "NULL" }
    if(is.null(conf$pc)){ conf$pc <- "NULL" }
    
    responseTypes <- foreach::foreach(ii = seq_along(smoothType), 
                                      .combine = c) %do% {
                                        
                                        conf.ii <- paste0(",k = ", conf$k, ", bs = ", conf$bs, ", m = ", conf$m,
                                                          ", d = ", conf$d, ", by = ", conf$by, ", fx = ", conf$fx,
                                                          ", np = ", conf$np, ", xt = ", conf$xt, ", id = ", conf$id,
                                                          ", sp = ", conf$sp, ", mc = ", conf$mc, ", pc = ", conf$pc)
                                        
                                        f.ii <- ifelse(smoothType[ii] == 1, v[ii], 
                                                       paste0(smoothType[ii], "(", v[ii], conf.ii, ")"))
                                        
                                        return(f.ii)
                                      }
    
    f1 <- paste0("~ ", paste("", responseTypes, collapse = " + "))
    return(f1)
  }

  class(formulas) <- "gamforms"

  return(formulas) 
}


