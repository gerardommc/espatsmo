#' @title Identify all the subsets of covariate combinations meeting a pre-defined correlation threshold
#' @description This function takes three arguments, the covariates must be a SpatRaster or data.frame containing only continuous variables
#' a tolerated correlation threshold and the maximum number of variables to be contained in each subset
#' @param covariates The set of covariates from which to identify sets of non-collinear axes.
#' @param thres Numeric positive between 0 - 1, representing the pearson correlation coefficient threshold of tolerated collinearity between covariates in each set of combinations
#' @param max.comb Numeric integer, representing the number of covariates for the combinations.
#' @return A data.frame with max.comb columns. The number of rows depends on the number combinations possible given the correlatin threshold specified.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' @export

findCompatibles <- function(covariates = NULL, 
                            thres = 0.6, 
                            max.comb = 3){

  if(is.null(covariates)){
    stop("Cannot compute compatibilities without data, please provide either a multi-band SpatRaster, data.frame or imList containing only continuous variables")
  }

  if(!class(covariates) %in% c("SpatRaster", "data.frame", "imList")){
    stop(paste0("Cannot compute compatibilities with a ", class(covariates), " please provide either a multi-band SpatRaster, data.frame or imList containing only continuous variables"))
  }

  if(inherits(covariates, "SpatRaster") & dim(covariates)[3] < 3){
    stop("Cannot compute compatibilities with less than three variables, please provide at least three covariates")
  }

  if(inherits(covariates, "data.frame") & ncol(covariates) < 3){
    stop("Cannot compute compatibilities with less than three variables, please provide at least three covariates")
  }

  if(inherits(covariates, "imList") & length(covariates) < 3){
    stop("Cannot compute compatibilities with less than three variables, please provide at least three covariates")
  } 

  if(inherits(covariates, "SpatRaster")){
    rdf <- as.data.frame(covariates, xy = FALSE) |> stats::na.omit()
  }

  if(inherits(covariates, "data.frame")){
    rdf <- covariates 
  }

  if(inherits(covariates, "imList")){
    r <- as.SpatRast(covariates)
    rdf <- as.data.frame(r, xy = FALSE)
  }

  cormat <- stats::cor(rdf)

  lt <- lower.tri(cormat)

  cormatl <- cormat
  cormatl[!lt] <- NA
  cormatl[!(cormatl[]<thres & cormatl[]>-thres)] <- F
  cormatl[cormatl!=0] <- T
  cormatl[is.na(cormatl)] <- F

  n <- names(rdf)

  combs <- utils::combn(n, max.comb)

  ref.size <- diag(max.comb) |> lower.tri() |> sum()

  lns <- numeric(ncol(combs))

  for(i in 1:ncol(combs)){
    lns[i] <- sum(cormatl[combs[,i], combs[,i]])
  }

  ids <- which(lns == ref.size)

  df <- data.frame(t(combs[,ids]))

  names(df) <- paste0("Variable_", 1:nrow(combs))

  return(df)
}