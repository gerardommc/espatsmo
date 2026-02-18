#' @title Mask unsampled regions of a set of covariates given a spatially variable sampling effort and presence points
#' @description
#' The function takes at least three arguments: (1) the set of covariates as SpatRaster or list of images
#' where pixels will be eliminated; (2) the bias data, which may be an externally generated raster layer by 
#' any means considered suitable, or a data frame with coordinates in the same CRS as covariates, which will be
#' used to interally generate a raster layer using spatstat Kernel density estimation density.ppp; and (3) the 
#' presence points to be modelled, to ensure that these pixels are included in the filtered covariates.
#' The remaining arguments are used to configure the density.ppp call.
#' @param covariates A set of covariates as SpatRaster or imList object. 
#' @param bias.data The data representing the variability of observation effort, which may be a two column data.frame of sampling localities a SpatRaster or imList.
#' @param points A two-column data.frame f x, y presence localitites for the species being modelled.
#' @param positive Logical. If bias data are a data.frame, the function will use spatstat.geom::density.ppp, and configures whether the  Kernel density will be forced to be positive.
#' @param kernel Character, to specify the type of smoothing kernel to be used by density.ppp.
#' @param sigma Numeric double, to specify the bandwidth in distance units of the kernel smoothing with density.ppp.
#' @param varcov A 2x2 covariance matrix used to configure the kernel smoothing with density.ppp.
#' @param weights Logical, whether weights willbe used in the kernel smoothing.
#' @param edge Logical, whether the kernel smoothing will have edge correction.
#' @param p.keep Numeric double in the interval 0-1, to specify the proportion of pixels that will be retained in the filtered covariates.
#' @param as.imList Logical, whether the returned object will be an imList or SpatRaster.
#' @return A SpatRaster or imList with the same covariates as the entered set.

maskBias <- function(covariates = NULL, 
                     bias.data = NULL,
                     points = NULL,
                     positive = TRUE, kernel = "gaussian",
                     sigma = NULL, varcov = NULL, 
                     weights = NULL, edge = TRUE,
                     p.keep = 0.5,
                     as.imList = F){
  require(terra)
  require(spatstat)

  if(is.null(covariates) | is.null(bias.data) | is.null(points)){
    stop("Cannot filter null data, please provide a sample SpatRaster, bias data and presence points")
  }
  
  if(class(covariates) == "SpatRaster"){
    imList <- imFromStack(covariates)
    w <- spatstat.geom::as.owin(imList[[1]])
  }
  
  if(class(covariates) == "imList"){
    imList <- covariates
    w <- spatstat.geom::as.owin(imList[[1]])
  }

  if(!class(bias.data) == "SpatRaster"){
    if(class(bias.data) == "data.frame"){
      pp <- spatstat.geom::ppp(x = bias.data$x, y = bias.data$y, window = w)
      Q <- spatstat.geom::pixelquad(pp)
    } 
    
   if(class(bias.data) == "ppp"){
      pp <- bias.data
      Q <- spatstat.geom::pixelquad(pp)
    }

    if(class(bias.data) == "quad"){
      pp <- spatstat.geom::ppp(x = bias.data$data$x, y = bias.data$data$y, window = w)
      Q <- bias.data
    }

    bias.layer <- spatstat.explore::density.ppp(pp, 
                              positive = positive, 
                              kernel = kernel,
                              sigma = sigma, varcov = varcov, 
                              weights = weights, edge = edge)
    
    bias.layer.r <- terra::rast(bias.layer) |> ZeroOneNorm()
  }
  
  if(class(bias.data) == "imList"){
    bias.layer.r <- terra::rast(bias.data)
  }
  
  if(class(bias.data) == "SpatRaster"){
    bias.layer.r <- bias.data
  }
  
  bias.df <- as.data.frame(bias.layer.r, xy = T)
  
  samp <- sample(1:nrow(bias.df), 
                 size = p.keep*nrow(bias.df), 
                 replace = F, prob = bias.df[, 3]) |> sort()
  
  if(class(points) == "ppp"){
    points <- data.frame(x = points$x, y = points$y)
  }

  if(class(points) == "quad"){
    points <- data.frame(x = points$data$x, y = points$data$y)
  }
  
  pres.r <- terra::rasterize(as.matrix(points), covariates[[1]])
  pres.r <-  terra::classify(pres.r, rcl = matrix(c(1, Inf, 1, 0, 0, NA), byrow = T, ncol = 3))
  
  samp.r <- terra::rasterize(bias.df[-samp, c("x", "y")], covariates[[1]])
  samp.r <- terra::classify(samp.r, rcl = matrix(c(1, Inf, 1, 0, 0, NA), byrow = T, ncol = 3))
  
  filt.r <- terra::merge(pres.r, samp.r)

  covariates.masked <- terra::mask(covariates, filt.r, inverse = T)

  if(as.imList){
    covariates.masked.im <- imFromStack(covariates.masked)
    return(covariates.masked.im)
  } else {
    return(covariates.masked)
  }
}

