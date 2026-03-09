#' @title Modify the area weights from a quadscheme to correct for observation bias
#' @description
#' If observation bias is believed to affect a point process, and sampling effort variability is understood,
#' this function provides a way of correcting for sampling bias, by modiffying the size of the area weights
#' in relation to an assumed uniform sampling scenario.
#' @param Q A quad object 
#' @param bias.data A two o three-column data frame containing the points of the sampling localitites, 
#' where the first two should be the x and y coordinates and the third the weights. Alternatively, a SpatRaster or
#' im object representing the spatial variability of observation effor.
#' @param im an im object with the same resolution of the quadrature scheme in the Q argument. 
#' @param nsim integer, the number of simulations used to derive the expected point intensity of 
#' observation localitites if bias.data is a set of point localities.
#' @param positive Logical, specciffying whether density of sampling localities should be positive.
#' @param kernel A character string, specifying the type of kernel for the sampling localities.
#' The defaiult value is "gaussian", but all the types of kernels in spatstat.explore::density.ppm are supported.
#' @param sigma A positive real-valued numeric, indicating the standard deviation in distance units to generate the kernel density 
#' @param varcov A two by two matrix for the covariance between x and y sampling localities if these are correlated
#' @param weights A numeric vector with the length of sampling localities with the weights for each observation.
#' @param edge Logical, whether to perform edge correction
#' @return A quadscheme but with modified, spatially variable area weights.

replaceQAreas <- function(Q = NULL, 
                          bias.data = NULL,
                          im = NULL, 
                          positive = TRUE, 
                          kernel = "gaussian",
                          sigma = NULL, 
                          varcov = NULL, 
                          weights = NULL, 
                          edge = TRUE){

  
  if(class(bias.data) == "data.frame"){
      window <- Q$dummy$window
    
    if(inherits(im, "SpatRaster")){r <- im}
    if(inherits(im, "im")){r <- terra::im(im)}
    if(inherits(im, "imList")){r <- terra::rast(im[[1]])}
    
    r.counts <- terra::rasterize(bias.data, r, fun = "count") |> terra::as.data.frame(xy = T)
    
    sample.ppp <- spatstat.geom::ppp(x = r.counts$x, y = r.counts$y, marks = r.counts$count, window = window)
    
    df.weights <- terra::as.data.frame(r, xy = T)
    df.weights[, 3] <- sum(r.counts$count)/nrow(df.weights)
    
    d.E <- terra::rast(df.weights) |> imFromStack()
    
    d.obs <- spatstat.explore::density.ppp(sample.ppp, positive = positive, kernel = kernel,
                                           sigma = sigma, varcov = varcov, 
                                           weights = marks, edge = edge)
    
    AreaWeights <- d.obs/d.E
    
    begin <- Q$w |> length() - length(AreaWeights[]) +1
    end <- Q$w |> length()
    
    NewAreas <- Q$w[begin:end] * AreaWeights[]
    Q$w[begin:end] <- NewAreas
    
    return(Q)}
  
  if(!class(bias.data) == "data.frame"){
    
    if(inherits(bias.data, "SpatRaster")){
      im.r <- terra::rast(im)
      crs(im.r) <- terra::crs(bias.data)
      bias.data <- terra::resample(bias.data, im.r)
      
      zo <- ZeroOneNorm(bias.data, as.imList = F)
      sum.zo <- terra::global(zo, "sum", na.rm = T)
      zo <- zo/sum.zo$sum
      
      SampleRatios <- imFromStack(zo/sum.zo$sum)
    }
    
    if(inherits(bias.data, "im")){
      im.r <- terra::rast(im)
      bias.data <- terra::rast(bias.data)
      bias.data <- terra::resample(bias.data, im.r)
      
      zo <- ZeroOneNorm(bias.data, as.imList = F)
      sum.zo <- terra::global(zo, "sum", na.rm = T)
      zo <- zo/sum.zo$sum
      
      SampleRatios <- imFromStack(zo/sum.zo$sum)
    }
    
    begin <- Q$w |> length() - length(SampleRatios[]) +1
    end <- Q$w |> length()
    
    TotalArea <- Q$w[begin:end] |> sum()
    
    AreaWeights <- SampleRatios * TotalArea
    
    NewAreas <- Q$w[begin:end] * AreaWeights[]
    Q$w[begin:end] <- NewAreas
    
    return(Q)
  }
}
