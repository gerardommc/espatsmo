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
#' observation localitites if bias.data is a set of point localities.
#' @param positive Logical, specciffying whether density of sampling localities should be positive.
#' @param kernel A character string, specifying the type of kernel for the sampling localities.
#' The defaiult value is "gaussian", but all the types of kernels in spatstat.explore::density.ppm are supported.
#' @param sigma A positive real-valued numeric, indicating the standard deviation in distance units to generate the kernel density 
#' @param varcov A two by two matrix for the covariance between x and y sampling localities if these are correlated
#' @param weights A numeric vector with the length of sampling localities with the weights for each observation.
#' @param edge Logical, whether to perform edge correction
#' @param zo.norm Logical, whether to coerce observation effort values to 0-1 scale
#' @return A quadscheme but with modified, spatially variable area weights.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |> terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |> read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |>  terra::rast()
#' 
#' iml <- imFromStack(r$bio1)
#' 
#' p.pp <- spatstat.geom::ppp(x = p$x, y = p$y, window = as.owin(iml))
#' 
#' Q <- spatstat.geom::pixelquad(p.pp)
#' 
#' QA <- replaceQAreas(Q = Q,
#'                     bias.data = bias,
#'                     im = iml)
#' @export

replaceQAreas <- function(Q = NULL, 
                          bias.data = NULL,
                          im = NULL, 
                          positive = TRUE, 
                          kernel = "gaussian",
                          sigma = NULL, 
                          varcov = NULL, 
                          weights = NULL, 
                          edge = TRUE,
                          zo.norm = FALSE){

  
  if(inherits(bias.data, "data.frame")){
    window <- Q$dummy$window
    
    if(inherits(im, "SpatRaster")){r <- im}
    if(inherits(im, "im")){r <- terra::rast(im)}
    if(inherits(im, "imList")){r <- terra::rast(im[[1]])}
    
    r.counts <- terra::rasterize(bias.data, r, fun = "count") |> terra::as.data.frame(xy = TRUE)
    
    sample.ppp <- spatstat.geom::ppp(x = r.counts$x, y = r.counts$y, marks = r.counts$count, window = window)
    
    df.weights <- terra::as.data.frame(r, xy = TRUE)
    df.weights[, 3] <- sum(r.counts$count)/nrow(df.weights)
    
    d.E <- terra::rast(df.weights) |> imFromStack()
    
    d.obs <- spatstat.explore::density.ppp(sample.ppp, positive = positive, kernel = kernel,
                                           sigma = sigma, varcov = varcov, 
                                           weights = marks, edge = edge)
    
    SampleRatios <- d.obs/d.E
    
    begin <- Q$w |> length() - length(SampleRatios[]) +1
    end <- Q$w |> length()
    
    NewAreas <- Q$w[begin:end] * SampleRatios[]
    Q$w[begin:end] <- NewAreas
    
    return(Q)
  }
  
  if(inherits(bias.data, "SpatRaster")){
      im.r <- terra::rast(im)
      terra::crs(im.r) <- terra::crs(bias.data)
    
      if(zo.norm){
        bias.data <- terra::resample(bias.data, im.r) |> ZeroOneNorm()
      } else {
        bias.data <- terra::resample(bias.data, im.r) 
      }
      
      sum.bias <- terra::global(bias.data, "sum", na.rm = TRUE)
      bias.data <- bias.data/sum.bias$sum 
      
      bias.data <- bias.data * Q$data$n
      
      ones <- terra::classify(bias.data, rcl = matrix(c(-Inf, Inf, 1), ncol = 3))
      
      n.pix <- terra::global(ones, "sum", na.rm = TRUE)
      
      SampleRatios <- imFromStack(bias.data/(ones/n.pix$sum))
    
      begin <- Q$w |> length() - length(SampleRatios[]) +1
      end <- Q$w |> length()
      
      NewAreas <- Q$w[begin:end] * SampleRatios[]
      Q$w[begin:end] <- NewAreas
      
      return(Q)
  }
    
  if(inherits(bias.data, "im")){
      im.r <- terra::rast(im)
      bias.data <- terra::rast(bias.data)
    
      if(zo.norm){
        bias.data <- terra::resample(bias.data, im.r) |> ZeroOneNorm()
      } else {
        bias.data <- terra::resample(bias.data, im.r) 
      }
      
      sum.bias <- terra::global(bias.data, "sum", na.rm = TRUE)
      bias.data <- bias.data/sum.bias$sum 
      
      bias.data <- bias.data * Q$data$n
      
      ones <- terra::classify(bias.data, rcl = matrix(c(-Inf, Inf, 1), ncol = 3))
      
      n.pix <- terra::global(ones, "sum", na.rm = TRUE)
      
      SampleRatios <- imFromStack(bias.data/(ones/n.pix$sum))

      begin <- Q$w |> length() - length(SampleRatios[]) +1
      end <- Q$w |> length()
      
      NewAreas <- Q$w[begin:end] * SampleRatios[]
      Q$w[begin:end] <- NewAreas
      
      return(Q)
  }
}
