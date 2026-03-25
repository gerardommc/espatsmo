#' @title Fit a point process model
#' @description The espatsmo interface to fit one point process model with spatstat but using generic raster and data.frame 
#' objects as inputs.
#' @param points A two-clumn data.frame, ppp or quad object containing the presence localities 
#' @param covariates A SpatRaster or imList containing the environmental covariates 
#' @param formula A character vector containing formula to be fitted. 
#' @param bias.data  A two-column data frame, containing the coordinates of the sampling localities, 
#' SpatRaster or imList object containing the spatial variability of the observation effort
#' @param bias.correction character, with values 'weights' or 'background', 
#' specifying the method used to control the effect of sampling bias
#' @param weight.bias.conf list, containing the following elements: 1) positive, 2) kernel, 3) sigma, 4) varcov, 5) weights, 6) edge and 7) p.keep. 
#' Where `p.keep` is relevant bias.correction if bias.correction = 'background' and specifies the proportion of
#' pixels to be retained after thinning the covariates. For the remaining elements of the list, please consult the help files of `spatstat.explore::density.ppm`.
#' @param ppm.conf list, containing the elemens to configure `spatstat.model:ppm`: 1) correction, 2) use.gam=FALSE,
#' 3) method, 4) forcefit, 5) improve.type, 6) improve.args, 7) prior.mean, 8) prior.var. Plase consult the help files for
#' `spatstat.model:ppm`
#' @param as.ppmSingle logical, whether to return a ppmSingle class object or a spatstat ppm model.
#' @return A list of class `ppmSingle` containing the specified number of models to bee returned and the configuration of the call.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' model <- ppmSingleFit(points= p, 
#'                      covariates = r, 
#'                      formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                      bias.data = bias,
#'                      bias.correction = "weights",
#'                      as.ppmSingle = F)
#' @export


ppmSingleFit <- function(points= NULL, 
                         covariates = NULL, 
                         formula = NULL, 
                         bias.data = NULL, #Data frame with sampling localities or raster layer
                         bias.correction = NULL,
                         as.ppmSingle = T,
                         weight.bias.conf = list(positive = TRUE, kernel = "gaussian",
                                        sigma = NULL, varcov = NULL,
                                        weights = NULL, edge = TRUE ,
                                        p.keep = 0.5),
                ppm.conf = list(correction="border",
                                use.gam=FALSE,
                                method="logi",
                                forcefit=FALSE,
                                improve.type = "none",
                                improve.args=list(),
                                prior.mean = NULL,
                                prior.var = NULL)){

  
  source("Spatstat-functions/imFromStack.R")
  source("Spatstat-functions/replaceQAreas.R")

  if(names(points)[1] != "x" & names(points)[2] != "y"){
    stop("Please change column names of presence points to \"x\" and \"y\" for long and lat respectively")
  }
  
  if(is.null(points) | is.null(formula) | is.null(covariates)){
    stop("Please provide a valid set of points, covariates and model formula")
  }

  if(!is.null(bias.correction)){
    if(bias.correction == "background"){
      imList <- maskBias(covariates = covariates, 
                         bias.data = bias.data,
                         points = points,
                         positive = weight.bias.conf$positive, 
                         kernel = weight.bias.conf$kernel,
                         sigma = weight.bias.conf$sigma, 
                         varcov = weight.bias.conf$sigma, 
                         weights = weight.bias.conf$weights, 
                         edge = weight.bias.conf$edge,
                         p.keep = weight.bias.conf$p.keep,
                         as.imList = T)
      
      w <- spatstat.geom::as.owin(imList[[1]])
      
      if(inherits(points, "data.frame")){
        pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
        Q <- spatstat.geom::pixelquad(pp)
      } 
      
      if(inherits(points, "ppp")) {
        Q <- spatstat.geom::pixelquad(points)
      }

      if(inherits(points, "quad")){
        Q <- points
      }
    }
    
    if(bias.correction == "weights"){
      if(inherits(bias.data, "data.frame")){

        imList <- imFromStack(covariates)
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = imList[[1]],
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma, 
                            weights = weight.bias.conf$weights, 
                            edge = weight.bias.conf$edge)
        Q <- Qa
      }
      
      if(inherits(bias.data, "SpatRaster")){
        
        imList <- imFromStack(covariates)
        
        bias.data <- bias.data |> ZeroOneNorm()
        
        sum.bias <- terra::global(bias.data, sum, na.rm = T, ID = F)[, 1]
        
        bias.data <- bias.data/sum.bias * nrow(points)
        
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = imList[[1]], 
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma, 
                            weights = weight.bias.conf$weights, 
                            edge = weight.bias.conf$edge)
        Q <- Qa
        
        
      }
    }
  } 
  
  if(is.null(bias.correction)){

    if(inherits(covariates, "SpatRaster")){
      imList <- imFromStack(covariates)
      w <- spatstat.geom::as.owin(imList[[1]])
    } 
    
    if(inherits(covariates, "imList")){
      imList <- covariates
      w <- spatstat.geom::as.owin(imList[[1]])
    }
    
    if(inherits(points, "data.frame")){
      pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
      Q <- spatstat.geom::pixelquad(pp)
    }
    
    if(inherits(points, "ppp")){
      Q <- spatstat.geom::pixelquad(points)
    }  

    if(inherits(points, "quad")){
      Q <- points
    } 
  }
  
  m <- spatstat.model::ppm(Q, 
                           trend = stats::as.formula(formula[i]), 
                           covariates = imList,
                           correction = ppm.conf$correction,
                           use.gam = ppm.conf$use.gam,
                           method = ppm.conf$method,
                           forcefit = ppm.conf$forcefit,
                           improve.type = ppm.conf$improve.type,
                           improve.args=ppm.conf$improve.args,
                           prior.mean = ppm.conf$prior.mean,
                           prior.var = ppm.conf$prior.var)
  
  if(as.ppmSingle){
        ret.list <- list(model = m,
                   call = list(bias.data = bias.data, #Data frame with sampling localities or raster layer
                               bias.correction = bias.correction,
                               weight.bias.conf = weight.bias.conf))
        class(ret.list) <- "ppmSingle"
  
        return(ret.list)
  } else {
    return(m)
  }
}