#' @title Fit one single point process model
#' @description 
#' @param points A two-clumn data.frame, ppp or quad object containing the presence localities 
#' @param covariates A SpatRaster or imList containing the environmental covariates 
#' @param formula A character vector containing formula to be fitted. 
#' @param bias.data  A two-column data frame, containing the coordinates of the sampling localities, 
#' SpatRaster or imList object containing the spatial variability of the observation effort
#' @param bias.correction character, with values 'weights' or 'background', 
#' specifying the method used to control the effect of sampling bias
#' @param weight.bias.conf list, containing the following elements: 
#' 1) nsim, 2) positive, 3) kernel, 4) sigma, 5) varcov, 6) weights, 7) edge and 8) p.keep. 
#' Where `p.keep` is relevant bias.correction if bias.correction = 'background' and specifies the proportion of
#' pixels to be retained after thinning the covariates. For the remaining elements of the list, please consult the help files of `spatstat.explore::density.ppm`.
#' @param ppm.conf list, containing the elemens to configure `spatstat.model:ppm`: 1) na.action, 2) correction, 3) use.gam=FALSE,
#' 4) method, 5) forcefit, 6) improve.type, 7) improve.args, 8) prior.mean, 9) prior.var. Plase consult the help files for
#' `spatstat.model:ppm`
#' @return A list of class `ppmSingle` containing the specified number of models to bee returned and the configuration of the call.



ppmSingleFit <- function(points= NULL, 
                covariates = NULL, 
                formula = NULL, 
                bias.data = NULL, #Data frame with sampling localities or raster layer
                bias.correction = NULL,
                weight.bias.conf = list(nsim = 39,
                                        positive = TRUE, kernel = "gaussian",
                                        sigma = NULL, varcov = NULL,
                                        weights = NULL, edge = TRUE ,
                                        p.keep = 0.5),
                ppm.conf = list(na.action = na.exclude,
                                correction="border",
                                use.gam=FALSE,
                                method="logi",
                                forcefit=FALSE,
                                improve.type = "none",
                                improve.args=list(),
                                prior.mean = NULL,
                                prior.var = NULL)){
  require(spatstat)
  require(terra)
  
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
      
      if(class(points) == "data.frame"){
        pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
        Q <- spatstat.geom::pixelquad(pp)
      } 
      
      if(class(points) == "ppp") {
        Q <- spatstat.geom::pixelquad(points)
      }

      if(class(points) == "quad"){
        Q <- points
      }
    }
    
    if(bias.correction == "weights"){
      if(class(bias.data) == "data.frame"){

        imList <- imFromStack(covariates)
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = imList[[1]], 
                            nsim =  weight.bias.conf$nsim,
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma, 
                            weights = weight.bias.conf$weights, 
                            edge = weight.bias.conf$edge)
        Q <- Qa
      }
      
      if(class(bias.data) == "SpatRaster"){
        
        imList <- imFromStack(covariates)
        
        bias.data <- bias.data |> ZeroOneNorm()
        
        sum.bias <- terra::global(bias.data, sum, na.rm = T, ID = F)[, 1]
        
        bias.data <- bias.data/sum.bias * nrow(points)
        
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = imList[[1]], 
                            nsim =  weight.bias.conf$nsim,
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

    if(class(covariates) == "SpatRaster"){
      imList <- imFromStack(covariates)
      w <- spatstat::as.owin(imList[[1]])
    } 
    
    if(class(covariates) == "imList"){
      imList <- covariates
      w <- spatstat::as.owin(imList[[1]])
    }
    
    if(class(points) == "data.frame"){
      pp <- spatstat::ppp(x = points$x, y = points$y, window = w)
      Q <- spatstat::pixelquad(pp)
    }
    
    if(class(points) == "ppp"){
      Q <- spatstat::pixelquad(points)
    }  

    if(class(points) == "quad"){
      Q <- points
    } 
  }
  
  m <- spatstat.model::ppm(Q, 
                           trend = as.formula(formula[i]), 
                           covariates = imList,
                           na.action = ppm.conf$na.action,
                           correction = ppm.conf$correction,
                           use.gam = ppm.conf$use.gam,
                           method = ppm.conf$method,
                           forcefit = ppm.conf$forcefit,
                           improve.type = ppm.conf$improve.type,
                           improve.args=ppm.conf$improve.args,
                           prior.mean = ppm.conf$prior.mean,
                           prior.var = ppm.conf$prior.var)
  
  ret.list <- list(model = m,
                   call = list(bias.data = bias.data, #Data frame with sampling localities or raster layer
                               bias.correction = bias.correction,
                               weight.bias.conf = weight.bias.conf))
  
  class(ret.list) <- "ppmSingle"
  
  return(ret.list)
}