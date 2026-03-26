#' @title Non parametric smoothed plot of the response of a point process to the values of the supplied 
#' covariates.
#' @description
#' This function takes at least two arguments: (1) a point process in the form of a two-column data.frame or ppp object, 
#' and (2) a set of covariates of class SpatRaster or list of images. The remaining arguments are used to configure the 
#' bias correction and therefore plot the responses to model according to the bias scenario. The bias correction methods
#' are either "weights" or "background". The "weights" method, provided in the argument bias.correction, changes the 
#' size of the areas in relation to a uniform sampling scenario, whereas "background" uses maskBias to eliminate unsampled pixels
#' from covariates regions. the various arguments to configure bias correction procedures are provided in the weight.bias.conf argument
#' which are passed to either replaceQAreas or maskBias. The funciton's output is a pdf file with the default name "ResponsePlot.pdf"
#' saved in the working directory, but can be changed using save.plot and plot.pars
#' @param points A two-column data.frame f x, y presence localitites for the species being modelled.
#' @param covariates A set of covariates as SpatRaster or imList object.
#' @param bias.data The data representing the variability of observation effort, which may be a two column data.frame of sampling localities a SpatRaster or imList.
#' @param bias.correction Character string, indicating whether to correct bias with background filtering of by altering the area weight.
#' @param weight.bias.conf A list arguments to pass to density.ppp to configure weight bias correction methos.
#' @param save.plot Logical to indicate if a plot of the responses is to be saved.
#' @param plot.pars A list to configure plot size and file name if save.plot is TRUE.
#' @param p.keep Numeric double in the interval 0-1, to specify the proportion of pixels that will be retained in the filtered covariates.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' plotResponses(points = p,
#'               covariates = r,
#'               bias.data = bias,
#'               bias.correction = "weights",
#'               save.plot = TRUE,
#'               plot.pars = list(name = "inst/extdata/ResponsePlot.pdf", 
#'                                width = 5, 
#'                                height = 5)) 
#' @export

plotResponses <- function(points = NULL,
                          covariates = NULL,
                          bias.data = NULL, #Data frame with sampling localities or raster layer
                          bias.correction = NULL,
                          weight.bias.conf = list(positive = TRUE, kernel = "gaussian",
                                                  sigma = NULL, varcov = NULL, 
                                                  weights = NULL, edge = TRUE,
                                                  p.keep = 0.5),
                          save.plot = FALSE,
                          plot.pars = list(name = "ResponsePlot.pdf", width = 5, height = 5)){

  
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
      Q <- spatstat.geom::pixelquad(pp)
    }
  }

  if(!is.null(bias.correction)){
    if(bias.correction == "weights"){

      if(inherits(covariates, "imList")){
        imList <- covariates
        w <- spatstat.geom::as.owin(imList[[1]])
      }

      if(inherits(covariates, "SpatRaster")){
        imList <- imFromStack(covariates)
        w <- spatstat.geom::as.owin(imList[[1]])
      }

      if(inherits(points, "data.frame")){
        pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
        Q <- spatstat.geom::pixelquad(pp)
      }
  
      if(inherits(points, "ppp")){
        Q <- spatstat.geom::pixelquad(pp)
      }

      Q <- replaceQAreas(Q = Q,
                         bias.data = bias.data,
                         im = imList[[1]],
                         positive = weight.bias.conf$positive,
                         kernel = weight.bias.conf$kernel,
                         sigma = weight.bias.conf$sigma,
                         varcov = weight.bias.conf$varcov,
                         weights = weight.bias.conf$weights,
                         edge = weight.bias.conf$edge)
    }
      
     if(bias.correction == "background"){
       
       if(is.null(weight.bias.conf$p.keep)){
         stop("Please provide a value for p.keep between 0 and 1 (proportion of pixels to be skimmed)")
       } else {
         if(weight.bias.conf$p.keep > 1 | weight.bias.conf$p.keep < 0){
                    stop("Please provide a value for p.keep between 0 and 1 (proportion of pixels to be skimmed)")
         }
       }
    
       imList <- maskBias(covariates = covariates,
                            bias.data = bias.data,
                            points = points,
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma, 
                            varcov = weight.bias.conf$varcov, 
                            weights = weight.bias.conf$weights,
                            edge = weight.bias.conf$edge,
                            p.keep = weight.bias.conf$p.keep,
                            as.imList = TRUE)
        
      w <- spatstat.geom::as.owin(imList[[1]])
      
      if(inherits(points, "data.frame")){
        pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
        Q <- spatstat.geom::pixelquad(pp)
      }
      
      if(inherits(points, "ppp")){
        pp <- points
        Q <- spatstat.geom::pixelquad(pp)
      }
        
      if(inherits(points, "quad")){
        pp <- spatstat.geom::ppp(x = points$data$x, y = points$data$y, window = w)
        Q <- points
      }
     }
  }
  
  if(save.plot){
    grDevices::pdf(plot.pars$name, width = plot.pars$width, height = plot.pars$height)
    for(i in seq_along(imList)){
      spatstat.explore::rhohat(Q, covariate = imList[[i]]) |> plot(main = names(imList)[i])
    }
    grDevices::dev.off()
  }else{
    for(i in seq_along(imList)){
      spatstat.explore::rhohat(Q, covariate = imList[[i]]) |> plot(main = names(imList)[i])
    }
  }
}