#' @title Fit batches of point process models
#' @description Fit multiple point process models and select a specified number on the basis of the AIC
#' @param points A two-clumn data.frame, ppp or quad object containing the presence localities 
#' @param covariates A SpatRaster or imList containing the environmental covariates 
#' @param formulas A character vector containing the linear or polynomial formulas generated with getPolyFormulas, to be fitted 
#' @param bias.data  A two-column data frame, containing the coordinates of the sampling localities, 
#' SpatRaster or imList object containing the spatial variability of the observation effort
#' @param bias.correction character, with values 'weights' or 'background', 
#' specifying the method used to control the effect of sampling bias
#' @param weight.bias.conf list, containing the following elements: 1) positive, 2) kernel, 3) sigma, 4) varcov, 5) weights, 6) edge and 7) p.keep. 
#' Where `p.keep` is relevant bias.correction if bias.correction = 'background' and specifies the proportion of
#' pixels to be retained after thinning the covariates. For the remaining elements of the list, please consult the help files of `spatstat.explore::density.ppm`.
#' @param parallel logical, to specify whether models will be fitted in parallel or in series. 
#' If TRUE, the doParallel package should be installed.
#' @param cores numeric integer, used to specify the number threads if parallel = TRUE 
#' @param top.models numeric integer, used to specify the number of fitted models to return on the basis of the AIC.
#' @param ppm.conf list, containing the elemens to configure `spatstat.model:ppm`: 1) correction, 2) use.gam=FALSE,
#' 3) method, 4) forcefit, 5) improve.type, 6) improve.args, 7) prior.mean, 8) prior.var. Plase consult the help files for
#' `spatstat.model:ppm`
#' @return A list of class `ppmBatch` containing the specified number of models to bee returned and the configuration of the call.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#'
#' resp <-  system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getPolyFormulas(respDF = resp, 
#'                          compatMat = compat)
#' 
#' models <- ppmBatchFit(points = p,
#'                       covariates = r,
#'                       formulas = forms,
#'                       bias.data = bias,
#'                       bias.correction = "weights",
#'                       parallel = FALSE,
#'                       top.models = 3)
#' @export


ppmBatchFit <- function(points= NULL, 
                        covariates = NULL, 
                        formulas = NULL, 
                        bias.data = NULL, #Data frame with sampling localities or raster layer
                        bias.correction = NULL,
                        goodness.fit = "AIC",
                        weight.bias.conf = list(positive = TRUE, kernel = "gaussian",
                                                sigma = NULL, varcov = NULL, 
                                                weights = NULL, edge = TRUE,
                                                p.keep = 0.5, zo.norm = FALSE), #Only relevant for Background bias correction
                        parallel = TRUE, 
                        cores = 2, 
                        top.models = 10,
                        ppm.conf = list(correction="border",
                                        use.gam=TRUE,
                                        method="mpl",
                                        forcefit=TRUE,
                                        improve.type = "none",
                                        improve.args=list(),
                                        prior.mean = NULL,
                                        prior.var = NULL)){ 
  
    `%do%` <- foreach::`%do%`
    `%dopar%` <- foreach::`%dopar%`
    
    if(length(formulas) == 1){
      stop("ppmBatchFit requires more than 1 formula, please use ppmSingleFit")
    }
  
    if(inherits(formulas, "gamforms")){ppm.conf$use.gam <- TRUE}
  
    if(!inherits(points, c("data.frame", "ppp", "quad"))){
      stop("Please provide presence points as either a two-column data.frame, ppp, or quad object")
    }
  
    if(inherits(points, "data.frame")){
      if(names(points)[1] != "x" & names(points)[2] != "y"){
        stop("Please change column names of presence points to \"x\" and \"y\" for long and lat respectively")
      }
    }
    
    #Methods for bias correction
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
                               as.imList = TRUE)
        
        w <- spatstat.geom::as.owin(imList[[1]])
        
        if(inherits(points, "data.frame")){
          pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
          Q <- spatstat.geom::pixelquad(pp)
        } 
        
        if(inherits(points, "ppp")){
          pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
          Q <- spatstat.geom::pixelquad(pp)
        }
        
        if(inherits(points, "quad")){
          pp <- spatstat.geom::ppp(x = points$data$x, y = points$data$x, window = w)
          Q <- spatstat.geom::pixelquad(pp)
        }
      }
      
      if(bias.correction == "weights"){

        if(inherits(bias.data, "data.frame")){
  
          imList <- imFromStack(covariates)
  
          w <- spatstat.geom::as.owin(imList[[1]])
          
          if(inherits(points, "data.frame")){
            pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          } 
          
          if(inherits(points, "ppp")){
            pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          }
          
          if(inherits(points, "quad")){
            pp <- spatstat.geom::ppp(x = points$data$x, y = points$data$x, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          }
          
          Qa <- replaceQAreas(Q = Q,
                              bias.data = bias.data, 
                              im = imList[[1]], 
                              positive = weight.bias.conf$positive,
                              kernel = weight.bias.conf$kernel,
                              sigma = weight.bias.conf$sigma,
                              varcov = weight.bias.conf$sigma, 
                              weights = weight.bias.conf$weights, 
                              edge = weight.bias.conf$edge,
                              zo.norm = weight.bias.conf$zo.norm)
          Q <- Qa
        }
        
        if(inherits(bias.data, "SpatRaster")){
          
          imList <- imFromStack(covariates)
          
          w <- spatstat.geom::as.owin(imList[[1]])
          
          if(inherits(points, "data.frame")){
            pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          } 
          
          if(inherits(points, "ppp")){
            pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          }
          
          if(inherits(points, "quad")){
            pp <- spatstat.geom::ppp(x = points$data$x, y = points$data$x, window = w)
            Q <- spatstat.geom::pixelquad(pp)
          }
          
          Qa <- replaceQAreas(Q = Q,
                              bias.data = bias.data, 
                              im = imList[[1]],
                              positive = weight.bias.conf$positive,
                              kernel = weight.bias.conf$kernel,
                              sigma = weight.bias.conf$sigma,
                              varcov = weight.bias.conf$sigma, 
                              weights = weight.bias.conf$weights, 
                              edge = weight.bias.conf$edge,
                              zo.norm = weight.bias.conf$zo.norm)
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
        pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = w)
        Q <- spatstat.geom::pixelquad(pp)
      }
      
      if(inherits(points, "quad")){
        pp <- spatstat.geom::ppp(x = points$data$x, y = points$data$x, window = w)
        Q <- Q <- spatstat.geom::pixelquad(pp)
      }
    }    
  
  if(parallel){
  
    doParallel::registerDoParallel(cores = cores)
    
    models <- foreach::foreach(i = seq_along(formulas)) %dopar% {
      
      form <- paste0("Q", formulas[i]) |> stats::as.formula()
      
      m <- spatstat.model::ppm.formula(form, 
                                      covariates = imList,
                                      correction = ppm.conf$correction,
                                      use.gam = ppm.conf$use.gam,
                                      method = ppm.conf$method,
                                      forcefit = ppm.conf$forcefit,
                                      improve.type = ppm.conf$improve.type,
                                      improve.args=ppm.conf$improve.args,
                                      prior.mean = ppm.conf$prior.mean,
                                      prior.var = ppm.conf$prior.var)
      
      return(m)
    }
  } else {
  
    models <- foreach::foreach(i = seq_along(formulas)) %do% {
  
      form <- paste0("Q", formulas[i]) |> stats::as.formula()
  
      m <- spatstat.model::ppm.formula(form, 
                                    covariates = imList,
                                    correction = ppm.conf$correction,
                                    use.gam = ppm.conf$use.gam,
                                    method = ppm.conf$method,
                                    forcefit = ppm.conf$forcefit,
                                    improve.type = ppm.conf$improve.type,
                                    improve.args=ppm.conf$improve.args,
                                    prior.mean = ppm.conf$prior.mean,
                                    prior.var = ppm.conf$prior.var)
      
      return(m)
    }
  }

  if(goodness.fit == "AIC"){
    perf <- sapply(models, spatstat.model::AIC.ppm) 
    perf.sort <- perf |> sort(decreasing = FALSE)
    top <- perf.sort[1:top.models]

    ids <- which(perf %in% top) |> rev()
  }
  
  if(goodness.fit == "BIC"){
    perf <- sapply(models, stats::BIC)
    perf.sort <- perf |> sort(decreasing = FALSE)
    top <- perf.sort[1:top.models]

    ids <- which(perf %in% top) |> rev()
  }

  if(goodness.fit == "logLik"){
    perf <- sapply(models, stats::logLik)
    perf.sort <- perf |> sort(decreasing = TRUE)
    top <- perf.sort[1:top.models]

    ids <- which(perf %in% top) |> sort(decreasing = TRUE)
  }
  
  gc(reset = TRUE)
  
  ret.list <- list(
    models = models[ids],
    call = list(formulas = formulas,
                covariates = covariates,
                bias.correction = bias.correction,
                bias.data = bias.data,
                top.models = top.models,
                goodness.fit = goodness.fit))
  
  class(ret.list) <- "ppmBatch"
  
  return(ret.list)
}
