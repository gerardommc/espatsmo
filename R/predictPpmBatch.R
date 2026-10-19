#' @title Generate predictions from ppmBatch object
#' @description
#' Project a single model onto geographic space. To do this, the user should provide
#' a batch of fitted models, and the covariates in SpatRaster format with the same names as
#' used to fit the model.
#' @param object A model object of class ppmBatch
#' @param newdata A SpatRaster object with covariates in the same order as they were provided in the fitted model
#' @param se Logical to specify if prediction includes standard errors.
#' @param level Numeric fractional value to specify the confidence of the prediction standard error intervals
#' @param crs A character string to specify the CRS of the output raster
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of the covariates
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
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
#' models <- ppmBatchFit(points = pr,
#'                       covariates = r,
#'                       formulas = forms,
#'                       parallel = FALSE,
#'                       top.models = 3)
#' 
#' predictions <- predict(object = models, newdata = r)
#' 
#' @export
#' @method predict ppmBatch

predict.ppmBatch <- function(object = NULL,
                              newdata = NULL,
                              se = FALSE,
                              level = 0.95,
                              crs = NULL){
      
      if(is.null(newdata)){
            newdata <- object$call$covariates
      }
      
      if(is.null(crs)){
            crs <- terra::crs(newdata)
      }
      
      m <- object$models
      
      iml <- imFromStack(newdata)
      
      if(terra::nlyr(r) > 1){
            locs <- iml[[1]]
            w <- spatstat.geom::as.owin(locs)
      } 
      
      if(terra::nlyr(r) == 1){
            locs <- iml
            w <- spatstat.geom::as.owin(locs)
      }
      
      preds <- lapply(m, function(x){
            pred <- spatstat.model::predict.ppm(object = x,
                                                covariates = iml,
                                                locations = locs,
                                                window = w,
                                                type = "trend",
                                                se = se,
                                                level = level)
            
            if(se){
                  class(pred) <- "imList"
                  pred.r <- as.SpatRast(pred, target.crs = crs)
                  names(pred.r) <- names(pred)
                  return(pred.r)
            }
            
            if(!se){
                  pred.r <- terra::rast(pred)
                  terra::crs(pred.r) <- crs
                  return(pred.r)
            }
      })
      
      class(preds) <- "stackList"
      return(preds)
}
