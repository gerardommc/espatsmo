#' @title Generate predictions from ppmSingle object
#' @description
#' Project a single model onto geographic space. To do this, the user should provide
#' a fitted model, and the covariates in SpatRaster format with the same names as
#' used to fit the model.
#' @param object A model object of class ppmSingle
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
#' #Model with random sampling
#' 
#' model <- ppmSingleFit(points = pr,
#'                       covariates = r,
#'                       formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)",
#'                       as.ppmSingle = FALSE)
#' 
#' predictions <- predict(object = model, newdata = r)
#' 
#' plot(predictions)
#' 
#' @export
#' @method predict ppmSingle

predict.ppmSingle <- function(object = best,
                              newdata = r,
                              se = TRUE,
                              level = 0.95,
                              crs = NULL){
      
      if(is.null(newdata)){
            stop("Please provide a SpatRaster with equal covariate names and the supplied object")
      }
      
      if(is.null(crs)){
            crs <- terra::crs(newdata)
      }
      
      
      m <- object$model
      
      iml <- imFromStack(newdata)
      
      if(terra::nlyr(r) > 1){
            locs <- iml[[1]]
            w <- spatstat.geom::as.owin(locs)
      } 
      
      if(terra::nlyr(r) == 1){
            locs <- iml
            w <- spatstat.geom::as.owin(locs)
      }
      
      pred <- spatstat.model::predict.ppm(object = m, 
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
      
}
