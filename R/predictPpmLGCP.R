#' @title Generate predictions from a ppmLGCP fixed terms.
#' @description
#' Project a ppmLGCP model onto geographic space from the fixed terms only. To do this, the user should provide
#' a fitted model, and the covariates in SpatRaster format with the same names as
#' used to fit the model.
#' @param object = A model object of class ppmLGCP
#' @param newdata = A SpatRaster object with covariates in the same order as they were provided in the fitted model
#' @param probs = The posterior probability quantiles to be returned by predict.ppmve
#' @param crs A character string to specify the CRS of the output raster
#' @return Returns a single or multiple band SpatRaster object, representing point intensity as a function of the covariates
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#'
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#'
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#'
#' model <- ppmLGCP(points= p, 
#'                 covariates = r, 
#'                 formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                 dist.ar = FALSE,
#'                 weight.units = "km",
#'                 coordinates = "m")
#' 
#' predictions <- predict(object = model, newdata = r, probs = c(0.0275, 0.5, 0.975))
#' 
#' plot(predictions)
#' 
#' @export
#' @method predict ppmLGCP

predict.ppmLGCP <- function(object = NULL,
                            newdata = NULL,
                            probs = c(0.025, 0.5, 0.975),
                            crs = NULL){

  
        `%do%` <- foreach::`%do%`      
  
        if(is.null(newdata)){
          newdata <- object$covariates
        }
        
        newdat.df <-  terra::as.data.frame(newdata, xy = T)
        
        nas.ids <- apply(newdat.df, 2, function(x){which(is.na(x))}) |> unlist() |> unique()
        
        if(length(nas.ids) > 0){
          newdat.df <- newdat.df[-nas.ids, ]
        }
        
        mod.mat <- stats::model.matrix(object = stats::as.formula(model$call$formula),
                                       data = newdat.df)
        
        means <- object$model$summary.fixed$mean
        sds <- object$model$summary.fixed$sd
        
        coefs <- sapply(probs, function(x){
                                 stats::qnorm(p = x,
                                              mean = means,
                                              sd = sds)})
        
        preds <- foreach::foreach(i = seq_along(probs), .combine = cbind) %do% {
          apply(mod.mat, 1, function(x){
            sum(x * coefs[, i])
          })
        }
        
        pred.r <- data.frame(newdat.df[, c("x", "y")], preds) |> terra::rast() |> exp()
        names(pred.r) <- paste0(probs)
        crs(pred.r) <- crs
        
        return(pred.r)
}
