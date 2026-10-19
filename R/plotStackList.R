#' @title Plot a stackList of batch model predictions
#' @description
#' Plot a batch of predictions froma ppmBatch object
#' @param x A stackList object
#' @param model Numeric value indicating the model number to plot
#' @return A r graphics plot
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
#' plot(predictions)
#' 
#' @export
#' @method plot stackList

plot.stackList <- function(x, model = NULL){

      if(is.null(model)){
            for(i in seq_along(x)){
                  terra::plot(x[[i]])
            }
      }

      if(!is.null(model)){

            if(model %in% seq_along(x$models)){
                  stop("Cannot extract model ids which are not contained within the list of fitted models, please provide valid id numbers")
            }
            
            if(length(model) > 1){
                  for(i in model){
                        terra::plot(x[[i]])
                  }
            }

            if(length(model) == 1){
                  terra::plot(x[[model]])
            }
      }
}
