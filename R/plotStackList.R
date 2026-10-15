#' @title Plot a stackList of batch model predictions
#' @description
#' Plot a batch of predictions froma ppmBatch object
#' @param x A stackList object
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

plot.stackList <- function(x){
      for(i in seq_along(x)){
            plot(x[[i]])
      }
}
