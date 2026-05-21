#' @title Calculate the Partial ROC statistic proposed by (Peterson et al. 2008).
#' @description Validate the predictive capacity of a spatial model to predict a set of independent locations.
#' @param raster A raster layer object in format terra:SpatRaster representing the predicitions of a statistical model
#' @param points A two-column data frame where the first two columns have to be the x and y coordinates respectively
#' @param p.points Numeric, double representing the proportion of validation points used in each iteration
#' @param iterations Numeric, integer the number of times the random sampling is repeated
#' @param buf Numeric, the radius around testing presence points which will eliminate all areas further away (redundant if omission > 0)
#' @param omission Numeric, double, representing the quantile of suitability values which will be used to exclude presence training points (my personal interpretation of Peterson et al. 2008).
#' @param save.plot Logical, to indicate whether to sace a pdf file of the simulated PartialROC analysis
#' @param plot.pars A list with entries 1) name (Character string), 2) width (plot width in inches) and 3) height (plot height in inches), to save the test plot.
#' @return A data.frame containing the simulated area ratios, and saves the plot of the areas and its median.
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
#'                      bias.data = bias, #Data frame with sampling localities or raster layer
#'                      bias.correction = "weights",
#'                      as.ppmSingle = FALSE)
#' 
#' spatstat.model::summary.ppm(model)
#' 
#' r.im <- imFromStack(r)
#' 
#' pred <- spatstat.model::predict.ppm(model, covariates = r.im, locations = r.im[[1]]) |> terra::rast()
#' 
#' valid.points <- system.file("extdata", "Valid_points.csv", package = "espatsmo") |> read.csv()
#' 
#' proc <- partialROC(raster = pred,
#'                    points = valid.points,
#'                    plot.pars = list(name = "PartialROC.pdf", 
#'                                     width = 5, 
#'                                     height = 5),
#'                    save.plot = FALSE)
#' 
#' @export

partialROC <- function(raster, 
                       points, 
                       p.points = 0.5, 
                       iterations = 39,
                       buf = NULL, 
                       omission = 0, 
                       save.plot = TRUE, 
                       plot.pars = list(name = "PartialROC.pdf", width = 5, height = 5)){
  
  if(!is.null(buf)){
    p <- terra::vect(as.matrix(points[, -3]))
    bu <- terra::buffer(p, width = buf)
    
    bur <- terra::rasterize(bu, raster)
    
    raster <- terra::mask(raster, bur)
  }
  
  raster <- log(raster + 0.1)
  
  if(omission > 0){
    vals <- terra::extract(raster, points[, 1:2])[,2]
    q <- stats::quantile(vals, 1-omission, na.rm = TRUE)
    om.r <- raster < q
    om.r <- terra::classify(om.r, rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = TRUE))
    raster <- terra::mask(raster, om.r)
    
    points <- points[vals < q, ]
  }
  
  r <- ZeroOneNorm(raster)
  r <- round(r, 2)
  
  thres <- seq(0, 1, len = 101)
  
  #Thresholding suitability layer
  r.thr <- r >= thres
  area.pred <- terra::global(r.thr, mean, na.rm = TRUE)$mean
  dArea <- area.pred[1:100] - area.pred[2:101]
  
  points.r <- as.data.frame(r.thr, xy = TRUE)[, c("x", "y")]
  
  mp <- matrix(0, nrow = iterations, ncol = length(area.pred))
  mr <- matrix(0, nrow = iterations, ncol = length(area.pred))
  
  areas <- matrix(0, nrow = iterations, ncol = 4)
  colnames(areas) <- c("iteration", 
                       "Area.pres",
                       "Area.rand",
                       "PartialROC")
  
  for(i in 1:iterations){
    samp.pres <- sample(1:nrow(points), size = nrow(points) * p.points, replace = FALSE)
    samp.rand <- sample(1:nrow(points.r), size = nrow(points) * p.points, replace = FALSE)
    
    pres.thr <- terra::extract(r.thr, points[samp.pres, 1:2])[, -1]
    rand.thr <- terra::extract(r.thr, points.r[samp.rand, ])[, -1]
    
    
    #Calculating proportion of predicted points  
    means.pres <- colMeans(pres.thr)
    means.rand <- colMeans(rand.thr)
    
    #Calculating areas
    
    rects.pres <- dArea * means.pres[-101]
    rects.rand <- dArea * means.rand[-101]
    
    Area.pres <- sum(rects.pres)
    Area.rand <- sum(rects.rand)
    Area.ratio <- Area.pres / Area.rand
    
    areas[i, ] <- c(iteration = i, 
                    Area.pres = Area.pres,
                    Area.rand = Area.rand,
                    PartialROC = Area.ratio)
    mp[i, ] <- means.pres
    mr[i, ] <- means.rand
  }
  
  areas  <- as.data.frame(areas)
  
  rat <- round(mean(areas$PartialROC), 2)
  P <- with(areas, length(which(PartialROC < 1))/iterations)
  
  if(save.plot){
  grDevices::pdf(paste0(plot.pars$name), width = plot.pars$width, height = plot.pars$height)
    plot(area.pred, colMeans(mp), xlab = "% Area predicted", 
         ylab = "1 - Omission error", col = "grey95", type = "l", 
         xlim = c(0, 1), ylim = c(0, 1), main = paste0("AUC ratio = ", rat, "\n P = ", P))
    graphics::lines(area.pred, colMeans(mr), col = "grey95", lty = 2, type = "l")
    for(j in 1:iterations){
      graphics::lines(area.pred, mp[j, ], col = "grey95", lwd = 0.25, type = "l")
      graphics::lines(area.pred, mr[j, ], col = "grey95", lty = 2, lwd = 0.25, type = "l")
    }
    graphics::lines(area.pred, colMeans(mp, na.rm = TRUE), col = "red", lwd = 1.5, type = "s")
    graphics::lines(area.pred, colMeans(mr, na.rm = TRUE), type = "s")
  grDevices::dev.off()
}

  
  areas <- as.data.frame(areas)
  
  return(areas)
}
