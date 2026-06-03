#' @title Calculate the Partial ROC statistic proposed by (Peterson et al. 2008).
#' @description Validate the predictive capacity of a spatial model to predict a set of independent locations.
#' @param raster A raster layer object in format terra:SpatRaster representing the predicitions of a statistical model.
#' @param points A two-column data frame where the first two columns have to be the x and y coordinates respectively.
#' @param p.points Numeric, double representing the proportion of validation points used in each iteration.
#' @param r.points Numeric, an integer representing the number of random points to be drawn from the raster predictions to calculate predicted areas.
#' @param n.thresholds Numeric, an integer indicating the number of model prediction threholds.
#' @param thres.criteria Character string with values "regular" or "quantiles" to indicate how model predictions will be thresholded.
#' @param iterations Numeric, integer the number of times the random sampling is repeated.
#' @param buf Numeric, the radius around testing presence points which will eliminate all areas further away (redundant if omission > 0).
#' @param log.transform A logical value to indicate whether to log-transform the raster values.
#' @param omission Numeric, double, representing the quantile of suitability values which will be used to exclude presence training points (my personal interpretation of Peterson et al. 2008).
#' @param save.plot Logical, to indicate whether to sace a pdf file of the simulated PartialROC analysis.
#' @param plot.pars A list with entries 1) name (Character string), 2) width (plot width in inches) and 3) height (plot height in inches), to save the test plot.
#' @param set.seed A logical value to specify whther to use a random seed generator.
#' @param seed An integer value to use as a seed.
#' @return A list containing a data.frame with the simulated area ratios, the mean average ratio and the P-value (probability that area-ratios < 1).
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
#'                    r.points = 5000,
#'                    n.thresholds = 100,
#'                    thres.criteria = "regular",
#'                    plot.pars = list(name = "PartialROC-regular.pdf", 
#'                                     width = 5, 
#'                                     height = 5),
#'                    save.plot = FALSE)
#' 
#' proc1 <- partialROC(raster = pred,
#'                    points = valid.points,
#'                    r.points = 5000,
#'                    n.thresholds = 100,
#'                    thres.criteria = "quantiles",
#'                    plot.pars = list(name = "PartialROC-quantiles.pdf", 
#'                                     width = 5, 
#'                                     height = 5),
#'                    save.plot = FALSE)
#' 
#' @export

partialROC <- function(raster, 
                       points, 
                       p.points = 0.5,
                       r.points = 1000,
                       iterations = 39,
                       buf = NULL, 
                       log.transform = TRUE,
                       n.thresholds = 100,
                       thres.criteria = "regular",
                       omission = 0, 
                       save.plot = TRUE, 
                       plot.pars = list(name = "PartialROC.pdf", width = 5, height = 5),
                       set.seed = FALSE,
                       seed = 432
                      ){
  
  `%do%` <- foreach::`%do%`

  if(set.seed){
    set.seed(seed)
  }

  if(!thres.criteria %in% c("regular", "quantiles")){
    stop("Please specify a valid value for thres.criteria, either \"regular\", or \"quantiles\" ")
  }
  
  if(!is.null(buf)){
    p <- terra::vect(as.matrix(points[, -3]))
    bu <- terra::buffer(p, width = buf)
    
    bur <- terra::rasterize(bu, raster)
    
    raster <- terra::mask(raster, bur)
  }
  
  if(omission > 0){
    vals <- terra::extract(raster, points[, 1:2])[,2]
    q <- stats::quantile(x = vals, 
                         probs = 1-omission, 
                         na.rm = TRUE)
    om.r <- raster < q
    om.r <- terra::classify(om.r, rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = TRUE))
    raster <- terra::mask(raster, om.r)
    
    points <- points[vals < q, ]
  }
  
  if(log.transform){
    raster <- log(raster + 0.01)
  }
  
  r <- ZeroOneNorm(raster)
  r <- round(r, ceiling(log10(n.thresholds)))
  
  #Thresholding suitability layer

  r.xy <- r |> as.data.frame(xy = TRUE)

  if(thres.criteria == "regular"){
    thres <- seq(0, 1, len = n.thresholds + 1)
  }

  if(thres.criteria == "quantiles"){
    thrs <- seq(0, 1, len = n.thresholds + 1)
    thres <- stats::quantile(x = r.xy[, 3], 
                             probs = thrs, 
                             na.rm = TRUE)
  }
  
  samp <- base::sample(x = 1:nrow(r.xy),
                       size = r.points,
                       replace = FALSE,
                       prob = NULL) |> sort()
  
  sample.points <- r.xy[samp, 3]
  
  pred.thrs <- foreach::foreach(i = seq_along(thres), .combine = cbind) %do% {
    sample.points >= thres[i]
  }  #Verificar que r.thr sea 
  
  pres.values <- terra::extract(r, points, ID = FALSE)[, 1]
  
  pres.thrs <- foreach::foreach(i = seq_along(thres), .combine = cbind) %do% {
    pres.values >= thres[i]
  } 
  
  area.pred <- colMeans(pred.thrs)
  
  dArea <- area.pred[1:n.thresholds] - area.pred[2:(n.thresholds + 1)]
  
  mp <- matrix(0, nrow = iterations, ncol = length(area.pred))
  mr <- matrix(0, nrow = iterations, ncol = length(area.pred))
  
  areas <- matrix(0, nrow = iterations, ncol = 4)
  colnames(areas) <- c("iteration", 
                       "Area.pres",
                       "Area.rand",
                       "PartialROC")
  
  for(i in 1:iterations){
    samp.pres <- sample(1:nrow(pres.thrs), size = nrow(points) * p.points, replace = FALSE) |> sort()
    samp.rand <- sample(1:nrow(pred.thrs), size = nrow(points) * p.points, replace = FALSE) |> sort()
    
    pres.thr <- pres.thrs[samp.pres, ]
    rand.thr <- pred.thrs[samp.rand, ]
    
    
    #Calculating proportion of predicted points  
    means.pres <- colMeans(pres.thr)
    means.rand <- colMeans(rand.thr)
    
    #Calculating areas
    
    rects.pres <- dArea * means.pres[-(n.thresholds  + 1)]
    rects.rand <- dArea * means.rand[-(n.thresholds  + 1)]
    
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
  
  ret.list <- list(simulations = areas,
                   Ratio = rat,
                   P = P)

  return(ret.list)
}
