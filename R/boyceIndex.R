#' @title Calculate the Boyce index
#' @description Measure predictive performance of a raster layer of evaluation points using the Boyce index.
#' @param raster A raster layer object in format terra:SpatRaster representing the predicitions of a statistical model.
#' @param points A two-column data frame where the first two columns have to be the x and y coordinates respectively.
#' @param r.points Numeric, an integer representing the number of random points to be drawn from the raster predictions to calculate predicted areas.
#' @param n.thresholds Numeric, an integer indicating the number of model prediction threholds.
#' @param thres.criteria Character string with values "regular" or "quantiles" to indicate how model predictions will be thresholded.
#' @param buf Numeric, the radius around testing presence points which will eliminate all areas further away (redundant if omission > 0).
#' @param log.transform A logical value to indicate whether to log-transform the raster values.
#' @param omission Numeric, double, representing the quantile of suitability values which will be used to exclude presence training points (my personal interpretation of Peterson et al. 2008).
#' @param save.plot Logical, to indicate whether to sace a pdf file of the simulated PartialROC analysis
#' @param plot.pars A list with entries 1) name (Character string), 2) width (plot width in inches) and 3) height (plot height in inches), to save the test plot.
#' @param set.seed A logical value to specify whther to use a random seed generator
#' @param seed An integer value to use as a seed.
#' @param remove.outer.points Logical, used to indicate iff points ling outside the study window are going to be removed from model run.
#' @return A list containing a data.frame with the estimated and observed frequencies and averae suitability.
#' values per class, along with the estimated Boyce index ad P-value (estimated from the correlation test).
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ] 
#' 
#' model <- ppmSingleFit(points= p, 
#'                      covariates = r, 
#'                      formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)",
#'                      as.ppmSingle = FALSE)
#' 
#' spatstat.model::summary.ppm(model)
#' 
#' v <- system.file("extdata", "ValidationSamples.csv", package = "espatsmo") |> read.csv()
#' 
#' pv <- p[v$Samples, ]
#' 
#' r.im <- imFromStack(r)
#' 
#' pred <- spatstat.model::predict.ppm(model, 
#'                                     covariates = r.im, 
#'                                     locations = r.im[[1]]) |> terra::rast() 
#' 
#' bi.re <- boyceIndex(raster = pred,
#'                    points = pv,
#'                    r.points = 5000,
#'                    n.thresholds = 25,
#'                    thres.criteria = "regular",
#'                    plot.pars = list(name = "BoyceIndex-regular.pdf", 
#'                                     width = 5, 
#'                                     height = 5),
#'                    save.plot = FALSE)
#' 
#' bi.qu <- boyceIndex(raster = pred,
#'                    points = pv,
#'                    r.points = 5000,
#'                    n.thresholds = 25,
#'                    thres.criteria = "quantiles",
#'                    plot.pars = list(name = "BoyceIndex-quantiles.pdf", 
#'                                     width = 5, 
#'                                     height = 5),
#'                    save.plot = FALSE)
#' 
#' @export

boyceIndex <- function(raster = NULL, 
                  points = NULL,
                  r.points = 5000,
                  n.thresholds = 25,
                  thres.criteria = "regular",
                  buf = NULL,
                  log.transform = FALSE,
                  omission = 0.05,
                  save.plot = TRUE, 
                  plot.pars = list(name = "BoyceIndex.pdf", width = 5, height = 5),
                  set.seed = FALSE,
                  seed = 432,
                  remove.outer.points = TRUE){
  
  if(set.seed){
    set.seed(seed)
  }

  if(!thres.criteria %in% c("regular", "quantiles")){
    stop("Please specify a valid value for thres.criteria, either \"regular\", or \"quantiles\" ")
  }

  if(remove.outer.points){
    nas <- terra::extract(raster, points[, 1:2], ID = FALSE)[, 1]
    points <- points[!is.na(nas), ]
  }
  
  if(!is.null(buf)){
    p <- terra::vect(as.matrix(points[, -3]))
    bu <- terra::buffer(p, width = buf)
    
    bur <- terra::rasterize(bu, raster)
    
    raster <- terra::mask(raster, bur)
  }
  
  if(omission > 0){
    vals <- terra::extract(raster, points[, 1:2])[,2]
    q <- stats::quantile(vals, 1-omission, na.rm = TRUE)
    om.r <- raster < q
    om.r <- terra::classify(om.r, rcl = matrix(c(-Inf, 0, NA), ncol = 3, byrow = TRUE))
    raster <- terra::mask(raster, om.r)
    
    points <- points[vals < q, ]
  }
  
  if(log.transform){
    raster <- log(raster + 0.01)
  }
  
  r <- ZeroOneNorm(raster)
  #r <- round(r, log10(n.thresholds) + 1)
  
  thres <- seq(0, 1, len = n.thresholds + 1)
  
  #Thresholding suitability layer
  
  r.xy <- r |> terra::as.data.frame(xy = TRUE)
  
  if(thres.criteria == "regular"){
    thres <- seq(0, 1, len = n.thresholds + 1)
  }
  
  if(thres.criteria == "quantiles"){
    thrs <- seq(0, 1, len = n.thresholds + 1)
    thres <- stats::quantile(x = r.xy[, 3], 
                             probs = thrs, 
                             na.rm = TRUE)
  }
  
  rt1 <- (r >= thres[1]) * (r <= thres[2])
  
  for(i in 2:n.thresholds){
    rt1 <- c(rt1, (r >= thres[i]) * (r <= thres[i+1]))
  }
  
  rt.rec <- terra::classify(rt1, matrix(c(0, NA)))
  
  r.bins <- r * rt.rec
  
  r.means <- terra::global(r.bins, mean, na.rm = TRUE)
  
  r.areas <- terra::global(rt1, mean, na.rm = TRUE)
  
  samp <- sample(x = 1:nrow(r.xy),
                       size = r.points,
                       replace = FALSE,
                       prob = NULL) |> sort()
  
  sample.points <- r.xy[samp, c("x", "y")]
  
  P.i <- (terra::extract(rt1, points, ID = FALSE) |> colSums())/nrow(points)
  
  E.i <- (terra::extract(rt1, sample.points, ID = FALSE) |> colSums())/nrow(sample.points)
  
  F.i <- P.i/E.i
  
  Bo.data <- data.frame(F.i = F.i, P.i = P.i, E.i = E.i, Suit.Class.Mean = r.means$mean) |> stats::na.omit()#Xtraer valores promedio de cada región
  Bo.data <- Bo.data[is.finite(Bo.data$F.i), ]

  cr <- tryCatch(
                stats::cor.test(x = Bo.data$F.i, 
                                y = Bo.data$Suit.Class.Mean, 
                                method = "spearman"),
                error = function(x){paste0("I cannot calculate the Boyce Index with the given sample size")})
  
  
  if(inherits(cr, "character")){
      Boyce <- cr
      P.value <- NA
      save.plot <- FALSE
      warning("Skipping Boyce index plot")
  } else {
      Boyce <- cr$estimate
      P.value <- cr$p.value
  }
  
  if(save.plot){
    grDevices::pdf(plot.pars$name, width = plot.pars$width, height = plot.pars$height)
    graphics::plot(F.i ~ Suit.Class.Mean, 
                   data = Bo.data, 
                   col = "red",
                   pch = 20,
                   xlab = "Suitability Class",
                   ylab = "Predicted/Expected",
                   main = paste0("Boyce index = ", round(Boyce, 2), ", P = ", round(P.value, 2)))
    grDevices::dev.off()
  }
  
  ret.list <- list(Frequencies = Bo.data,
                   Boyce.index = Boyce,
                   P.value = P.value)
  
  return(ret.list)
}
