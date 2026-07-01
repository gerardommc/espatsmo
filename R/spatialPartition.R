#' @title Split a geographical area into random spatial blocks
#' @description
#' Generate partitions of a geographical area to split into two separate data sets to
#'  measure a fitted model's predictive capacity
#' @param covariates A SpatRaster object from the terra package
#' @param points A two-column data.frame containing the x and y coordinates of presente localities
#' @param no.blocks An integer numeric value indicating the number of blocks into which data will be divided
#' @param part.criteria A character string indicating whether the blocks will be generated from random samples or from regular intervals
#' @param mask.criteria A character string indicating whether the mask criteria to split the raster covariates will be generated 
#' from evenly spaced points or randomly generated localities.
#' @param seed A numeric integer value o be used as the seed for randomisation, only used if mask and partition criteria are set to be random.
#' @return An object of class spatPart containing a list of three elements. The first two elements are the partitions, each of which contains the 
#' masked covariates and presence points, The third element, contains the randomising seed.
#' @examples 
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#'
#' part <- spatialPartition(covariates = r,
#'                          points = pr, 
#'                          no.blocks = 100,
#'                          part.criteria = "random",
#'                          mask.criteria = "random",
#'                          seed = 432)
#'                          
#' plot(part, partition = 1)
#' plot(part, partition = 0)
#' @export


spatialPartition <- function(covariates = NULL,
                             points = NULL,
                             no.blocks = 50,
                             part.criteria = "regular",
                             mask.criteria = "regular",
                             seed = 432){
  set.seed(seed)
  
  r <- covariates
  p <- points
  
  if(part.criteria == "random" & mask.criteria == "random"){
    xy <- terra::crds(r) |> terra::as.data.frame(xy = TRUE)
    
    xy <- xy[sample(1:nrow(xy), size =  no.blocks, replace = FALSE, prob = NULL),]
    
    v <- xy |> as.matrix() |> terra::vect()
    terra::crs(v) <- terra::crs(r)
    
    vor <- terra::voronoi(v)
    
    vor$id <- 1:nrow(xy)
    
    zo <- rep(c(1, 0), length.out = nrow(xy))
    
    vor$zo<- zo[shuffle(length(zo), set.seed = TRUE, seed = seed)]
    
    vor.1 <- vor[vor$zo == 1]
    vor.0 <- vor[vor$zo == 0]
    
    r1 <- terra::mask(x = r, mask = vor.1, touches = TRUE)
    r0 <- terra::mask(x = r, mask = vor.0, touches = TRUE)
  }
  
  if(part.criteria == "regular" & mask.criteria == "random"){
    xy <- terra::crds(r) |> terra::as.data.frame(xy = TRUE)
    
    xy <- xy[sample(1:nrow(xy), size= no.blocks, replace = FALSE, prob = NULL),]
    
    v <- xy |> as.matrix() |> terra::vect()
    terra::crs(v) <- terra::crs(r)
    
    vor <- terra::voronoi(v)
    
    vor$id <- 1:nrow(xy)
    
    vor$zo <- vor$id %% 2
    
    vor.1 <- vor[vor$zo == 1]
    vor.0 <- vor[vor$zo == 0]
    
    r1 <- terra::mask(x = r, mask = vor.1, touches = TRUE)
    r0 <- terra::mask(x = r, mask = vor.0, touches = TRUE)
  }
  
  if(part.criteria == "regular" & mask.criteria== "regular"){
    
    ex <- terra::ext(r) 
    
    xmin <- ex$xmin
    xmax <- ex$xmax
    
    ymin <- ex$ymin
    ymax <- ex$ymax
    
    ratio <- (xmax - xmin)/(ymax - ymin)
    
    len.x <- (sqrt(no.blocks) * ratio) |> ceiling() 
    len.y <- (sqrt(no.blocks)) |> ceiling() 
    
    b <- terra::rast(nrows = len.y,
                      ncol = len.x,
                      xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax,
                      crs = terra::crs(r))
    
    ri <- terra::init(b, "row")
    ci <- terra::init(b, "col")
    
    board <- (ri + ci) %% 2
    
    vor <- terra::as.polygons(board)
    
    vor$zo<- vor$row
    
    vor.1 <- vor[vor$zo == 1]
    vor.0 <- vor[vor$zo == 0]
    
    r1 <- terra::mask(x = r, mask = vor.1, touches = TRUE)
    r0 <- terra::mask(x = r, mask = vor.0, touches = TRUE)
  }
  
  if(part.criteria == "random" & mask.criteria== "regular"){
    
    ex <- terra::ext(r) 
    
    xmin <- ex$xmin
    xmax <- ex$xmax
    
    ymin <- ex$ymin
    ymax <- ex$ymax
    
    ratio <- (xmax - xmin)/(ymax - ymin)
    
    len.x <- (sqrt(no.blocks) * ratio) |> ceiling() 
    len.y <- (sqrt(no.blocks)) |> ceiling() 
    
    xy <- expand.grid(x = seq(xmin, xmax, length.out = len.x),
                      y = seq(ymin, ymax, length.out = len.y))
    
    v <- xy |> as.matrix() |> terra::vect()
    terra::crs(v) <- terra::crs(r)
    
    vor <- terra::voronoi(v)
    
    zo <- rep(c(1, 0), length.out = len.x * len.y)
    
    vor$zo<- zo[shuffle(length(zo), set.seed = TRUE, seed = seed)]
    
    vor.1 <- vor[vor$zo == 1]
    vor.0 <- vor[vor$zo == 0]
    
    r1 <- terra::mask(x = r, mask = vor.1, touches = TRUE)
    r0 <- terra::mask(x = r, mask = vor.0, touches = TRUE)
  }
  
  nas1 <- terra::extract(r1, p, ID = F, raw = T) |> data.frame()
  nas0 <- terra::extract(r0, p, ID = F, raw = T) |> data.frame()
  
  p1 <- p[!is.na(nas1[,1]),]
  p0 <- p[!is.na(nas0[,1]),]
  
  retList <- list(part0 = list(covariates = r0,
                               points = p0),
                  part1 = list(covariates = r1,
                               points = p1),
                  seed = seed)
  
  class(retList) <- "spatPart"
  
  return(retList)
}
