#' @title Force minimum and maximum values of a raster layer to be in the 0-1 interval
#' @description
#' The function requires three arguments: (1) r, should be a SpatRaster, imList or im object to be processed; 
#' (2) as.imList, whether the generated layer should be returned as an imList or im object; and (3) the crs of
#' the final SpatRaster, if the desired output is a SpatRaster, which is only used iff the input data is an im or
#' imList
#' @param r A SpatRast or imList object.
#' @param as.imList Logical, whether to return the zero-one normalised raster as an imList object.
#' @param target.crs A character string specifying the CRS of the returned object.
#' @return A SpatRaster or imList object.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' r.zon <- ZeroOneNorm(r = r,
#'                      as.imList = F,
#'                      target.crs = "EPSG:6372")
#' @export 

ZeroOneNorm <- function(r = NULL, 
                        as.imList = F,
                        target.crs = NULL){

  if(is.null(r)){
    stop("Please provide a valid SpatRaster or imList")
  }

  if(inherits(r, "imList")){
    r <- as.SpatRast(r)
  }

  if(inherits(r, "im")){
    r <- terra::rast(r)
  }

  mi <- terra::global(r, min, na.rm = T)[, 1]
  ma <- terra::global(r, max, na.rm = T)[, 1]
  r.norm <- (r - mi)/(ma - mi)

  if(!as.imList){
    if(inherits(r, "imList") | inherits(r, "im")){
      terra::crs(r.norm) <- target.crs
    }
    return(r.norm)
  }

  if(as.imList){
    im.list <- imFromStack(r.norm)
    return(im.list)
  }
}
