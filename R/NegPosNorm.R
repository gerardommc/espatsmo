#' @title Force a SpatRaster or imList to have values between (-1, 1)
#' @description
#' This function takes three arguments, r ust be a SpatRaster or imList
#' target.crs, is the desired crs if the input object r is an imList and 
#' returns a SpatRaster or imList with values (-1, 1) whilst retaining all 
#' its internal variability
#' @param r A single or multi-band SpatRaster or imList
#' @param target.crs Character string, indicating the coordinate reference system of the returned object, only used when the class of r is imList
#' @param as.imList Logical, whether to return the object as an imList
#' @return A SpatRaster or imList with all its values compresed between 0-1.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' r.np <- NegPosNorm(r = r,
#'                    target.crs = "EPSG:6372",
#'                    as.imList = F)
#' 
#' @export

NegPosNorm <- function(r = NULL, 
                       target.crs = NULL,
                       as.imList = F){

  if(inherits(r, "imList")){
    r <- as.SpatRast(r)
    terra::crs(r) <- target.crs
  }

  if(!inherits(r, "SpatRaster")){
    stop("Cannot normalise with the supplied object, please provide a SpatRaster")
  }

  ma <- terra::global(r, max, ID = F)[, 1]
  r1 <- r - ma/2
  r2 <- r1/ma

  if(as.imList){
    im.list <- imFromStack(r2)
    return(im.list)
  }

  if(!as.imList){
    return(r2)
  }
}
