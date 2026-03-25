#' @title Change the CRS of a set of coordinates, which are returned in a data.frame
#' @description Change the CRS of a set of coordinates, which are returned in a data.frame
#' @param lon A numeric vector containing the logitude or x coordinates of presence localities.
#' @param lat A numeric vector containing the latitude or y coordinates of presence localities.
#' @param origin.crs A character string specifying the CRS of the given point localities.
#' @param target.crs A character string specifying the new desired CRS in which to return the given point localities.
#' @return A two column data frame with the transformed point localities in the new CRS.
#' @examples
#' p <-  system.file("extdata", "points_WGS84.csv", package = "espatsmo") |>  read.csv()
#' 
#' p.itrf <- projectPoints(lon = p$x,
#'                         lat = p$y,
#'                         origin.crs = "EPSG:4326",
#'                         target.crs = "EPSG:6372")
#' @export

projectPoints <- function(lon = NULL,
                          lat = NULL,
                          origin.crs = "EPSG:4326",
                          target.crs = NULL){

  if(is.null(target.crs)){stop("Please provide a valid target CRS")}
  
  v <- data.frame(lon = lon, lat = lat) |> terra::vect()
  terra::crs(v) <- origin.crs
  
  vp <- terra::project(v, target.crs) |> terra::crds() |> data.frame()
  
  return(vp)
}