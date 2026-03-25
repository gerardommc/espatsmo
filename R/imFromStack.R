#' @title Convert a raster stack to a list of image objects for spatstat use
#' @description
#' This function takes a single argument, which is a single or multiple band SpatRaster terra object to transform
#' into a list of im objects, and assigns an object class imList
#' @param r A terra SpatRaster, with 1:N bands
#' @return A list of im objects, assigned the class imList.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast()
#'  
#' r.iml <- imFromStack(r)
#' @export

imFromStack <- function(r = NULL){

    if(is.null(r)){
        stop("Cannot generate imList without a SpatRaster")
    }
    
    r.df <- as.data.frame(r, xy = T)
    
    ux = sort(unique(r.df$x)) #Extracting unique coordinates
    uy = sort(unique(r.df$y))
    nx = length(ux) #length of unique coordinates
    ny = length(uy)
    ref.cols = match(r.df$x, ux) #position of every data point
    ref.lines = match(r.df$y, uy)
    vec = rep(NA, max(ref.lines)*max(ref.cols)) # A vector with the length of data points
    ref.vec = (ref.cols - 1)*max(ref.lines) + ref.lines
    vec[ref.vec] = 1
    data.mask = matrix(vec, max(ref.lines), max(ref.cols))
    
    if(ncol(r.df) < 4){
        vec.all = rep(NA, max(ref.lines)*max(ref.cols))
        vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
        vec.all[ref.vec] = r.df[,3]
        lay <- spatstat.geom::im(matrix(vec.all, max(ref.lines), max(ref.cols),
                         dimnames = list(uy, ux)), xcol = ux, yrow = uy)
        return(lay)
    }
    
    if(ncol(r.df) > 3){

        im.list <- list()

        for(i in 3:ncol(r.df)){
            vec.all = rep(NA, max(ref.lines)*max(ref.cols))
            vec.ref = (ref.cols - 1)*max(ref.lines) + ref.lines
            vec.all[ref.vec] = r.df[,i]
            lay <- spatstat.geom::im(matrix(vec.all, max(ref.lines), max(ref.cols),
                             dimnames = list(uy, ux)), xcol = ux, yrow = uy)
            im.list[[i - 2]] <- lay
        }
        names(im.list) <- names(r)
        
        class(im.list) <- "imList"

        return(im.list)
    }
}