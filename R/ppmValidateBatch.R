#' @title Generate validation statistics
#' @description
#' Quantify the predictive capacity of a point process model by data partitioning 
#' @param model A model object of class ppmBatch
#' @param method A character string specifying whther the method to use is the partial ROC ("proc"),
#'  the Boyce index ("boyce") or both.
#' @param proc.pars A list containing the arguments passed to `partialROC`
#' @param boyce.pars A list containing the arguments passed to `boyceIndex`
#' @param part.pars A list containing the arguments passed to `boyceIndex`
#' @param crs A character string using authority codes for coordinate reference systems
#' @return A summary of the validation statistics generated
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' resp <- system.file("extdata", "Exponents.csv", package = "espatsmo") |> read.csv()
#' 
#' compat <- findCompatibles(covariates = r,
#'                           thres = 0.6,
#'                           max.comb = 3)
#' 
#' forms <- getPolyFormulas(respDF = resp, 
#'                          compatMat = compat)
#' 
#' models <- ppmBatchFit(points = p,
#'                       covariates = r,
#'                       formulas = forms,
#'                       parallel = FALSE,
#'                       top.models = 3)
#' 
#' val <- ppmValidate(model = models,
#'                    method = c("proc", "boyce"),
#'                    crs = "EPSG:6372")
#' 
#' @export
#' @method ppmValidate ppmBatch

ppmValidate.ppmBatch <- function(model = NULL,
                                 method = c("proc", "boyce"),
                                 proc.pars = list(p.points = 0.5,
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
                                                  seed = 432,
                                                  remove.outer.points = TRUE),
                                 boyce.pars = list(r.points = 5000,
                                                   n.thresholds = 25,
                                                   thres.criteria = "regular",
                                                   buf = NULL,
                                                   log.transform = FALSE,
                                                   omission = 0.05,
                                                   save.plot = TRUE,
                                                   plot.pars = list(name = "BoyceIndex.pdf", width = 5, height = 5),
                                                   set.seed = FALSE,
                                                   seed = 432,
                                                   remove.outer.points = TRUE),
                                part.pars = list(no.blocks = 50,
                                                 part.criteria = "regular",
                                                 mask.criteria = "regular",
                                                 seed = 432),
                                crs = NULL){
  
    ppms <- model$models
    
    vals <- lapply(seq_along(ppms), function(i){
      
      if("proc" %in% method){
      nproc0 <- proc.pars$plot.pars$name
      nproc <-  base::substr(x = nproc0, 
                             start = 1, 
                             stop = nchar(nproc0) - 4)
      
      proc.parsi <- proc.pars
      proc.parsi$plot.pars$name <- list(name = paste0(nproc, "-", i, ".pdf"))
      }
      
      if("boyce" %in% method){
        nboy0 <- boyce.pars$plot.pars$name
        nboy <-  base::substr(x = nboy0, 
                              start = 1, 
                              stop = nchar(nboy0) - 4)
        
        boyce.parsi <- boyce.pars
        boyce.parsi$plot.pars$name <- list(name = paste0(nboy, "-", i, ".pdf"))
      }
      

      
      val <- ppmValidate.ppm(model = ppms[[i]],
                             method = method,
                             proc.pars = proc.parsi,
                             boyce.pars = boyce.parsi,
                             part.pars = part.pars,
                             crs = crs)
      return(val)
    })
    
    names(vals) <- paste0("Model.", seq_along(vals))
    
    return(vals)
}
