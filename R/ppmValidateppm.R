#' @title Generate validation statistics
#' @description
#' Quantify the predictive capacity of a point process model by data partitioning 
#' @param model A model object of class ppm, ppmSingle or ppmBatch
#' @param method A character string specifying whther the method to use is the partial ROC ("proc"),
#'  the Boyce index ("boyce") or both.
#' @param proc.pars A list containing the arguments passed to `partialROC`
#' @param boyce.pars A list containing the arguments passed to `boyceIndex`
#' @param part.pars A list containing the arguments passed to `boyceIndex`
#' @param ppm.conf A list contating the arguentos passed to `ppm`, as in `ppmSingleFit`
#' @param crs A character string using authority codes for coordinate reference systems
#' @return A summary of the validation statistics generated
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' s <-  system.file("extdata", "RandomSamples.csv", package = "espatsmo") |>  read.csv()
#' 
#' pr <- p[s$Samples, ]
#' 
#' model <- ppmSingleFit(points = pr,
#'                       covariates = r,
#'                       formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)",
#'                       as.ppmSingle = FALSE)
#' 
#' val <- ppmValidate(model = model,
#'                    method = c("proc", "boyce"),
#'                    crs = "EPSG:6372")
#' 
#' @export
#' @method ppmValidate ppm

ppmValidate.ppm <- function(model = NULL,
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
                            ppm.conf =  list(correction = "border", 
                                             use.gam = TRUE, 
                                             method = "mpl", 
                                             forcefit = TRUE,
                                             improve.type = "none", 
                                             improve.args = list(), 
                                             prior.mean = NULL, 
                                             prior.var = NULL),
                            crs = NULL){
  
  r <- model$covariates |> as.SpatRast(target.crs = crs)
  names(r) <- names(model$covariates)
  
  Q <- model$Q
  p <- data.frame(x = Q$data$x,
                  y = Q$data$y)
  
  part <- spatialPartition(covariates = r,
                           points = p,
                           no.blocks = part.pars$no.blocks,
                           part.criteria = part.pars$part.criteria,
                           mask.criteria = part.pars$mask.criteria,
                           seed = part.pars$seed)
  
  ## Configure bias.data here
  
  iml1 <- imFromStack(part$part1$covariates)
  if(inherits(iml1, "imList")){
    w1 <- spatstat.geom::as.owin(iml1[[1]])
  }
  if(inherits( iml1, "im")){
    w1 <- spatstat.geom::as.owin(iml1)
  }
  ppp1 <- spatstat.geom::ppp(x = part$part1$points$x,
                             y = part$part1$points$y,
                             window = w1)
  Q1 <- spatstat.geom::pixelquad(ppp1)
  
  iml0 <- imFromStack(part$part0$covariates)
  if(inherits(iml0, "imList")){
    w0 <- spatstat.geom::as.owin(iml0[[1]])
  }
  if(inherits(iml0, "im")){
    w0 <- spatstat.geom::as.owin(iml0)
  }
  ppp0 <- spatstat.geom::ppp(x = part$part0$points$x,
                             y = part$part0$points$y,
                             window = w0)
  Q0 <- spatstat.geom::pixelquad(ppp0)
  
  m1 <- spatstat.model::ppm.quad(Q = Q1,
                                 trend = formula(model),
                                 covariates = iml1,
                                 correction = ppm.conf$correction,
                                 use.gam = ppm.conf$use.gam, 
                                 method = ppm.conf$method, 
                                 forcefit = ppm.conf$forcefit,
                                 improve.type = ppm.conf$improve.type, 
                                 improve.args = ppm.conf$improve.args, 
                                 prior.mean = ppm.conf$prior.mean, 
                                 prior.var = ppm.conf$prior.var)
  
  m0 <-   m1 <- spatstat.model::ppm.quad(Q = Q0,
                                         trend = formula(model),
                                         covariates = iml0,
                                         correction = ppm.conf$correction,
                                         use.gam = ppm.conf$use.gam, 
                                         method = ppm.conf$method, 
                                         forcefit = ppm.conf$forcefit,
                                         improve.type = ppm.conf$improve.type, 
                                         improve.args = ppm.conf$improve.args, 
                                         prior.mean = ppm.conf$prior.mean, 
                                         prior.var = ppm.conf$prior.var)
  
  pred1 <- spatstat.model::predict.ppm(m1, locations = iml0[[1]], covariates = iml0) |> terra::rast()
  pred0 <- spatstat.model::predict.ppm(m0, locations = iml1[[1]], covariates = iml1) |> terra::rast()
  
  retList <- list()
  
  if("proc" %in% method){

    nproc1 <- proc.pars$plot.pars$name
    nproc1 <-  base::substr(x = nproc1, 
                            start = 1, 
                            stop = nchar(nproc1) - 4)
    
    plot.pars.proc1 <- proc.pars$plot.pars
    plot.pars.proc1$name <- nproc1

    proc1 <- partialROC(raster = pred0,
                        points = part$part1$points,
                        p.points = proc.pars$p.points,
                        r.points =  proc.pars$r.points,
                        iterations =  proc.pars$iterations,
                        buf =  proc.pars$buf,
                        log.transform =  proc.pars$log.transform,
                        n.thresholds =  proc.pars$n.thresholds,
                        thres.criteria =  proc.pars$thres.criteria,
                        omission =  proc.pars$omission,
                        save.plot =  proc.pars$save.plot,
                        plot.pars =  plot.pars.proc1,
                        set.seed =  proc.pars$set.seed,
                        seed =  proc.pars$seed)
    
    nproc0 <- proc.pars$plot.pars$name
    nproc0 <-  base::substr(x = nproc0, 
                            start = 1, 
                            stop = nchar(nproc0) - 4)
    
    plot.pars.proc0 <- proc.pars$plot.pars
    plot.pars.proc0$name <- nproc0
    
    proc0 <- partialROC(raster = pred1,
                        points = part$part0$points,
                        p.points = proc.pars$p.points,
                        r.points =  proc.pars$r.points,
                        iterations =  proc.pars$iterations,
                        buf =  proc.pars$buf,
                        log.transform =  proc.pars$log.transform,
                        n.thresholds =  proc.pars$n.thresholds,
                        thres.criteria =  proc.pars$thres.criteria,
                        omission =  proc.pars$omission,
                        save.plot =  proc.pars$save.plot,
                        plot.pars =  plot.pars.proc0,
                        set.seed =  proc.pars$set.seed,
                        seed =  proc.pars$seed)
    
    retList$proc <- list(part1 = proc1,
                         part0 = proc0)
  }
  
  if("boyce" %in% method){

    nboy1 <- boyce.pars$plot.pars$name
    nboy1 <- base::substr(x = nboy1, 
                            start = 1, 
                            stop = nchar(nboy1) - 4)
    
    plot.pars.boy1 <- boyce.pars$plot.pars
    plot.pars.boy1$name <- nboy1

    boy1 <- boyceIndex(raster = pred0,
                       points = part$part1$points,  
                       r.points = boyce.pars$r.points,
                       n.thresholds = boyce.pars$n.thresholds,
                       thres.criteria = boyce.pars$thres.criteria,
                       buf = boyce.pars$buf,
                       log.transform = boyce.pars$log.transform,
                       omission = boyce.pars$omission,
                       save.plot = boyce.pars$save.plot,
                       plot.pars = plot.pars.boy1,
                       set.seed = boyce.pars$set.seed,
                       seed = boyce.pars$seed)
    
    nboy0 <- boyce.pars$plot.pars$name
    nboy0 <- base::substr(x = nboy0, 
                            start = 1, 
                            stop = nchar(nboy0) - 4)
    
    plot.pars.boy0 <- boyce.pars$plot.pars
    plot.pars.boy0$name <- nboy1
    
    boy0 <-  boyceIndex(raster = pred1,
                        points = part$part0$points,  
                        r.points = boyce.pars$r.points,
                        n.thresholds = boyce.pars$n.thresholds,
                        thres.criteria = boyce.pars$thres.criteria,
                        buf = boyce.pars$buf,
                        log.transform = boyce.pars$log.transform,
                        omission = boyce.pars$omission,
                        save.plot = boyce.pars$save.plot,
                        plot.pars = plot.pars.boy0,
                        set.seed = boyce.pars$set.seed,
                        seed = boyce.pars$seed)
    
    retList$boyce <- list(boyce1 = boy1,
                          boyce0 = boy0)
  }
  
  retList$partition <-  list(part1 = list(points = ppp1, 
                                          covariates = iml1),
                             part0 = list(points = ppp0, 
                                          covariates = iml0))
  
  return(retList)
  
}
