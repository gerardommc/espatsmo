#' @title Fit a log-Gaussian Cox point process model
#' @description Fit a log-Gaussian Cox point process model using INLA, with the posibility of controlling observer bias
#' by modifying area weights or the distribution of background (quadrature) points.
#' @param points A two-clumn data.frame, ppp or quad object containing the presence localities 
#' @param covariates A SpatRaster or imList containing the environmental covariates 
#' @param  formula character, speciftying the right hand side of a linear combination of covariate predictors
#' @param  offset character, specifying the name of the covariate to use as the population at risk or the population density from which presence points were drawn.
#' @param  no.bkgd numeric integer, specifying the number of background or quadratre points used to fit the model 
#' @param  group.bins numeric integer, specifying the numbr off bins into which the distance to presence points will be divided in order to fit the autoregressive random effect for distance to presence points
#' @param bias.data  A two-column data frame, containing the coordinates of the sampling localities, 
#' SpatRaster or imList object containing the spatial variability of the observation effort
#' @param bias.correction character, with values 'weights' or 'background', 
#' specifying the method used to control the effect of sampling bias
#' @param covariate.zero Character, indicating the variable name representing observer bias within the model formula. Must be a variable contained in the covariates raster stack
#' @param covariance.func character, specifying the name of the covariance function for the gaussian process, either "pcmatern", or "matern". Please check how these options should be configured 
#' via the `prior.conf` argument.
#' @param weight.bias.conf list, containing the following elements: 1) positive, 2) kernel, 3) sigma, 4) varcov, 5) weights, 6) edge
#' @param prior.conf list, with two different specifications depending on the covariance function selected. 
#' If `covariance.func = "pcmatern"`, the list must provide values for: 1) alpha = 2,
#' 2) prior.range = c(100, 0.01), 3) prior.sigma = c(0.05, 0.01)).
#' If `covariance.func = "matern"`, the list must provide values for:
#' 1) sigma0 = 0.05, 2) range0 = 0.2, 3) alpha = 2, 4) kappa0 = sqrt(8)/20,
#' 5) B.tau = c(min = -1, max = 1), 6) B.kappa = c(max = 0, min = -1), 7)theta.prior.mean = c(0, 0),
#' 8)theta.prior.prec = c(0.1, 1).
#' For in-depth descrition of the meaning of these parameters, the reades should consult INLA guidelines
#' for both covariance functions and particularly adjust the distance-dependent arguments for the 
#' coordinate units of the working point process.
#' @param  dist.ar Logical, indicating whether to include an autorregressive random effect for distance to presence points
#' @param  mesh.par list, which must provide values for: 1) edge = 50, 
#' and 2) offset = 50; used to configure the projector mesh to integrate the random field.
#' @param  coordinates character, used to specify the x and y units in which the covariance 
#' matrix will be computed, and takes two possible values, "m" or "km".
#' @param  dist.units character, used to specify the distance units in which the covariance 
#' matrix will be computed, and takes two possible values, "m" or "km".
#' @param  weight.units = "km",
#' @param  verbose logical, to configure whether the model fitting process will bee printed on screen in real time.
#' @param  inla.mode character, to specify whether the "experimental" of NULL model will be used
#' @return list of class `ppmLGCP` containing the fitted model, call, mesh and predictions for fixed and random effects.
#' @examples
#' r <- system.file("extdata", "ChelsaBio.tif", package = "espatsmo") |>  terra::rast() |> scale()
#' 
#' p <- system.file("extdata", "points.csv", package = "espatsmo") |>  read.csv()
#' 
#' bias <- system.file("extdata", "Target-group.tif", package = "espatsmo") |> terra::rast()
#' 
#' model <- ppmLGCP(points= p, 
#'                  covariates = r, 
#'                  formula = "~ bio1 + bio2 + bio12 + I(bio1^2) + I(bio2^2) + I(bio12^2)", 
#'                  dist.ar = FALSE,
#'                  weight.units = "km",
#'                  coordinates = "m")
#' @export

ppmLGCP <- function(points = NULL,
                    covariates = NULL,
                    formula = NULL,
                    offset = NULL,
                    no.bkgd = NULL,
                    group.bins = 50,
                    bias.data = NULL,
                    bias.correction = NULL, #options = "background", "weights",
                    covariate.zero = NULL,
                    covariance.func = "pcmatern", #options = "matern", "pcmatern"
                    weight.bias.conf = list(positive = TRUE, kernel = "gaussian",
                                            sigma = NULL, varcov = NULL, 
                                            weights = NULL, edge = TRUE),
                    prior.conf = list(alpha = 2,
                                      prior.range = c(10000, 0.01),
                                      prior.sigma = c(0.05, 0.01)),
                    mesh.par = list(edge = 5000, offset = 10000), # for the distance between points in the spde, units are x the grain size in raster layers
                    coordinates = "m",
                    dist.units = "km",
                    weight.units = "km",
                    dist.ar = FALSE,
                    verbose = FALSE,
                    inla.mode = "experimental"){
  
  INLA::inla.setOption(inla.mode = inla.mode)

  `%do%` <- foreach::`%do%`

  if(is.null(points) | is.null(covariates) | is.null(formula)){
    stop("Please provide a valid set of points, covariates and one model formula")
  }

  if(names(points)[1] != "x" | names(points)[2] != "y"){
    stop("Please provide a set of points as a two column data.frame with names \"x\" and \"y\", or a valid ppp or quad object")
  }
  
  if(inherits(points, "ppp")){
    points <- data.frame(x = points$x, y = points$y)
  }

  if(inherits(points, "quad")){
    points <- data.frame(x = points$data$x, y = points$data$y)
  }
  
  #Function to calculate haversine distance if polar coordinates are used
  dist.coords <- function(x, y, units = dist.units){ 
    if(dist.units == "km"){
      r <- 6378.1370
    }
    
    if(dist.units == "m"){
      r <- 63781370
    }
    
    xrad <- x * pi/180
    yrad <- y * pi/180
    
    xrad <- x * pi/180
    yrad <- y * pi/180
    
    dx <- xrad[2] - xrad[1]
    dy <- yrad[2] - yrad[1]
    
    a <- sin(dx / 2)^2 + sin(dy / 2)^2 * cos(xrad[1]) * cos(xrad[2])
    c <- 2 * r * asin(sqrt(a))
    return(c)
  }
  
  cov.df <- as.data.frame(covariates, xy = TRUE)
  
  zeroes <- terra::classify(covariates[[1]], rcl = matrix(c(-Inf, Inf, 0), ncol = 3))
  ones <- zeroes + 1
  
  #Conffiguring background locations
  iml <- imFromStack(covariates)
  win <- spatstat.geom::as.owin(iml[[1]])
  p.pp <- spatstat.geom::ppp(x = points$x, y = points$y, window = win)
  
  #Calculating weights
  Q <- spatstat.geom::pixelquad(p.pp)
  
  beg <- Q$w |> length() - iml[[1]][] |>length() + 1
  en <- Q$w |> length()
  
  Q$w[beg:en] <- max(Q$w[beg:en])
  
  Area.size <- Q$w[beg:en] |> sum()
  
  no.units <- (beg:en) |> length()
  
  if(!is.null(bias.correction)){
    if(bias.correction == "background" & is.null(no.bkgd)){
      stop("Background bias correction is incompatible without a specific number of background points")
    }
    
    if(is.null(bias.data)){
      stop("Please provide an observation bias data set")
    }
  }
  
  if(dist.ar & is.null(group.bins)){
    stop("Autorregressive distance effects requires a number of distance bins to be provided")
  }
  
  #Creating spde mesh  
  hull <- cov.df[, c("x", "y")] |> as.matrix() |> fmesher::fm_nonconvex_hull()
  
  ux <- cov.df$x |> sort() |> unique()
  dx <- ux[2] - ux[1]
  
  #Background locations without bias correction
  if(is.null(bias.correction)){
    area.weights <- iml[[1]]
    area.weights[] <- Q$w[beg:en]
    weights.r <- area.weights |> terra::rast()
    
    if(!is.null(no.bkgd)){
      s <- sample(1:nrow(cov.df), no.bkgd) |> sort()
      wei <- Area.size/no.bkgd
      
      quad.xy <- cov.df[s, c("x", "y")]
      
    } else {
      
      mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                       max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                       boundary = hull,
                                       cutoff = dx,
                                       offset = mesh.par$offset)
      
      quad.xy <- mesh$loc[, 1:2] |> data.frame()
      names(quad.xy) <- c("x", "y")
      
      nas.quad <- ! terra::extract(covariates[[1]], quad.xy, ID = FALSE)[, 1] |> is.na()
      
      quad.xy <- quad.xy[nas.quad, ]
      
      wei <- Area.size/nrow(quad.xy)
    }
  }
  
  if(!is.null(bias.correction)){
    #Background data with bias correction
    #Bias correction based on are weights
    if(bias.correction == "weights"){
      if(inherits(bias.data, "data.frame")){
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = iml[[1]],
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma, 
                            weights = weight.bias.conf$weights, 
                            edge = weight.bias.conf$edge)
        
        area.weights <- iml[[1]]
        area.weights[] <- Qa$w[beg:en]
        weights.r <- area.weights |> terra::rast()
        
        # Adjust here for the number of background points and total area size
        
        if(!is.null(no.bkgd)){
          s <- sample(1:nrow(cov.df), no.bkgd) |> sort()
          
          wei <- terra::extract(weights.r, cov.df[s, c("x", "y")])[,2]
          
          nas <- wei |> is.na()
          
          s <- s[!nas]
          
          mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                           max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                           boundary = hull,
                                           cutoff = dx,
                                           offset = mesh.par$offset)
          
          quad.xy <- cov.df[s, c("x", "y")]
          
          wei <- wei[!nas] * nrow(cov.df)/nrow(quad.xy)
        } else {
          mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                           max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                           boundary = hull,
                                           cutoff = dx,
                                           offset = mesh.par$offset)
          
          quad.xy <- mesh$loc[, 1:2] |> data.frame()
          names(quad.xy) <- c("x", "y")
          
          nas.quad <- ! terra::extract(covariates[[1]], quad.xy, ID = FALSE)[, 1] |> is.na()
          
          quad.xy <- quad.xy[nas.quad, ]
          
          wei <- terra::extract(weights.r, quad.xy)[, 2] * nrow(cov.df)/nrow(quad.xy)
        }
        
      }
      
      if(inherits(bias.data, "SpatRaster")){
        
        bias.data <- terra::resample(bias.data, covariates[[1]])
        
        Qa <- replaceQAreas(Q = Q,
                            bias.data = bias.data, 
                            im = iml[[1]], 
                            positive = weight.bias.conf$positive,
                            kernel = weight.bias.conf$kernel,
                            sigma = weight.bias.conf$sigma,
                            varcov = weight.bias.conf$sigma, 
                            weights = weight.bias.conf$weights, 
                            edge = weight.bias.conf$edge)
        
        area.weights <- iml[[1]]
        area.weights[] <- Qa$w[beg:en]
        weights.r <- area.weights |> terra::rast()
        
        if(!is.null(no.bkgd)){
          s <- sample(1:nrow(cov.df), no.bkgd) |> sort()
          
          wei <- terra::extract(weights.r, cov.df[s, c("x", "y")])[,2]
          
          nas <- wei |> is.na()
          
          s <- s[!nas]
          
          mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                           max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                           boundary = hull,
                                           cutoff = dx,
                                           offset = mesh.par$offset)
          
          quad.xy <- cov.df[s, c("x", "y")]
          
          wei <- wei[!nas] * nrow(cov.df)/nrow(quad.xy)
        } else {
          mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                           max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                           boundary = hull,
                                           cutoff = dx,
                                           offset = mesh.par$offset)
          
          quad.xy <- mesh$loc[, 1:2] |> data.frame()
          names(quad.xy) <- c("x", "y")
          
          nas.quad <- ! terra::extract(covariates[[1]], quad.xy, ID = FALSE)[, 1] |> is.na()
          
          quad.xy <- quad.xy[nas.quad, ]
          
          wei <- terra::extract(weights.r, quad.xy)[, 2] * nrow(cov.df)/nrow(quad.xy)
        }
      }
    }
    
    #Bias correction based on location of  background data
    if(bias.correction == "background"){
      if(inherits(bias.data, "data.frame")){
        bias.ppp <- spatstat.geom::ppp(x = bias.data$x, y = bias.data$y, window = win)
        dens.r <- spatstat.explore::density.ppp(bias.ppp, 
                                                positive = weight.bias.conf$positive,
                                                kernel = weight.bias.conf$kernel,
                                                sigma = weight.bias.conf$sigma,
                                                varcov = weight.bias.conf$varcov, 
                                                weights = weight.bias.conf$weights, 
                                                edge = weight.bias.conf$edge) |> terra::rast()
        
        dens.df <- as.data.frame(dens.r, xy = TRUE)
        
        area.weights <- iml[[1]]
        area.weights[] <- Q$w[beg:en]
        weights.r <- area.weights |> terra::rast()
        
        s <- sample(1:nrow(dens.df), no.bkgd, prob = dens.df[, 3]) |> sort()
        
        mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                         max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                         boundary = hull,
                                         cutoff = dx,
                                         offset = mesh.par$offset)
        
        quad.xy <- cov.df[s, c("x", "y")]
        
        wei <- max(Q$w) * nrow(cov.df)/nrow(quad.xy)
      }
      
      if(inherits(bias.data, "SpatRaster")){
        bias.data <- terra::resample(bias.data, covariates[[1]]) |> ZeroOneNorm()
        
        bias.df <- as.data.frame(bias.data, xy = TRUE)
        
        area.weights <- iml[[1]]
        area.weights[] <- Q$w[beg:en]
        weights.r <- area.weights |> terra::rast()
        
        s <- sample(1:nrow(bias.df), no.bkgd, prob = bias.df[, 3]) |> sort()
        
        mesh <- fmesher::fm_mesh_2d_inla(loc.domain = cov.df[, c("x", "y")], 
                                         max.edge = c(mesh.par$edge, (mesh.par$edge)*5),
                                         boundary = hull,
                                         cutoff = dx,
                                         offset = mesh.par$offset)
        
        quad.xy <- cov.df[s, c("x", "y")]
        
        wei <- max(Q$w) * nrow(cov.df)/nrow(quad.xy)
      }
    }
  }
  
  if(coordinates == "deg" & weight.units == "km"){
    wei <- wei * 1.0E4
  }
  
  if(coordinates == "deg" & weight.units == "m"){
    wei <- wei * 1.0E7
  }
  
  if(coordinates == "m" & weight.units == "km"){
    wei <- wei/1.0E6
  }
  
  if(dist.ar){
    pres.d <- foreach::foreach(i = 1:nrow(points), .combine = c) %do% {
      if(coordinates == "deg"){
        x1 <- points$x[i]
        x2 <- points$x[-i]
        
        y1 <- points$y[i]
        y2 <- points$y[-i]
        pd <- sapply(seq_along(y2), function(ii){dist.coords(x = c(x1, x2[ii]),
                                                             y = c(y1, y2[ii]))})
        return(pd[which.min(pd)[1]])
      }else{
        if(coordinates == "m" & dist.units == "km"){
          pd <- sqrt((points$x[i]/1000 - points$x[-i]/1000)^2 + (points$y[i]/1000 - points$y[-i]/1000)^2)
          return(pd[which.min(pd)[1]])
        }else{      
          pd <- sqrt((points$x[i] - points$x[-i])^2 + (points$y[i] - points$y[-i])^2)
          return(pd[which.min(pd)[1]])
        }
      }
    }
    
    quad.d <- foreach::foreach(i = 1:nrow(quad.xy), .combine = c) %do% {
      if(coordinates == "deg"){
        x1 <- quad.xy$x[i]
        x2 <- points$x
        
        y1 <- quad.xy$y[i]
        y2 <- points$y
        qd <- sapply(seq_along(y2), function(ii){dist.coords(x = c(x1, x2[ii]),
                                                             y = c(y1, y2[ii]))})
        return(qd[which.min(qd)[1]])
      }else{
        if(coordinates == "m" & dist.units == "km"){
          qd <- sqrt((quad.xy$x[i]/1000 - points$x/1000)^2 + (quad.xy$y[i]/1000 - points$y/1000)^2)
          return(qd[which.min(qd)[1]])
        }else{
          qd <- sqrt((quad.xy$x[i] - points$x)^2 + (quad.xy$y[i] - points$y)^2)
          return(qd[which.min(qd)[1]])
        }
      }
    }
    
    gc(reset = TRUE)
  }

  #Distance to quadrature points

  points.v <- data.frame(lon = points$x, lat = points$y) |> as.matrix()  |> terra::vect()
  terra::crs(points.v) <- terra::crs(covariates)

  points.r <- terra::rasterize(points.v, covariates, fun = "count")
  points.df <- as.data.frame(points.r, xy = TRUE)
  
  #Covariance priors
  if(covariance.func == "pcmatern"){
      spde <- INLA::inla.spde2.pcmatern(mesh,
                                  alpha = prior.conf$alpha,
                                  prior.range = prior.conf$prior.range,
                                  prior.sigma = prior.conf$prior.sigma)
  }

  if(covariance.func == "matern"){
    spde <- INLA::inla.spde2.matern(mesh, 
                                    alpha = prior.conf$alpha,
                                    B.tau = prior.conf$B.tau,
                                    B.kappa = prior.conf$B.kappa, 
                                    theta.prior.mean = prior.conf$theta.prior.mean,
                                    theta.prior.prec = prior.conf$theta.prior.prec)
  }

  #Putting data together
  #Presence A matrix
  presence.df <- terra::extract(covariates, points.v, ID = FALSE)
  if(dist.ar){
    presence.df <- data.frame(points[, c("x", "y")], presence.df, dist = pres.d)
  } else {
    presence.df <- data.frame(points[, c("x", "y")], presence.df)
  }
  nas.pres<- apply(presence.df, 2, function(x){which(is.na(x))}) |> unlist() |> unique()
  if(length(nas.pres) > 0){
    presence.df <- presence.df[-nas.pres, ]
  }
  presence.sp <- presence.df |> sf::st_as_sf(coords = c("x", "y")) |> sf::as_Spatial()
  A_point <- INLA::inla.spde.make.A(mesh, loc = presence.sp)
  
  #quadrature A matrix
  quad.data <- terra::extract(covariates, quad.xy, ID = FALSE)
  if(dist.ar){
    quad.data <- data.frame(quad.xy, quad.data, dist = quad.d)
  } else {
    quad.data <- data.frame(quad.xy, quad.data)
  }
  nas.quad <- apply(quad.data, 2, function(x){which(is.na(x))}) |> unlist() |> unique()
  if(length(nas.quad) > 0){
    quad.data <- quad.data[-nas.quad, ]
  }
  quad.sp <- quad.data |> sf::st_as_sf(coords = c("x", "y")) |> sf::as_Spatial()
  A_quad <- INLA::inla.spde.make.A(mesh, loc = quad.sp)
  
  #Predictors matrix
  covs.df <- as.data.frame(covariates, xy = TRUE) |> stats::na.omit()
  
  if(dist.ar){
    if(dist.units == "km"){
      bkgd.dist <- terra::distance(points.r, method = "geo", unit = "m") / 1000
    }
    if(dist.units == "m"){
      bkgd.dist <- terra::distance(points.r, method = "geo", unit = "m")
    }
    covs.df$dist <- terra::extract(bkgd.dist, covs.df[, c("x", "y")])$count
    final.form <- paste0(formula, " + dist")
  } else {
    final.form <- formula
  }

  covs.sp <- covs.df |> sf::st_as_sf(coords = c("x", "y")) |> sf::as_Spatial()
  A_pred <- INLA::inla.spde.make.A(mesh, loc = covs.sp)
  

  if(!is.null(offset)){
    if(length(offset)==1){
      final.form <- paste0(final.form, " + ", paste0("offset(log(", offset, "))"))
    }else{
      final.form <- paste0(final.form, " + ", paste(paste0("offset(log(", offset, "))"), collapse = " + "))
    } #Revisar cómo evitar que dist esté en los efectos fijos
  }
  
  sp.mm <- stats::model.matrix(formula(final.form), presence.df) |> data.frame()
  quad.mm <- stats::model.matrix(formula(final.form), quad.data) |> data.frame()
  
  if(!is.null(covariate.zero)){
    covs.df1 <- covs.df
    covs.df1[[covariate.zero]] <- 0
    pred.mm <- stats::model.matrix(formula(final.form), covs.df1) |> data.frame()
  }else{
    pred.mm <- stats::model.matrix(formula(final.form), covs.df) |> data.frame()
  }
  
  if(!is.null(offset)){
    sp.mm[, offset] <- presence.df[, offset]
    quad.mm[, offset] <- quad.data[, offset]
    pred.mm[, offset] <- covs.df[, offset]
  }
  
  #Building data stacks
  stack.sp <- INLA::inla.stack(data = list(y = 1, e = 0),
                         A = list(A_point, 1),
                         tag = "sp",
                         effects =list(
                           list(i = 1:mesh$n),
                           sp.mm)
  )
  
  if(length(nas.quad) > 0 & length(wei) > 1){
    stack.quad <- INLA::inla.stack(data = list(y = 0, e = wei[-nas.quad]),
                             A = list(A_quad, 1),
                             effects =list(
                               list(i = 1:mesh$n),
                               quad.mm))
  }else{
    stack.quad <- INLA::inla.stack(data = list(y = 0, e = wei),
                             A = list(A_quad, 1),
                             effects =list(
                               list(i = 1:mesh$n),
                               quad.mm))
  }
  
  stack.pred <- INLA::inla.stack(data = list(y = 1, e = 1),
                           A = list(A_pred, 1), tag= "pred",
                           effects = list(
                             list(i = 1:mesh$n),
                             pred.mm)
  )
  
  stack.all <- INLA::inla.stack(stack.sp, stack.quad, stack.pred)
  
  which.dist <- which(names(sp.mm) == "dist")
  
  #Formulas for models with autorregressive distance random effects
  if(dist.ar){
    
    which.dist <- which(names(sp.mm) == "dist")
    
    if(!is.null(offset)){
      which.offset <- which(names(sp.mm) == offset)
      if(length(offset) == 1){
        formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm)[-c(which.dist, which.offset)], collapse = " + "), paste0(" + offset(log(", offset,"))"), " + f(inla.group(dist, n = ", group.bins, ", method = \"quantile\"), 
                     model = \"rw1\", scale.model = TRUE) + f(i, model = spde)")  |> stats::as.formula()
      }else{
        formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm)[-c(which.dist, which.offset)], collapse = " + "), " + ", paste(paste0("offset(log(", offset,"))"), collapse = " + "), " + f(inla.group(dist, n =", group.bins, ", method = \"quantile\"), 
                     model = \"rw1\", scale.model = TRUE) + f(i, model = spde)") |> stats::as.formula()
      }
    }else{
      formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm)[-which.dist], collapse = " + "), " + f(inla.group(dist, n = ", group.bins, ", method = \"quantile\"), 
                     model = \"rw1\", scale.model = TRUE) + f(i, model = spde)") |> stats::as.formula()
    }
  } 
  
  #Formulas for models without autorregressive distance random effects
  if(!dist.ar){
    if(!is.null(offset)){
      which.offset <- which(names(sp.mm) == offset)
      if(length(offset) == 1){
        formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm)[-c(which.offset)], collapse = " + "), paste0(" + offset(log(", offset,"))", "+ f(i, model = spde)"))  |> stats::as.formula()
      }else{
        formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm)[-c(which.offset)], collapse = " + "), " + ", paste(paste0("offset(log(", offset,"))"), collapse = " + "), " + f(i, model = spde)") |> stats::as.formula()
      }
    }else{
      formula.inla <- paste0("y ~ 0 + ", paste(names(sp.mm), collapse = " + "), " + f(i, model = spde)") |> stats::as.formula()
    }
  }
  
  data.stack <- INLA::inla.stack.data(stack.all)
  fit <- INLA::inla(formula = formula.inla,
                    family = "poisson",
                    data = data.stack,
                    control.predictor = list(A = INLA::inla.stack.A(stack.all), compute = TRUE),
                    E = data.stack$e, control.compute = list(dic = TRUE),
                    control.fixed = list(expand.factor.strategy = "inla"),
                    verbose = verbose)
  
  #xtracting results to return raster layer
  idx <- INLA::inla.stack.index(stack.all, 'pred')$data
  
  model.preds <- fit$summary.fitted.values[idx, c("mean", "sd", "0.025quant",  "0.975quant", "mode")]
  
  if(!is.null(covariate.zero)){covs.df[[covariate.zero]] <- 0}
  
  #Model with covariates
  model.preds.r <- terra::rast(data.frame(covs.df[, c("x", "y")], model.preds))
  
  #Model predictions without spatial effects and offsets
  mm.pred <- stats::model.matrix(stats::as.formula(formula), covs.df)
  preds.mean <- apply(mm.pred, 1, function(x){
    sum(x * fit$summary.fixed$mean)
  })
  
  nas.cov <- apply(covs.df, 2, function(x){which(is.na(x))}) |> unlist() |> unique()
  
  if(length(nas.cov) > 0){
    preds.mean.r <- terra::rast(data.frame(covs.df[-nas.cov, c("x", "y")], preds.mean)) |> exp()
  } else {
    preds.mean.r <- terra::rast(data.frame(covs.df[, c("x", "y")], preds.mean)) |> exp()
  }
  
  ret.data <- list(model = fit,
                   predictions = model.preds.r,
                   base.predictions = preds.mean.r,
                   points = points,
                   quad.points = quad.xy,
                   covariates = covariates,
                   mesh = mesh,
                   call = list(formula = formula,
                               offset = offset,
                               no.bkgd = no.bkgd,
                               group.bins = group.bins,
                               bias.correction = bias.correction,
                               weight.bias.conf = weight.bias.conf,
                               prior.conf = prior.conf,
                               mesh.par = mesh.par, 
                               coordinates = coordinates, #the alternative is m
                               dist.units = dist.units, #the alternative is "m" for m
                               weight.units = weight.units,
                               inla.mode = inla.mode))
  
  class(ret.data) <- c("ppmLGCP", covariance.func)
  
  return(ret.data)
}