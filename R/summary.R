#'@title Generate summaries from model objects
#'@description These functions require a range of inputs depending on the class of the supplied
#'objects. For objects of class ppmLGCP, the summary is for either fixed or random effects
#'

summary <- function(model, effects = NULL){
  UseMethod("summary", model)
}

summary.ppmLGCP <- function(model = NULL, effects = NULL){
  if(is.null(effects)){
    pl <- list(`Fixed effects` = model$model$summary.fixed,
               `Random effects` = m$model$summary.hyperpar)
    print(pl, quote = F)
  }
  
  if(!is.null(effects)){
      if(effects == "fixed"){
          print(model$model$summary.fixed, quote = F)
      }
      
      if(effects == "random"){
          print(model$model$summary.hyperpar, quote = F)
      }
    }
}

summary.ppmBatch <- function(model = NULL, id = NULL){
  spatstat.model::summary.ppm(object = model$models[[id]])
}

summary.ppmSingle <- function(model = NULL){
  spatstat.model::summary.ppm(object = model$model)
}


