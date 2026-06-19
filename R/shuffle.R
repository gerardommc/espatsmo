#' @title Shuffle a number sequence or a vector
#' @description
#' Randomise the order of a 1:n sequence of numbers 
#' @param n A numeric integer value
#' @param intCoercion A character string to inicate how non-integer values will be coerced to an integer, 
#' should be "ceiling of "floor"
#' @param set.seed Logical, used to state whether to use the system-wide seed or 
#' the one provided via the seed argument
#' @param seed A numeric value used as the randomising seed
#' @return A numeric integer vector with values 1:n but in random order 
#' @examples 
#' shuffle(10)
#' 
#' x <- rnorm(20)
#' 
#' shuffle(x)
#' @export

shuffle <- function(n = NULL,
                    intCoercion = "ceiling",
                    set.seed = FALSE,
                    seed = 432){
  
  if(set.seed){
    set.seed(seed)
  }
  
  if(!is.numeric(n)){
    stop("Please provide a numeric object")
  }
  
  if(length(n) == 1){
    
    if(!is.integer(n)){
      if(intCoercion == "ceiling"){
        n <- ceiling(n) |> as.integer()
      }
      
      if(intCoercion == "floor"){
        n <- floor(n) |> as.integer()
      }
    }    
    
    x <- 1:n
    x.1 <- numeric(n)
    
    for(i in 1:(n-1)){
      x.1[i] <- sample(x, 1)
      x <- x[-which(x == x.1[i])]
    }
    x.1[n] <- x
    return(x.1)
  }
  
  if(length(n) > 1){
    
    n1 <- length(n)
    
    x <- 1:n1
    x.1 <- numeric(n1)
    for(i in 1:(n1-1)){
      x.1[i] <- sample(x, 1)
      x <- x[-which(x == x.1[i])]
    }
    
    x.1[n1] <- x
    return(n[x.1])
  }
  
}

