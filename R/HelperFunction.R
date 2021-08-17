## bound number to (0, 1)
bound01 <- function(x, r = 1e-10){
  xx <- x
  xx[x < r] <- r
  xx[x > 1-r] <- 1-r
  return(as.numeric(xx))
}

## bound number to (0, Inf)
bound <- function(x, r = 0.001){
  xx <- x
  xx[x < r] <- r
  return(as.numeric(xx))
}
