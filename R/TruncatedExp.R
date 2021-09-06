#' The Truncated Exponential Distribution
#'
#' Density, distribution function and random generation for the exponential
#' distribution with rate \code{rate} truncated above by \code{c}
#'
#' @param x,q vector of quantiles
#' @param n number of observations
#' @param rate vector of rates
#' @param c truncation cutoff
#'
#' @export

dexp.trunc <- function(x, c = Inf, rate = 1){
  #Density function
  dexp <- rate * exp(-rate * x) / (1 - exp(-rate * c))
  dexp
}

#' @rdname dexp.trunc

pexp.trunc <- function(q, c = Inf, rate = 1){
  #Distribution function
  pexp <- (1 - exp(-rate * q)) / (1 - exp(-rate*c))
  pexp
}

#' @rdname dexp.trunc

rexp.trunc <- function(n, c = Inf, rate = 1){
  #i.i.d. sample of size n
  rexp <- runif(n)
  rexp <- - (1/rate) * log(1 - (1 - exp(-rate * c)) * rexp)
  rexp
}
