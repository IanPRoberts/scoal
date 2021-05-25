#' Vector sampling
#'
#' Samples a vector in the same way as base::sample except if the vector has
#' length 1, the single value is always returned
#'
#' @inheritParams base::sample

sample.vector <- function(x, size, replace = FALSE, prob = NULL){
  #Samples a vector. If the vector has length 1, always returns the single value
  if (length(x) == 1){
    x
  }
  else{
    sample(x, size, replace, prob)
  }
}
