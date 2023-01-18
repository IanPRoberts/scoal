#' Forward-in-time Migration Matrix
#'
#' Computes a forward-in-time migration matrix given a backwards-in-time migration matrix and effective population sizes
#'
#' @param migration.matrix Backwards-in-time migration matrix
#' @param effective.population Vector of effective population sizes
#'
#' @return Forward-in-time migration matrix
#'
#' @export

forward.migration.matrix <- function(migration.matrix, effective.population = rep(1, dim(migration.matrix)[1])){
  n.deme <- dim(migration.matrix)[1]
  f.mig.mat <- matrix(0, n.deme, n.deme)
  for (i in 1 : n.deme){
    for (j in (1 : n.deme)[-i]){
      f.mig.mat[i,j] <- effective.population[j] * migration.matrix[j,i] / effective.population[i]
    }
  }
  ###### Matrix multiplication by diag(effective.pop) better?? (left and right mult by same diag matrix)
  return(f.mig.mat)
}
