dta.likelihood <- function(ED, effective_population, gen_length, migration_matrix, node_indices){
  root_row <- which(is.na(ED[,2]))

  f_migration_matrix <- forward.migration.matrix(migration_matrix, effective_population)
  f_migration_matrix_rowsums <- rowSums(f_migration_matrix)
  log_likelihood <- 0

  for (i in (1 : dim(ED)[1])[-root_row]){
    node_parent_row <- node_indices[ED[i,2]]
    time_increment <- ED[i, 6] - ED[node_parent_row, 6]

    parent_deme <- ED[node_parent_row, 5]
    node_deme <- ED[i, 5]

    if (parent_deme == node_deme){
      log_likelihood <- log_likelihood - f_migration_matrix_rowsums[parent_deme] * time_increment
    } else{
    log_likelihood <- log_likelihood - f_migration_matrix_rowsums[parent_deme] * time_increment + log(f_migration_matrix[parent_deme, node_deme])
    }
  }
  return(list(log.likelihood = log_likelihood, likelihood = exp(log_likelihood)))
}
