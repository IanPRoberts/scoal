library(ape)

test_that("Migration pair birth adds two migrations",{
  data <- matrix(c(1:10, rep(2022, 10), rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  homo.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(homo.phylo)
  node.indices <- 1 : dim(ED)[1]
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.pair.birth(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_equal(M+2, M.2)
    M <- M.2
  }

  data <- matrix(c(1:10, rep(2022, 4), rep(2020, 5), 2019, rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  hetero.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(hetero.phylo)
  node.indices <- 1 : dim(ED)[1]
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.pair.birth(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_equal(M+2, M.2)
    M <- M.2
  }
})

test_that("Migration pair death remove zero or two migrations",{
  data <- matrix(c(1:10, rep(2022, 10), rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  homo.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(homo.phylo)
  node.indices <- 1 : dim(ED)[1]
  for(i in 1 : 20){
    proposal <- ed.mig.pair.birth(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
  }
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for (j in 1 : 5){
    proposal <- ed.mig.pair.death(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M-2))
    M <- M.2
  }

  data <- matrix(c(1:10, rep(2022, 4), rep(2020, 5), 2019, rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  hetero.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(hetero.phylo)
  node.indices <- 1 : dim(ED)[1]
  for(i in 1 : 20){
    proposal <- ed.mig.pair.birth(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
  }
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.pair.death(ED, 3, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M-2))
    M <- M.2
  }
})

test_that("Migration birth adds zero or one migration",{
  data <- matrix(c(1:10, rep(2022, 10), rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  homo.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(homo.phylo)
  node.indices <- 1 : dim(ED)[1]
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.birth(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M+1))
    M <- M.2
  }

  data <- matrix(c(1:10, rep(2022, 4), rep(2020, 5), 2019, rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  hetero.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(hetero.phylo)
  node.indices <- 1 : dim(ED)[1]
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.birth(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M+1))
    M <- M.2
  }
})

test_that("Migration death removes zero or one migration",{
  data <- matrix(c(1:10, rep(2022, 10), rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  homo.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(homo.phylo)
  node.indices <- 1 : dim(ED)[1]
  for (i in 1 : 20){
    proposal <- ed.mig.birth(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
  }
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.death(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M-1))
    M <- M.2
  }

  data <- matrix(c(1:10, rep(2022, 4), rep(2020, 5), 2019, rep(1, 4), rep(2,4), rep(3, 2)),nrow = 10, ncol = 3)
  hetero.phylo <- Structured.sim(data, 1, 1, 3, matrix(c(0, rep(1/10, 3), 0, rep(1/10, 3), 0),3,3), FALSE)
  ED <- phylo.to.ed(hetero.phylo)
  for (i in 1 : 20){
    proposal <- ed.mig.birth(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
  }
  M <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
  for(i in 1 : 5){
    proposal <- ed.mig.death(ED, 3, TRUE, node.indices)
    ED <- proposal$ED
    node.indices <- proposal$node.indices
    M.2 <- sum((is.na(ED[,4])) & (!is.na(ED[,3])))
    expect_true(M.2 %in% c(M, M-1))
    M <- M.2
  }
})
