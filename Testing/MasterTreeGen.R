#### Read BEAST nexus
n.deme <- 2
n <- 150
data <- matrix(c(1:n, rep_len(2022, n), rep_len(1:n.deme, n)), n, 3)
mig_mat <- matrix(0.01, n.deme, n.deme); diag(mig_mat) <- 0 #matrix(c(0, 0.01, 0.02, 0.03, 0, 0.04, 0.05, 0.06, 0), 3, 3)
eff_pop <- rep(120, n.deme) #10 * 1 : n.deme
scoal::master.xml(effective.pop = eff_pop,
                  migration.matrix = mig_mat,
                  leaf.data = data,
                  save.path = ".",
                  run.xml = TRUE)

file.path <- "./SCMasterSimTree_out.nexus"

treedata <- treeio::read.beast(file.path)
node.deme <- 1 + as.numeric(treedata@data$location[order(as.numeric(treedata@data$node))])

phylo <- treedata@phylo
phylo$node.deme <- node.deme
structured.plot(phylo)
# saveRDS(phylo.to.ed(phylo), "./SCMasterED.RDS")
