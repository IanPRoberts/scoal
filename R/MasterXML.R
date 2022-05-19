#' Generates xml for structured coalescent simulation using MASTER
#'
#' Generates a .xml file to be input to BEAST 2 for structured coalescent simulations using MASTER
#'
#' @param effective.pop effective population size from which the sample is taken
#' @param migration.matrix matrix of migration rates between demes
#' @param data  nx3 matrix with first column giving the tip labels, second column the time at which the sample was taken and third column the initial deme of each sample point
#' @param save.path file path to output .xml file to
#' @param run.xml if TRUE, run generated .xml through BEAST
#'
#' @return
#'
#' @export

master.xml <- function(effective.pop, migration.matrix, leaf.data, save.path = ".", run.xml = FALSE){
  if (save.path == "."){
    save.path <- getwd()
  }

  n.deme <- length(effective.pop)
  out <- paste0("<beast version='2.0' namespace='master:master.model:master.steppers:master.conditions:master.postprocessors:master.outputs'> \n",
  "\t <run spec='InheritanceTrajectory' \n",
  "\t\t verbosity='2'> \n\n",
  "\t\t <model spec='Model'> \n",
  "\t\t\t <populationType spec='PopulationType' typeName='L' id='L' dim='", n.deme, "'/> \n",
  "\t\t\t <reactionGroup spec='ReactionGroup' reactionGroupName='Coalescence'> \n")

  for (i in 1 : n.deme){
    out <- paste0(out,
                       "\t\t\t\t <reaction spec='Reaction' rate='", 1/effective.pop[i], "'> \n",
                       "\t\t\t\t\t 2L[", i-1, "]:1 -> L[", i-1, "]:1 \n",
                       "\t\t\t\t </reaction> \n")
  }

  out <- paste0(out,
                     "\t\t\t </reactionGroup> \n\n",
                     "\t\t\t <reactionGroup spec='ReactionGroup' reactionGroupName='Migration'> \n")

  for (i in 1 : n.deme){
    for (j in (1:n.deme)[-i]){
      out <- paste0(out,
                         "\t\t\t\t <reaction spec='Reaction' rate='", migration.matrix[i,j], "'> \n",
                         "\t\t\t\t\t L[", i-1, "] -> L[", j-1, "] \n",
                         "\t\t\t\t </reaction> \n")
    }
  }

  out <- paste0(out,
                     "\t\t\t </reactionGroup> \n\n",
                     "\t\t </model> \n\n",
                     "\t\t <initialState spec='InitState'> \n")

  leaf.times <- unique(data[,2])
  for (time in leaf.times){
    for (deme in 1 : n.deme){
      leaf.set <- (1 : dim(data)[1])[(data[,2] == time) & (data[,3] == deme)]
      if (length(leaf.set) > 0){
        out <- paste0(out,
                      "\t\t\t <lineageSeedMultiple spec='MultipleIndividuals' copies='", length(leaf.set), "' time='", time, "'> \n",
                      "\t\t\t\t <population spec='Population' type='@L' location='", deme - 1, "'/> \n",
                      "\t\t\t </lineageSeedMultiple> \n")
      }
    }
  }

  out <- paste0(out,
                "\t\t </initialState> \n\n",
                "\t\t <lineageEndCondition spec='LineageEndCondition' nLineages='1'/> \n\n",
                "\t\t <output spec='NexusOutput' fileName='", save.path, "/SCMasterSimTree_out.nexus' reverseTime='true'/> \n",
                  "\t </run> \n",
                  "</beast>")

  file.create(paste0(save.path, "/SCMasterSim.xml"))
  writeLines(out, con = paste0(save.path, "/SCMasterSim.xml"))

  if (run.xml == TRUE){
    system(paste0('java -jar "C:/Program Files/BEAST/lib/launcher.jar" "', save.path, '/SCMasterSim.xml"'))
  }
}
