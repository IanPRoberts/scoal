#' Generates xml for structured coalescent simulation using MASTER
#'
#' Generates a .xml file to be input to BEAST 2 for structured coalescent simulations using MASTER
#'
#' @param effective.pop effective population size from which the sample is taken
#' @param migration.matrix matrix of migration rates between demes
#' @param leaf.data  nx3 matrix with first column giving the tip labels, second column the time at which the sample was taken and third column the initial deme of each sample point
#' @param save.path file path to output .xml file to
#' @param run.xml if TRUE, run generated .xml through BEAST
#'
#' @return xml file
#'
#' @export

master.xml <- function(effective.pop, migration.matrix, leaf.data, save.path = ".", run.xml = FALSE){
  if (save.path == "."){
    save.path <- getwd()
  }

  n.deme <- length(effective.pop)
  out <- paste0("<beast version='2.0' namespace='master:master.model:master.steppers:master.conditions:master.postprocessors:master.outputs'> \n",
  "\t <run spec='InheritanceTrajectory' \n",
  "\t\t verbosity='0'> \n\n",
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

  leaf.times <- unique(leaf.data[,2])
  for (time in leaf.times){
    for (deme in 1 : n.deme){
      leaf.set <- (1 : dim(leaf.data)[1])[(leaf.data[,2] == time) & (leaf.data[,3] == deme)]
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


#' Generates BEAST2 xml for MultiType Tree
#'
#' Generates a .xml file to run BEAST2 using either MultiTypeTree or BASTA with a fixed tree.
#' Initialises with all coalescent rates the same and all backward-in-time migration rates the same
#'
#' @param strphylo structured phylo object giving initialisation condition for fixed tree run with MultiType Tree
#' @param coal_rate initial estimate for coalescent rates
#' @param bit_mig_rate initial estimate for backward-in-time migration matrix
#' @param N total number of MCMC iterations (including burn-in)
#' @param thin thinning increment for the MCMC
#' @param con A connection object or a character string
#' @param BEAST2_package Select package to prepare xml file for, either MTT (MultiTypeTree) or BASTA
#' @param run_name Name for logger files to be saved as
#'
#' @return output file or file content on screen
#'
#' @export

fixed_tree_xml <- function(strphylo, n_deme, coal_rate, bit_mig_mat, N=1e7, thin=1e3, con = stdout(), BEAST2_package = "MTT", run_name = "scoal"){
  n_tip <- length(strphylo$tip.label)
  node_ages <- ape::node.depth.edgelength(strphylo)

  # Format integer initial values as integer.0 instead of just integer
  # if (isTRUE(bit_mig_rate == floor(bit_mig_rate))){
  #   bit_mig_rate <- sprintf("%.1f", bit_mig_rate)
  # } else {
  #   bit_mig_rate <- sprintf("%f", bit_mig_rate)
  # }
  # if (isTRUE(coal_rate == floor(coal_rate))){
  #   coal_rate <- sprintf("%.1f", coal_rate)
  # } else {
  #   coal_rate <- sprintf("%f", coal_rate)
  # }

  out <- paste("<beast version='2.7'",
	"namespace='beast.base.evolution.alignment:",
	"beast.pkgmgmt:",
	"beast.base.core:",
	"beast.base.inference:",
	"beast.pkgmgmt:",
	"beast.base.core:",
	"beast.base.inference.parameter:",
	"beast.base.evolution.tree:",
	"beast.base.evolution.tree.coalescent:",
	"beast.base.util:",
	"beast.base.math:",
	"beast.base.evolution.operator:",
	"beast.base.inference.operator:",
	"beast.base.evolution.sitemodel:",
	"beast.base.evolution.substitutionmodel:",
	"beast.base.evolution.likelihood:",
	"beast.evolution.migrationmodel:",
	"beast.base.inference.distribution:",
	"multitypetree.distributions:",
	"multitypetree.operators:",
	"multitypetree.util'>\n",
	"<!-- Alignment -->",
	"<alignment spec=\"beast.base.evolution.alignment.Alignment\" id=\"alignment\" dataType=\"nucleotide\">", sep ="\n\t")

  out <- paste(out,
               paste0("<sequence taxon='", strphylo$tip.label, "' value='?'/>", collapse = "\n\t\t"), sep = "\n\t\t")

  out <- paste(out,
                "</alignment> \n",
  ###### Leaf demes
                "<typeTraitSet",
                "\t spec='TraitSet'",
                "\t id='typeTraitSet'",
                "\t traitname='type'",
                "\t value='", sep = "\n\t")

  ###### Tip demes
  out <- paste(out,
               paste0(strphylo$tip.label, "=", strphylo$node.deme[1:n_tip], collapse = ",\n\t\t\t"), sep = "\n\t\t\t")

  out <- paste(out,
               "'>",
               "\t <taxa spec='TaxonSet' alignment='@alignment'/>",
               " </typeTraitSet> \n",
  ##### Leaf times
               " <timeTraitSet",
               "\t spec='TraitSet'",
               "\t id='timeTraitSet'",
               "\t traitname='date-backward'",
                "\t value='\n", sep = "\n\t")

  out <- paste(out,
               paste0(strphylo$tip.label, "=", node_ages[1:n_tip], collapse = ",\n\t\t\t"), sep = "\n\t\t\t")

  out <- paste(out,
               "'>",
               "\t <taxa spec='TaxonSet' alignment='@alignment'/>",
               " </timeTraitSet> \n", sep = "\n\t")

  ##### HKY Substitution Model
  # Fixed tree method, substitution model should be fairly irrelevant??
  out <- paste(out,
                " <!-- Substitution model (HKY) -->",
                " <siteModel spec=\"SiteModel\" id=\"siteModel\">",
                "\t <mutationRate spec='RealParameter' id=\"mutationRate\" value=\"1.0\"/>",
                "\t <substModel spec=\"HKY\">",
                "\t\t <kappa spec='RealParameter' id=\"hky.kappa\" value=\"1.0\"/>",
                "\t\t <frequencies estimate=\"false\" spec='Frequencies'>",
                "\t\t\t <frequencies spec='RealParameter' id=\"hky.freq\" value=\"0.25 0.25 0.25 0.25\"/>",
                "\t\t </frequencies>",
                "\t </substModel>",
                " </siteModel> \n", sep = "\n\t")

  ##### Set up migration model
  out <- paste(out,
                " <!-- Migration model -->",
                " <migrationModel spec='multitypetree.evolution.tree.SCMigrationModel' id='migModel'>",
                paste0("\t <rateMatrix spec='RealParameter' dimension='", n_deme * (n_deme - 1), "' id=\"rateMatrix\">"),
                paste0("\t\t ", paste(bit_mig_mat[-(1 + 0:(n_deme - 1) * (n_deme + 1))], collapse = " ")),
                "\t </rateMatrix>",
                paste0("\t <popSizes spec='RealParameter' dimension=\"", n_deme, "\" id=\"popSizes\">"),
                paste0("\t\t ", paste(1/coal_rate, collapse = " ")),
                "\t </popSizes>",
               "\t <typeSet id=\"typeSet\" spec='multitypetree.evolution.tree.TypeSet' typeTraitSet=\"@typeTraitSet\"/>",
                " </migrationModel> \n\n", sep = "\n\t")

  ##### Prior setup
  out <- paste(out,
                "<!-- Parameter priors -->",
                "<input spec='CompoundDistribution' id='parameterPriors'>",
                # Mutation rate
                "\t <distribution spec='beast.base.inference.distribution.Prior' x=\"@mutationRate\">",
                "\t\t <distr spec='LogNormalDistributionModel' M=\"0.0\" S=\"4.0\"/>",
                "\t </distribution> \n",
                # HKY Kappa
                "\t <distribution spec='beast.base.inference.distribution.Prior' x=\"@hky.kappa\">",
                "\t\t <distr spec='LogNormalDistributionModel' M=\"0.0\" S=\"4.0\"/>",
                "\t </distribution> \n",
                # Migration rates
                "\t <distribution spec='beast.base.inference.distribution.Prior' x=\"@rateMatrix\">",
                "\t\t <distr spec='LogNormalDistributionModel' M=\"0.0\" S=\"4.0\"/>",
                "\t </distribution> \n",
                # Coalescent rates
                "\t <distribution spec='beast.base.inference.distribution.Prior' x=\"@popSizes\">",
                "\t\t <distr spec=\"LogNormalDistributionModel\"  M=\"0.0\" S=\"4.0\"/>",
                "\t </distribution>",
                "</input> \n", sep = "\n\t")

  ##### Structured Coalescent Density
  out <- paste(out,
                "<!-- Probability of tree given migration rates and population sizes -->",
                "<input spec='StructuredCoalescentTreeDensity' id='treePrior'>",
                "\t <multiTypeTree idref=\"tree\"/>",
                "\t <migrationModel idref=\"migModel\"/>",
                "</input> \n",
               ##### MCMC setup
               paste0("<run spec=\"MCMC\" id=\"mcmc\" chainLength=\"", sprintf("%i", N), "\" storeEvery=\"", sprintf("%i", thin) ,"\">"),
               sep = "\n\t")


  # Convert strphylo to BEAST2 metacommented Newick
  phylo <- strphylo
  class(phylo) <- "phylo"
  phylo$node.deme <- phylo$log.likelihood <- phylo$likelihood <- NULL
  treedata <- tidytree::as.treedata(phylo)
  treedata@data <- tidytree::tibble(type = paste0("\"", strphylo$node.deme, "\""), node = 1:length(strphylo$node.deme))

  newick_tree <- treeio::write.beast.newick(treedata)

  out <- paste(out,
               "<!-- Initialize tree from Newick string -->",
               "<init spec='multitypetree.evolution.tree.MultiTypeTreeFromNewick'",
               "\t id='tree' >",
               paste0("\t <![CDATA[ \n\t\t", newick_tree, "\n\t ]]>"),
               "\t <typeSet idref='typeSet'/>",
               "\t <typeTrait idref='typeTraitSet'/>",
               "\t <trait idref='timeTraitSet'/>",
               "</init> \n",
               ##### Initial condition
               "<state> ",
               "\t <stateNode idref=\"tree\"/>",
               "\t <stateNode idref=\"rateMatrix\"/>",
               "\t <stateNode idref=\"popSizes\"/>",
               "\t <stateNode idref=\"mutationRate\"/>",
               "\t <stateNode idref=\"hky.kappa\"/>",
               "\t <stateNode idref=\"hky.freq\"/>",
               "</state>\n",
               ##### Posterior distribution
               "<distribution spec='CompoundDistribution' id='posterior'>",
               "\t <distribution idref='treePrior'/>",
               "\t <distribution idref=\"parameterPriors\"/>",
               "</distribution> \n",
               sep = "\n\t")


  if (BEAST2_package == "MTT"){
    ##### MTT operators
    out <- paste(out,
                 "<!-- parameter scaling operators -->",
                 # Migration rates scaler
                 "<operator spec='ScaleOperator'",
                 "\t id='RateScaler'",
                 "\t parameter=\"@rateMatrix\"",
                 "\t scaleFactor=\"0.8\"",
                 "\t weight=\"1\"/>",
                 #Population sizes scaler
                 "<operator spec='ScaleOperator'",
                 "\t id='PopSizeScaler'",
                 "\t parameter=\"@popSizes\"",
                 "\t scaleFactor=\"0.8\"",
                 "\t weight=\"1\"/>",
                 #HKY model parameters (not updated, commented out)
                 "<!--",
                 "<operator spec='ScaleOperator'",
                 "\t id='muRateScaler'",
                 "\t parameter=\"@mutationRate\"",
                 "\t scaleFactor=\"0.8\"",
                 "\t weight=\"1\"/>",
                 ##
                 "<operator spec='ScaleOperator'",
                 "\t id='kappaScaler'",
                 "\t parameter=\"@hky.kappa\"",
                 "\t scaleFactor=\"0.8\"",
                 "\t weight=\"0.1\"/>",
                 ##
                 "<operator spec='DeltaExchangeOperator'",
                 "\t id='freqExchanger'",
                 "\t parameter=\"@hky.freq\"",
                 "\t delta=\"0.01\"",
                 "\t weight=\"0.1\"/>",
                 "-->",
                 # Multi-type tree operators
                 "<!-- Multi-type tree operators -->",
                 "<!--",
                 ## Subtree exchange
                 "<operator",
                 "\t spec='TypedSubtreeExchange'",
                 "\t id='STX'",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t migrationModel=\"@migModel\"/>",
                 #Wilson-Balding
                 "<operator",
                 "\t spec=\"TypedWilsonBalding\"",
                 "\t id=\"TWB\"",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t migrationModel=\"@migModel\"",
                 "\t alpha=\"0.2\"/>",
                 "-->",
                 # Node retype
                 "<operator",
                 "\t spec=\"NodeRetype\"",
                 "\t id=\"NR\"",
                 "\t weight=\"1\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t migrationModel=\"@migModel\"/>",
                 # Node shift and retype 1 (root)
                 "<!--",
                 "<operator",
                 "\t spec=\"NodeShiftRetype\"",
                 "\t id=\"NSR1\"",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t rootScaleFactor=\"0.8\"",
                 "\t migrationModel=\"@migModel\"",
                 "\t rootOnly=\"true\"/>",
                 # Node shift & retype 2 (non-root)
                 "<operator",
                 "\t spec=\"NodeShiftRetype\"",
                 "\t id=\"NSR2\"",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t migrationModel=\"@migModel\"",
                 "noRoot=\"true\"/>",
                 # Multi-type uniform
                 "<operator",
                 "\t spec=\"MultiTypeUniform\"",
                 "\t id=\"MTU\"",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t includeRoot=\"true\"",
                 "\t rootScaleFactor=\"0.9\"/>",
                 # Multi-type tree scale
                 "<operator",
                 "\t spec=\"MultiTypeTreeScale\"",
                 "\t id=\"MTTS1\"",
                 "\t weight=\"10\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t scaleFactor=\"0.98\"",
                 "\t useOldTreeScaler=\"true\">",
                 "\t <parameter idref=\"popSizes\"/>",
                 "\t <parameterInverse idref=\"rateMatrix\"/>",
                 "\t <parameterInverse idref=\"mutationRate\"/>",
                 "</operator>",
                 #Multi-type tree scale 2
                 "<operator",
                 "\t spec=\"MultiTypeTreeScale\"",
                 "\t id=\"MTTS2\"",
                 "\t weight=\"1\"",
                 "\t multiTypeTree=\"@tree\"",
                 "\t migrationModel=\"@migModel\"",
                 "\t scaleFactor=\"0.98\"",
                 "\t useOldTreeScaler=\"true\">",
                 "</operator>",
                 "-->"
                 ,sep = "\n\t")
  } else if (BEAST2_package == "BASTA"){
    ##### BASTA operators
  } else {
    stop("Invalid BEAST2_package input; selected from MTT or BASTA")
  }

  out <- paste(out,
               "<!-- Loggers -->",
               # Tracer log file
               "<logger",
               paste0("\t logEvery=\"", thin, "\""),
               "\t fileName=\"$(filebase).log\">",
               "\t <model idref='posterior'/>",
               "\t <log idref=\"posterior\"/>",
               "\t <log idref=\"treePrior\"/>",
               "\t <log spec='TreeHeightLogger' tree='@tree'/>",
               "<log id=\"migModelLogger\" spec=\"MigrationModelLogger\" migrationModel=\"@migModel\" multiTypeTree=\"@tree\"/>",
               "\t <!--",
               "\t <log idref=\"mutationRate\"/>",
               "\t <log idref=\"hky.kappa\"/>",
               "\t <log idref=\"hky.freq\"/>",
               "\t -->",
               "</logger>",
               # Trees output
               "<logger",
               paste0("\t logEvery=\"", thin, "\""),
               "\t fileName=\"$(filebase).$(tree).trees\"",
               "\t mode=\"tree\">",
               "\t <log idref=\"tree\"/>",
               "</logger>",
               #
               "<logger logEvery=\"1000\">",
               "\t <model idref='posterior'/>",
               "\t <log idref=\"posterior\"/>",
               "\t<log idref=\"treePrior\"/>",
               "\t <log spec='TreeHeightLogger' tree='@tree'/>",
               "\t <ESS spec='beast.base.inference.util.ESS' name='log' arg=\"@treePrior\"/>",
               "\t <ESS spec='beast.base.inference.util.ESS' name='log' arg=\"@rateMatrix\"/>",
               "\t <ESS spec='beast.base.inference.util.ESS' name='log' arg=\"@popSizes\"/>",
               "</logger>",
               "</run>",
               sep = "\n\t")

  out <- paste0(out, "\n </beast>")

  writeLines(out, con)
}
