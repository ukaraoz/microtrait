#!/usr/bin/Rscript

# R implementation of consenTRAIT
# Credits: Nick Youngblut
# from https://github.com/nick-youngblut/consenTRAIT

# init
rm(list=ls())

# opt parsing
suppressPackageStartupMessages(library(docopt))

'Usage: consenTRAIT.r [options] <tree> <root> <trait>

options:
  <tree>       Newick tree file (multitree).
  <root>       Name of root taxon.
  <trait>      Trait table (no headers).
  -b=<b>       Number trait bootstraps per bootstrap tree.
               [Default: 10]
  -s=<s>       Percent shared trait cutoff.
               [Default: 90]
  -c=<c>       Cluster size file name prefix.
               [Default: cluster_size]
  -d=<d>       Cluster distance file name prefix.
               [Default: cluster_dist]
  -t=<m>       Tau_D table name.
               [Default: Tau_D.txt]
  -u=<u>       Tau_D bootstrap table name.
               [Default: Tau_D_boot.txt]
  -x=<x>       Create a test tree & trait file with `x` taxa.
               Output files: "consentrait_TEST*"
  -r=<r>       Model rate parameters for the example traits.
               Lower values mean more conservatism.
               [Default: 0.5,1000]
  -n=<b>       Number of simulated trees to make.
               NOTE: this is only for the `x` option.
               [Default: 10]
  -p=<p>       Number of processors.
               [Default: 1]
  -h           Help

description:
  The script requires two input files as arguments:
    1) a Newick Tree
    2) a tab delimited text file with names of
       each taxon in the first column and then 0 or 1
       values for each trait in the following columns
       (and no headers)

  The newick tree file should consist of multiple bootstrap
  trees, one bootstrap tree per line in the file.

  Tau_D is average consensus sequence distance (branch length)
  between trait values.

  This version of consenTRAIT will conduct the non-parametric
  bootstrapping (in parallel) for calculating p-values.

  The `x` option:
    You will still need to provide something for the `tree` `root` `trait`
    commands. Example: `consentrait.r -x 20 tmp tmp tmp`

  Dependencies: data.table, adephylo, ape, docopt, parallel
  (if using >1 processor).

  For more info: http://www.ess.uci.edu/group/amartiny/research/consentrait

  OUTPUT:
    Many accessory tables are written during the run, but the main
    table with the tauD values for each trait are written to STDOUT.
    So, to save your table, redirect the output to a file.

' -> doc
opts = docopt(doc)
opts[['-r']] = unlist(strsplit(opts[['-r']], split=','))
opts[['-r']] = sapply(opts[['-r']], as.numeric)

# packages
pkgs <- c('data.table', 'adephylo', 'ape')
for(x in pkgs){
  suppressPackageStartupMessages(library(x, character.only=TRUE))
}
if(opts[['-p']] > 1){
  suppressPackageStartupMessages(library('parallel', character.only=TRUE))
}


#-- functions --#
format.means = function(x, table){
  x = as.data.frame(t(x))
  col.n = sapply(1:(ncol(table)-1), function(x) paste(c('t', x), collapse=''))
  colnames(x) = col.n
  x$tree = 1:nrow(x)
  x = x[,c('tree', col.n)]
  return(x)
}


init_data_files = function(cluster_size_file, cluster_dist_file){
  cat(c("trait","tree","subtree", "distance","cluster_size"),
      '\n', file = cluster_size_file,
      sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
  cat(c("trait","tree","subtree", "cluster","distance"),
      '\n', file = cluster_dist_file,
      sep = "\t", fill = FALSE, labels = NULL,append = FALSE)
}

cons.trait = function(table, tree_all, opts, boot=FALSE){
  if(boot==TRUE){
    cat('Analyzing bootstrap replicate...\n', file = stderr())
  }
  # init
  n.trees = length(tree_all)
  Mean_all = matrix(nrow=ncol(table)-1,ncol=n.trees)

  for (m in 1:length(tree_all)) {
    if(boot==FALSE){
      cat("Analyzing tree: ",m,"\n", file = stderr())
    }

    # testing if table and tree contain the same entries - else drop tips
    tree = tree_all[[m]]
    z = subset(tree$tip.label,!(tree$tip.label %in% table[,1]))
    if (length(z) > 0) {
      warning("  Dropping ", length(z), " tips")
      drop.tip(tree,z)
    }

    #rooting tree with first taxon - change if different root
    message('Rooting tree on ', opts[['<root>']])
    tree = multi2di(tree)
    root_tree = root(tree,opts[['<root>']]) #,resolve.root=T)
    #replacing negative branch lengths - e.g., from PHYLIP
    message('Checking branch lengths')
    root_tree$edge.length[root_tree$edge.length <= 0] =  0.000001
    message('Getting all subtrees')
    subtree = subtrees(root_tree, wait=FALSE)

    cluster_mean = numeric(length=0)
    # loop through all traits
    message('Processing ', ncol(table)-1, 'traits')
    for (j in 2:ncol(table)) {
      if(boot==FALSE){
        message("  Analyzing trait: ",j-1)
      }
      #Loading trait table
      table_tmp = table[,c(1,j)]
      colnames(table_tmp)[1] = "ID"
      colnames(table_tmp)[2] = "Trait"
      table_tmp$Trait = as.numeric(as.character(table_tmp$Trait))

      # removing all entries not in tree
      message('Filtering trait table')
      table_tmp2 = data.table(table_tmp)
      setkey(table_tmp2,ID)
      table2 = table_tmp2[intersect(table_tmp2$ID,root_tree$tip.label)]
      setkey(table2,ID)

      #initializing result vectors and file names
      positives = vector(mode="list",length=0)
      cluster_size = numeric(length=0)
      cluster_dist = numeric(length=0)
      # output files
      cluster_size_file = paste(opts[['-c']],'_t',j-1,".txt",sep="")
      cluster_dist_file = paste(opts[['-d']],'_t',j-1,".txt",sep="")

      # Init cluster size & distance files
      if ((m == 1) & (boot == FALSE)){
        init_data_files(cluster_size_file, cluster_dist_file)
      }

      #loop through all subtrees and determining if any subtrees have >P% positives
      message('Processing all subtrees (n=', length(subtree), ')')
      for (i in 1:length(subtree)){
        tip_names = subtree[[i]]$tip.label
        # apply % shared cutoff
        if (mean(table2[tip_names][,Trait]) > perc.share.cutoff) {
          match_test = match(tip_names,positives)
          if (all(is.na(match_test))) {
            positives = c(positives,tip_names)
            cluster_dist = distRoot(subtree[[i]],tip_names, method=c("p"))
            cluster_size = append(cluster_size,mean(cluster_dist))

            # printing to files, if not a bootstrap
            if(boot == FALSE){
              cat(j-1,m,i,mean(cluster_dist),length(cluster_dist),
                  '\n', file = cluster_size_file,
                  sep = "\t", fill = FALSE, labels = NULL,append = TRUE)

              for(cdl in 1:length(cluster_dist)){
                cat(j-1,m,i,cdl,cluster_dist[cdl], '\n', file = cluster_dist_file,
                    sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
              }
            }
          }
          else if (any(is.na(match_test))) {
            cat("Assertion error: NAs present\n", file = stderr())
          }
          else {

          }
        }
      }

      ##### find singletons ######
      a = table2[table2$Trait == 1,][,ID]
      g = as.character(a)

      singletons_names = setdiff(g,positives)
      if (length(singletons_names) > 0) {
        for (h in 1:length(singletons_names)){
          # weigh singletons with half
          we = which.edge(root_tree,singletons_names[h])
          singleton_edges = 0.5*root_tree$edge.length[we]
          cluster_size = append(cluster_size,singleton_edges)

          if(boot == FALSE){
            cat(j-1,m,NA,singleton_edges,1, '\n', file=cluster_size_file,
                sep = "\t",fill = FALSE, labels = NULL,append = TRUE)
          }
        }
      }
    # means of cluster sizes
    Mean_all[j-1,m] = mean(cluster_size)
    }
  }
  return(Mean_all)
}

sim_trait = function(phy, model, tries=1000){
  # simulating a binary trait.
  # phy = phylogeny, model = rate_model, tries=number of tries to have trait evolve at least once
  try_cnt = 1
  while(1==1){
    trait = rTraitDisc(phy, model=model, states=c(0, 1), root.value=1)
    trait = ifelse(trait == 1, 0, 1)
    if(sum(trait)> 0){
      return(trait)
    }
    try_cnt = try_cnt + 1
    if(try_cnt > tries){
      stop('Exceeded tries to evolve the trait. Raise the rate parameter values')
    }
  }
}



#-- main --#
# Test files
if(!is.null(opts[['-x']])){
  out.tree = 'consentrait_TEST.nwk'
  out.trait = 'consentrait_TEST.txt'
  out.tree.pic1 = 'consentrait_TEST_trait1.pdf'
  out.tree.pic2 = 'consentrait_TEST_trait2.pdf'

  # tree(s)
  n.taxa = as.numeric(opts['-x'])
  tree = rcoal(n.taxa)
  write.tree(tree, out.tree)

  ## traits
  taxa = tree$tip.label
  model1 = matrix(c(0,0,1,0), 2) * opts[['-r']][1]
  trait1 = sim_trait(tree, model1)
  model2 = matrix(c(1,1,1,1), 2) * opts[['-r']][2]
  trait2 = sim_trait(tree, model2)

  df = data.frame('taxa' = taxa,
    'trait1' = trait1,
    'trait2' = trait2)
  write.table(df, out.trait, sep='\t', quote=FALSE,
              row.names=FALSE, col.names=FALSE)

  # plotting tree
  pdf(out.tree.pic1, bg='white')
  plot(tree, show.tip.label=FALSE, adj=1)
  trait1 = ifelse(trait1==0, 'blue', 'red')
  tiplabels(pch=22, col=NULL, bg=trait1)
  dev.off()
  pdf(out.tree.pic2, bg='white')
  plot(tree, show.tip.label=FALSE, adj=1)
  trait2 = ifelse(trait2==0, 'blue', 'red')
  tiplabels(pch=22, col=NULL, bg=trait2)
  dev.off()

  # status
  msg = paste(c('Test files written: ', out.tree,
    out.trait, out.tree.pic1, out.tree.pic2), collapse='\n  ')
  cat(paste0(msg,'\n'), file = stderr())

  opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()
}


#-- consenTRAIT run --#
# Params
perc.share.cutoff = as.numeric(opts[['-s']])
stopifnot((perc.share.cutoff >= 0) &  (perc.share.cutoff <= 100))
perc.share.cutoff = perc.share.cutoff / 100


# Import
## Newick tree (multitree) - replace to read.nexus if using nexus tress
tree_all = read.tree(opts[['<tree>']],keep.multi = TRUE)
## Trait table w. no headers
table = read.table(opts[['<trait>']], sep = "\t", header=FALSE)

# Tau_D for each bootstrap tree
Mean_all = cons.trait(table, tree_all, opts)
Mean_all = format.means(Mean_all, table)
## writting
write.table(Mean_all,opts[['-t']], sep = "\t", quote=FALSE, row.names=FALSE)


# non-paramtric bootstrapping
## making randomly arranged trait tables
#random.traits = function(df){
table.l = list()
for(i in 1:opts[['-b']]){
  df.rand = apply(as.data.frame(table[,2:ncol(table)]), 2,
    function(x) sample(c(0,1), length(x), replace=TRUE))
  df.rand = as.data.frame(df.rand)
  tmp = colnames(df.rand)
  df.rand$V1 = table[,1]
  df.rand = df.rand[,c('V1', tmp)]
  table.l[[i]] = df.rand
}

## (parallel) calling of cons.trait
if(opts[['-p']] > 1){
  cat('Bootstrapping in parallel\n', file = stderr())
  cl1 = makeCluster(opts[['-p']], type='FORK')
  mean_boots = parLapply(cl1, table.l, cons.trait,
    tree_all=tree_all, opts=opts, boot=TRUE)
  stopCluster(cl1)
} else {
  mean_boots = lapply(table.l, cons.trait,
    tree_all=tree_all, opts=opts, boot=TRUE)
}
## formatting output
for(i in 1:length(mean_boots)){
  tmp = format.means(mean_boots[[i]], table)
  tmp$bootstrap = i
  mean_boots[[i]] = tmp
}
mean_boots = do.call(rbind, mean_boots)
write.table(mean_boots,opts[['-u']], sep = "\t", quote=FALSE, row.names=FALSE)


# calcalating p-value
## getting mean Tau_D of real data
mean_tauD = apply(as.data.frame(Mean_all[,2:ncol(Mean_all)]), 2, mean)
if(length(mean_tauD) == 1){
  names(mean_tauD) = c('t1')
}
## determine p-value
cat('Trait\ttau_D\tp-value', '\n')
for(n in names(mean_tauD)){
  tau_D = mean_tauD[n]
  boot_tauD = mean_boots[,n]
  p = 1 - sum(tau_D < boot_tauD) / length(boot_tauD)
  line = paste(c(n, tau_D, p), collapse='\t')
  cat(line, '\n')
}
