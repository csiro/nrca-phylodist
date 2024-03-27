# to be honest, not sure what all the libraries were needed for
# some may have only been necessary for additional bells and whistles beyond what we need now

library(ape)
library(phytools)
library(adephylo)
library(MASS)

# function to recursively collect descendants of node; is used by degreesofsep function

descendantlist <- function(thistree, thisnode)
{
  if (thisnode <= length(thistree$tip.label))
  {
    return (thisnode)
  }
  else
  {
    wherenext <- which(thistree$edge[,1]==thisnode)  # get immediate descendants
    thislist <- NULL
    for (x in 1:length(wherenext))
    {
      thislist <- c(thislist, descendantlist(thistree, thistree$edge[wherenext[x],2]))
    }
    return(thislist)
  }
}

# function to calculate degrees of phylogenetic separation between tree tips
# seems to work, but results are, of course and in contrast to phylogenetic distances, asymmetric
# read resulting matrix line-wise: line shows degrees of separation of taxa in columns, from
# perspective of the taxon in the line

degreesofsep <- function(thistree)
{
  dosmatrix <- matrix(0, nrow=length(thistree$tip.label), ncol=length(thistree$tip.label))
  colnames(dosmatrix) <- thistree$tip.label
  rownames(dosmatrix) <- thistree$tip.label
  for (x in 1:length(thistree$tip.label))
  {
    prior_y <- x         # start at present terminal
    y <- thistree$edge[which(thistree$edge[,2]==x),1]   # get immediately ancestral node
    currentdist <- 0
    while (y != (length(thistree$tip.label)+1)) # move downtree until root node is found
    {
      currentdesc <- which(thistree$edge[,1]==y)
      for (z in 1:length(currentdesc))
      {
        if (thistree$edge[currentdesc[z],2]!=prior_y)
        {
          dosmatrix[x,descendantlist(thistree,thistree$edge[currentdesc[z],2])] <- currentdist
        }
      }
      prior_y <- y
      y <- thistree$edge[which(thistree$edge[,2]==y),1]   # get immediately ancestral node
      currentdist <- currentdist + 1
    }
    currentdesc <- which(thistree$edge[,1]==y)
    for (z in 1:length(currentdesc))
    {
      if (thistree$edge[currentdesc[z],2]!=prior_y)
      {
        dosmatrix[x,descendantlist(thistree,thistree$edge[currentdesc[z],2])] <- currentdist
      }
    }
  }
  return(dosmatrix)
}

# read phylogenetic tree
mytree <- read.tree("astereae_concatenated.tre")

# reroot, if necessary, using two terminals to define clade
# not sure how to do that for end-user
# one can, of course, define outgroups before phylogenetic analysis, if one were so inclined

myOG <- getMRCA(mytree, c("Dimorphotheca_pluvialis", "Ewartia_nubigena", "Abrotanella_nivigena","Cotula_coronopifolia"))
mytree <- root(mytree, node = myOG)

# infer matrix of pairwise patristic distances between all terminals
# this takes quite some time for the example tree, because it is very large

myPatristic <- distTips(mytree)
myPatristicM <- as.matrix(myPatristic)
myPatristicMordered <- myPatristicM[order(rownames(myPatristicM)), order(rownames(myPatristicM))]
write.matrix(myPatristicMordered, file="patristicdists.tsv", sep="\t")

# now calculate degrees of separation, i.e. counting nodes between any terminal and its ancestral lineage splits
# this takes quite some time for the example tree, because it is very large

myDegsep <- degreesofsep(mytree)
myDegsep <- myDegsep[order(rownames(myDegsep)),order(rownames(myDegsep))]
write.matrix(myDegsep, file="degsep.tsv", sep="\t")

# make data frame for one target species with both its degrees of separation and patristic distances
# this would probably make more sense as a function that gets myspecies and the data as parameters

myspecies <- "Erigeron_bonariensis"
Terminal <- row.names(myDegsep)
Degsep <- myDegsep[which(row.names(myDegsep)==myspecies),]
PatristicDist <- myPatristicMordered[which(row.names(myPatristicMordered)==myspecies),]
PhyloDists <- data.frame(Terminal[order(PatristicDist)],Degsep[order(PatristicDist)],PatristicDist[order(PatristicDist)])
write.table(PhyloDists, file="phylodists_Erigeron_bonariensis.tsv", sep="\t")
