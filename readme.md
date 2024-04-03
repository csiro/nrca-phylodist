---
Calculating phylogenetic distances for weed biological control

Authors of code: Nunzio Knerr, Stephanie Chen, Alexander Schmidt-Lebuhn
---

## What does this script do and what is it useful for?
This code was introduced in a submitted paper titled 'Phylogenomics-driven host test list selection for weed biological control'. It contains functions for calculating phylogenetic distance measures useful for creating host tests list in classical weed biological control.

The degree of relatedness between two taxa on a phylogeny is indicated by the number of nodes separating them. Here, we provide functions that calculate two distance measures, degree of separation i.e. node count and patristic distance, given an input phylogenetic tree.

## Descendant List Function
First define a function to recursively collect descendants of a node. This is used by the 'degreeofsep' function later on.

``` {.r .cell-code}
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
```

## Degree Of Separation Function
Function for calculating degrees of separation from a specified target weed.

``` {.r .cell-code}
# function to calculate degrees of phylogenetic separation between tree tips
# results are in contrast to phylogenetic distances, asymmetric
# read resulting matrix line-wise: line shows degrees of separation of taxa in columns, from the target taxon

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
```

## User Input Variables
Specify the inputs and outputs for use in the script. A tree file in newick format is required. The outgroup(s) may be specified. The target taxon i.e. target weed for biological control is also specified here so that the distance measures can be calculated in relation to the target.

``` {.r .cell-code}
# the tree to use
treeFileName <- "astereae_concatenated.tre"
# a list of taxa in the outgroup
taxonListForOutgroup <- c("Dimorphotheca_pluvialis", "Ewartia_nubigena", "Abrotanella_nivigena","Cotula_coronopifolia") 
# the target taxon to calculate the distances/separation from
myTargetTaxon <- "Erigeron_bonariensis"
# the output file name to use
outputFileName <- "phylodists_Erigeron_bonariensis.tsv"
```

## Example Usage

``` {.r .cell-code}
# load libraries
library(ape)
library(adephylo)
```

``` {.r .cell-code}
# read phylogenetic tree
mytree <- ape::read.tree(treeFileName)

# call the get Most Recent Common Ancesstors (MRCA)
myOG <- getMRCA(mytree, taxonListForOutgroup)
#root the tree based on the MRCA results
mytree <- root(mytree, node = myOG)

# infer matrix of pairwise patristic distances between all terminals
# this takes quite some time for larger trees

myPatristic <- distTips(mytree)
myPatristicM <- as.matrix(myPatristic)
myPatristicMordered <- myPatristicM[order(rownames(myPatristicM)), order(rownames(myPatristicM))]
#write.matrix(myPatristicMordered, file="patristicdists.tsv", sep="\t")

# now calculate degrees of separation, i.e. counting nodes between any terminal and its ancestral lineage splits
# this will take quite some time for larger trees

myDegsep <- degreesofsep(mytree)
myDegsep <- myDegsep[order(rownames(myDegsep)),order(rownames(myDegsep))]
#write.matrix(myDegsep, file="degsep.tsv", sep="\t")

# make data frame for one target species with both its degrees of separation and patristic distances
targetTaxon <- myTargetTaxon

Terminal <- row.names(myDegsep)

Degsep <- myDegsep[which(row.names(myDegsep)==targetTaxon),]

PatristicDist <- myPatristicMordered[which(row.names(myPatristicMordered) ==
                                            targetTaxon),]

PhyloDists <- data.frame(Terminal[order(PatristicDist)], 
                         Degsep[order(PatristicDist)], 
                         PatristicDist[order(PatristicDist)])

write.table(PhyloDists, file=outputFileName, sep="\t")

knitr::kable(PhyloDists)
```

The output is a tab separated file with columns for the scientific name, degree of separation and patristic distance. See the file phylodists_Erigeron_bonariensis.tsv for an example that was used as a case study in the published paper.

## Citation and contact information
Submitted paper: 'Phylogenomics-driven host test list selection for weed biological control'. The corresponding author of the paper is Alexander Schmidt-Lebuhn (alexander.s-l@csiro.au).
