# Calculating phylogenetic distances for weed biological control

Authors of code: Nunzio Knerr, Stephanie Chen, Alexander Schmidt-Lebuhn

## What does this script do and what is it useful for?

This code was introduced in a submitted paper titled 'Phylogenomics-driven host test list selection for weed biological control'. It contains functions for calculating phylogenetic distance measures useful for creating host tests list in classical weed biological control.

The degree of relatedness between two taxa on a phylogeny is indicated by the number of nodes separating them. Here, we provide functions that calculate two distance measures, degree of separation i.e. node count and patristic distance, given an input phylogenetic tree.

## Descendant List Function

First define a function to recursively collect descendants of a node. This is used by the 'degreeofsep' function later on.

``` {.r .cell-code}
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

Function for calculating degrees of separation i.e. node count from a specified target weed.

``` {.r .cell-code}
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
# phylogenetic tree as newick file
treeFileName <- "astereae_concatenated.tre"
# specify outgroup(s)
taxonListForOutgroup <- c("Dimorphotheca_pluvialis", "Ewartia_nubigena", "Abrotanella_nivigena","Cotula_coronopifolia") 
# the target taxon to calculate distances from
myTargetTaxon <- "Erigeron_bonariensis"
# the output file name
outputFileName <- paste0("phylodists_", myTargetTaxon, ".tsv")
```

## Example Usage

``` {.r .cell-code}
# load libraries
# load libraries
library(ape)
library(adephylo)

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
Degsep <- myDegsep[which(row.names(myDegsep)==targetTaxon), ]
PatristicDist <- myPatristicMordered[which(row.names(myPatristicMordered) ==
                                            targetTaxon), ]

PhyloDists <- data.frame(Species = Terminal[order(PatristicDist)], 
                         DegSep = Degsep[order(PatristicDist)], 
                         PatristicDist = PatristicDist[order(PatristicDist)],
                         row.names = NULL)

write.table(PhyloDists, file = outputFileName, sep = "\t", row.names = FALSE)

knitr::kable(PhyloDists)
```

The output is a tab separated file with columns for the scientific name, degree of separation and patristic distance. See the file phylodists_Erigeron_bonariensis.tsv for an example that was used as a case study in the paper.

## Copyright and license information
Copyright (C) 2024  CSIRO

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Citation and contact information

Please cite the following paper if you use these scripts:

Stephanie H. Chen, Ben Gooden, Michelle A. Rafter, Gavin C. Hunter, Alicia Grealy, Nunzio Knerr, Alexander N. Schmidt-Lebuhn. *Phylogenomics-driven host test list selection for weed biological control*. Biological Control, Volume 193, 2024, 105529, <https://doi.org/10.1016/j.biocontrol.2024.105529>[\
\
](#0)Contact the corresponding author of the paper, [Alexander Schmidt-Lebuhn](mailto:alexander.s-l@csiro.au), if you have any questions.
