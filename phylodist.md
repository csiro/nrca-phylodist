---
author:
- Nunzio Knerr
- Alexander Schmidt-Lebuhn
- Stephanie Chen
authors:
- name: Nunzio Knerr
- name: Alexander Schmidt-Lebuhn
- name: Stephanie Chen
title: Phylogenetic Distance
toc-title: Table of contents
---

## Phylogenetic Distance/Separation

This code contains the functions for calculating the phylogenetic
separation as counts of nodes from any given taxon to a target taxon in
a phylogenetic tree. The resulting output is a matrix that contains a
column with the node count from the target taxon. This is very useful
when comparing how closely related any given taxon is to the target
taxon in the tree.

## Descendant List Function

First define a function to recursively collect descendants of a node.
This is used buy the 'degreeofsep' function later on.

::: cell
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
:::

## Degree Of Separation Function

::: cell
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
:::

## User Input Variables

Specify the inputs and outputs for use in the script

::: cell
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
:::

## Example Usage

::: cell
``` {.r .cell-code}
# load libraries
library(ape)
library(adephylo)
```

::: {.cell-output .cell-output-stderr}
    Loading required package: ade4
:::

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

::: cell-output-display
  -------------------------------------------------------------------------------------------------------------------------------------
                                  Terminal.order.PatristicDist..     Degsep.order.PatristicDist..   PatristicDist.order.PatristicDist..
  ------------------------------- -------------------------------- ------------------------------ -------------------------------------
  Erigeron_bonariensis            Erigeron_bonariensis                                          0                             0.0000000

  Erigeron_sumatrensis            Erigeron_sumatrensis                                          0                             0.0261555

  Erigeron_bilbaoanus             Erigeron_bilbaoanus                                           1                             0.0291839

  Erigeron_canadensis             Erigeron_canadensis                                           1                             0.0330095

  Erigeron_primulifolius          Erigeron_primulifolius                                        3                             0.0375416

  Erigeron_conyzoides             Erigeron_conyzoides                                           2                             0.0404355

  Erigeron_karvinskianus          Erigeron_karvinskianus                                        2                             0.0408097

  Olearia_imbricata               Olearia_imbricata                                             7                             0.0543982

  Calotis_squamigera              Calotis_squamigera                                            7                             0.0550798

  Pappochroma_setosus             Pappochroma_setosus                                           7                             0.0555396

  Olearia_laciniifolia            Olearia_laciniifolia                                          7                             0.0559845

  Baccharis_pingraea              Baccharis_pingraea                                            5                             0.0562133

  Iotasperma_sessilifolium        Iotasperma_sessilifolium                                      7                             0.0564643

  Baccharis_halimifolia           Baccharis_halimifolia                                         5                             0.0567125

  Tetramolopium_sp_Mt_Bowen       Tetramolopium_sp_Mt_Bowen                                     7                             0.0576629

  Dichromochlamys_dentatifolia    Dichromochlamys_dentatifolia                                  7                             0.0583778

  Pappochroma_paludicola          Pappochroma_paludicola                                        7                             0.0588579

  Lagenophora_stipitata           Lagenophora_stipitata                                         7                             0.0590056

  Olearia_elaeophila              Olearia_elaeophila                                            7                             0.0590869

  Solidago_canadensis             Solidago_canadensis                                           4                             0.0592591

  Olearia_ramulosa                Olearia_ramulosa                                              7                             0.0592983

  Olearia_muricata                Olearia_muricata                                              7                             0.0593151

  Pappochroma_nitidus             Pappochroma_nitidus                                           7                             0.0593572

  Olearia_phlogopappa             Olearia_phlogopappa                                           7                             0.0595578

  Olearia_lirata                  Olearia_lirata                                                7                             0.0595767

  Achnophora_tatei                Achnophora_tatei                                              7                             0.0596219

  Pappochroma_tasmanicum          Pappochroma_tasmanicum                                        7                             0.0596471

  Erodiophyllum_acanthocephalum   Erodiophyllum_acanthocephalum                                 7                             0.0596626

  Olearia_iodochroa               Olearia_iodochroa                                             7                             0.0597435

  Olearia_canescens               Olearia_canescens                                             7                             0.0599336

  Solidago_chilensis              Solidago_chilensis                                            4                             0.0600133

  Lagenophora_huegelii            Lagenophora_huegelii                                          7                             0.0600429

  Camptacra_barbata               Camptacra_barbata                                             7                             0.0603631

  Pilbara_trudgenii               Pilbara_trudgenii                                             7                             0.0603928

  Olearia_astroloba               Olearia_astroloba                                             7                             0.0604746

  Vittadinia_muelleri             Vittadinia_muelleri                                           7                             0.0606977

  Olearia_rudis                   Olearia_rudis                                                 7                             0.0608804

  Olearia_ericoides               Olearia_ericoides                                             7                             0.0608955

  Solidago_sempervirens           Solidago_sempervirens                                         4                             0.0609066

  Solidago_altissima              Solidago_altissima                                            4                             0.0609625

  Pappochroma_bellidioides        Pappochroma_bellidioides                                      7                             0.0609645

  Olearia_algida                  Olearia_algida                                                7                             0.0609823

  Calotis_pubescens               Calotis_pubescens                                             7                             0.0610436

  Calotis_cuneifolia              Calotis_cuneifolia                                            7                             0.0612126

  Calotis_plumulifera             Calotis_plumulifera                                           7                             0.0612460

  Lagenophora_gracilis            Lagenophora_gracilis                                          7                             0.0612657

  Ixiochlamys_nana                Ixiochlamys_nana                                              7                             0.0612743

  Olearia_picridifolia            Olearia_picridifolia                                          7                             0.0613071

  Olearia_muelleri                Olearia_muelleri                                              7                             0.0613912

  Calotis_lappulacea              Calotis_lappulacea                                            7                             0.0614783

  Pappochroma_gunnii              Pappochroma_gunnii                                            7                             0.0614913

  Vittadinia_tenuissima           Vittadinia_tenuissima                                         7                             0.0614996

  Olearia_gordonii                Olearia_gordonii                                              7                             0.0615324

  Olearia_eremaea                 Olearia_eremaea                                               7                             0.0615434

  Callistephus_chinensis          Callistephus_chinensis                                        7                             0.0615576

  Lagenophora_montana             Lagenophora_montana                                           7                             0.0615746

  Olearia_dampieri                Olearia_dampieri                                              7                             0.0616293

  Olearia_arida                   Olearia_arida                                                 7                             0.0617296

  Olearia_frostii                 Olearia_frostii                                               7                             0.0617338

  Olearia_stellulata              Olearia_stellulata                                            7                             0.0617574

  Olearia_lanuginosa              Olearia_lanuginosa                                            7                             0.0618187

  Olearia_brachyphylla            Olearia_brachyphylla                                          7                             0.0618603

  Ceratogyne_obionoides           Ceratogyne_obionoides                                         7                             0.0619580

  Solenogyne_gunnii               Solenogyne_gunnii                                             7                             0.0620782

  Vittadinia_eremaea              Vittadinia_eremaea                                            7                             0.0621175

  Olearia_heterocarpa             Olearia_heterocarpa                                           7                             0.0622672

  Olearia_cordata                 Olearia_cordata                                               7                             0.0624459

  Olearia_burgessii               Olearia_burgessii                                             7                             0.0625578

  Vittadinia_humerata             Vittadinia_humerata                                           7                             0.0625778

  Calotis_erinacea                Calotis_erinacea                                              7                             0.0626284

  Olearia_humilis                 Olearia_humilis                                               7                             0.0626831

  Calotis_inermis                 Calotis_inermis                                               7                             0.0627982

  Olearia_lepidophylla            Olearia_lepidophylla                                          7                             0.0628116

  Ixiochlamys_cuneifolia          Ixiochlamys_cuneifolia                                        7                             0.0629836

  Olearia_heloderma               Olearia_heloderma                                             7                             0.0629884

  Olearia_decurrens               Olearia_decurrens                                             7                             0.0630876

  Calotis_porphyroglossa          Calotis_porphyroglossa                                        7                             0.0631434

  Calotis_hispidula               Calotis_hispidula                                             7                             0.0631677

  Olearia_exiguifolia             Olearia_exiguifolia                                           7                             0.0632015

  Olearia_microdisca              Olearia_microdisca                                            7                             0.0632349

  Olearia_adenolasia              Olearia_adenolasia                                            7                             0.0632539

  Olearia_revoluta                Olearia_revoluta                                              7                             0.0632651

  Olearia_lasiophylla             Olearia_lasiophylla                                           7                             0.0632988

  Olearia_ramosissima             Olearia_ramosissima                                           7                             0.0633501

  Olearia_stuartii                Olearia_stuartii                                              7                             0.0634107

  Dimorphocoma_minutula           Dimorphocoma_minutula                                         7                             0.0635802

  Olearia_suffruticosa            Olearia_suffruticosa                                          7                             0.0635967

  Olearia_adenophora              Olearia_adenophora                                            7                             0.0636735

  Vittadinia_gracilis             Vittadinia_gracilis                                           7                             0.0637866

  Olearia_asterotricha            Olearia_asterotricha                                          7                             0.0638010

  Olearia_calcarea                Olearia_calcarea                                              7                             0.0638715

  Vittadinia_sulcata              Vittadinia_sulcata                                            7                             0.0639150

  Symphyotrichum_novi-belgii      Symphyotrichum_novi-belgii                                    4                             0.0641152

  Olearia_glandulosa              Olearia_glandulosa                                            7                             0.0641304

  Vittadinia_sericea              Vittadinia_sericea                                            7                             0.0641346

  Vittadinia_obovata              Vittadinia_obovata                                            7                             0.0641518

  Calotis_glandulosa              Calotis_glandulosa                                            7                             0.0641648

  Heterotheca_grandiflora         Heterotheca_grandiflora                                       4                             0.0642155

  Vittadinia_scabra               Vittadinia_scabra                                             7                             0.0642202

  Vittadinia_nullarborensis       Vittadinia_nullarborensis                                     7                             0.0642270

  Olearia_tenuifolia              Olearia_tenuifolia                                            7                             0.0642515

  Calotis_spec                    Calotis_spec                                                  7                             0.0643719

  Olearia_ballii                  Olearia_ballii                                                7                             0.0643862

  Olearia_passerinoides           Olearia_passerinoides                                         7                             0.0643986

  Vittadinia_bicolor              Vittadinia_bicolor                                            7                             0.0644209

  Vittadinia_cuneata              Vittadinia_cuneata                                            7                             0.0645438

  Olearia_rotundifolia            Olearia_rotundifolia                                          7                             0.0645679

  Olearia_sp_Mt_Edward            Olearia_sp_Mt_Edward                                          7                             0.0646577

  Olearia_sp_Waterhouse_Range     Olearia_sp_Waterhouse_Range                                   7                             0.0647168

  Olearia_incana                  Olearia_incana                                                7                             0.0647630

  Olearia_subspicata              Olearia_subspicata                                            7                             0.0647790

  Vittadinia_blackii              Vittadinia_blackii                                            7                             0.0648817

  Olearia_magniflora              Olearia_magniflora                                            7                             0.0649962

  Calotis_cuneata                 Calotis_cuneata                                               7                             0.0650201

  Calotis_breviradiata            Calotis_breviradiata                                          7                             0.0650778

  Olearia_montana                 Olearia_montana                                               7                             0.0651623

  Calotis_cymbacantha             Calotis_cymbacantha                                           7                             0.0651632

  Olearia_stenophylla             Olearia_stenophylla                                           7                             0.0651938

  Calotis_moorei                  Calotis_moorei                                                7                             0.0652256

  Calotis_breviseta               Calotis_breviseta                                             7                             0.0652447

  Olearia_hygrophila              Olearia_hygrophila                                            7                             0.0652709

  Vittadinia_pterochaeta          Vittadinia_pterochaeta                                        7                             0.0653444

  Vittadinia_dissecta             Vittadinia_dissecta                                           7                             0.0653483

  Olearia_brevipedunculata        Olearia_brevipedunculata                                      7                             0.0653766

  Solenogyne_dominii              Solenogyne_dominii                                            7                             0.0655776

  Olearia_hookeri                 Olearia_hookeri                                               7                             0.0655927

  Olearia_homolepis               Olearia_homolepis                                             7                             0.0656471

  Calotis_anthemoides             Calotis_anthemoides                                           7                             0.0656992

  Calotis_latiuscula              Calotis_latiuscula                                            7                             0.0657039

  Olearia_paucidentata            Olearia_paucidentata                                          7                             0.0659540

  Ixiochlamys_filicifolia         Ixiochlamys_filicifolia                                       7                             0.0659836

  Vittadinia_australasica         Vittadinia_australasica                                       7                             0.0660363

  Vittadinia_constricta           Vittadinia_constricta                                         7                             0.0660453

  Olearia_tubuliflora             Olearia_tubuliflora                                           7                             0.0660838

  Calotis_scabiosifolia           Calotis_scabiosifolia                                         7                             0.0661289

  Calotis_multicaulis             Calotis_multicaulis                                           7                             0.0661511

  Olearia_axillaris               Olearia_axillaris                                             7                             0.0661861

  Olearia_quercifolia             Olearia_quercifolia                                           7                             0.0662368

  Calotis_dentex                  Calotis_dentex                                                7                             0.0662432

  Olearia_minor                   Olearia_minor                                                 7                             0.0663110

  Olearia_nernstii                Olearia_nernstii                                              7                             0.0664234

  Calotis_scapigera               Calotis_scapigera                                             7                             0.0664831

  Vittadinia_burbidgeae           Vittadinia_burbidgeae                                         7                             0.0664915

  Vittadinia_cervicularis         Vittadinia_cervicularis                                       7                             0.0665686

  Vittadinia_hispidula            Vittadinia_hispidula                                          7                             0.0667640

  Olearia_plucheacea              Olearia_plucheacea                                            7                             0.0668257

  Olearia_cassiniae               Olearia_cassiniae                                             7                             0.0671698

  Olearia_ciliata                 Olearia_ciliata                                               7                             0.0672107

  Olearia_tomentosa               Olearia_tomentosa                                             7                             0.0672827

  Vittadinia_pustulata            Vittadinia_pustulata                                          7                             0.0673213

  Olearia_propinqua               Olearia_propinqua                                             7                             0.0674967

  Isoetopsis_graminifolia         Isoetopsis_graminifolia                                       7                             0.0676389

  Olearia_gravis                  Olearia_gravis                                                7                             0.0677839

  Solenogyne_bellioides           Solenogyne_bellioides                                         7                             0.0678598

  Minuria_multiseta               Minuria_multiseta                                             7                             0.0679362

  Olearia_pimeleoides             Olearia_pimeleoides                                           7                             0.0679892

  Kippistia_suaedifolia           Kippistia_suaedifolia                                         7                             0.0680246

  Calotis_xanthosioidea           Calotis_xanthosioidea                                         7                             0.0680945

  Olearia_microphylla             Olearia_microphylla                                           7                             0.0684939

  Olearia_rugosa                  Olearia_rugosa                                                7                             0.0685557

  Olearia_floribunda              Olearia_floribunda                                            7                             0.0688581

  Calotis_ancyrocarpa             Calotis_ancyrocarpa                                           7                             0.0689911

  Vittadinia_condyloides          Vittadinia_condyloides                                        7                             0.0690272

  Olearia_ferresii                Olearia_ferresii                                              7                             0.0704175

  Calotis_suffruticosa            Calotis_suffruticosa                                          7                             0.0704827

  Minuria_cunninghamii            Minuria_cunninghamii                                          7                             0.0705773

  Elachanthus_pusillus            Elachanthus_pusillus                                          7                             0.0706080

  Minuria_tridens                 Minuria_tridens                                               7                             0.0710899

  Olearia_glutinosa               Olearia_glutinosa                                             7                             0.0715413

  Minuria_rigida                  Minuria_rigida                                                7                             0.0716783

  Vittadinia_megacephala          Vittadinia_megacephala                                        7                             0.0717702

  Minuria_gardneri                Minuria_gardneri                                              7                             0.0719858

  Minuria_leptophylla             Minuria_leptophylla                                           7                             0.0725579

  Pembertonia_latisquamea         Pembertonia_latisquamea                                       7                             0.0744824

  Allittia_cardiocarpa            Allittia_cardiocarpa                                          7                             0.0749899

  Brachyscome_tesquorum           Brachyscome_tesquorum                                         7                             0.0783672

  Brachyscome_procumbens          Brachyscome_procumbens                                        7                             0.0791969

  Brachyscome_obovata             Brachyscome_obovata                                           7                             0.0802079

  Brachyscome_petrophila          Brachyscome_petrophila                                        7                             0.0811156

  Brachyscome_barkerae            Brachyscome_barkerae                                          7                             0.0812040

  Brachyscome_graminea            Brachyscome_graminea                                          7                             0.0813418

  Brachyscome_dichromosomatica    Brachyscome_dichromosomatica                                  7                             0.0814912

  Brachyscome_basaltica           Brachyscome_basaltica                                         7                             0.0814979

  Brachyscome_linearifolia        Brachyscome_linearifolia                                      7                             0.0817175

  Brachyscome_tamworthensis       Brachyscome_tamworthensis                                     7                             0.0820061

  Brachyscome_wingello            Brachyscome_wingello                                          7                             0.0821678

  Brachyscome_papillosa           Brachyscome_papillosa                                         7                             0.0824993

  Brachyscome_tetrapterocarpa     Brachyscome_tetrapterocarpa                                   7                             0.0831951

  Brachyscome_melanocarpa         Brachyscome_melanocarpa                                       7                             0.0835014

  Brachyscome_scapigera           Brachyscome_scapigera                                         7                             0.0837817

  Brachyscome_spathulata          Brachyscome_spathulata                                        7                             0.0838394

  Brachyscome_stolonifera         Brachyscome_stolonifera                                       7                             0.0839217

  Brachyscome_nova-anglica        Brachyscome_nova-anglica                                      7                             0.0843874

  Brachyscome_willisii            Brachyscome_willisii                                          7                             0.0845908

  Brachyscome_paludicola          Brachyscome_paludicola                                        7                             0.0846877

  Brachyscome_sinclairii          Brachyscome_sinclairii                                        7                             0.0848743

  Brachyscome_ptychocarpa         Brachyscome_ptychocarpa                                       7                             0.0849860

  Brachyscome_ascendens           Brachyscome_ascendens                                         7                             0.0852916

  Brachyscome_chrysoglossa        Brachyscome_chrysoglossa                                      7                             0.0855822

  Brachyscome_nodosa              Brachyscome_nodosa                                            7                             0.0856969

  Brachyscome_diversifolia        Brachyscome_diversifolia                                      7                             0.0857456

  Brachyscome_salkiniae           Brachyscome_salkiniae                                         7                             0.0858459

  Brachyscome_rara                Brachyscome_rara                                              7                             0.0858727

  Brachyscome_dissectifolia       Brachyscome_dissectifolia                                     7                             0.0859569

  Brachyscome_segmentosa          Brachyscome_segmentosa                                        7                             0.0859746

  Brachyscome_breviscapis         Brachyscome_breviscapis                                       7                             0.0859940

  Brachyscome_aculeata            Brachyscome_aculeata                                          7                             0.0860449

  Brachyscome_microcarpa          Brachyscome_microcarpa                                        7                             0.0862575

  Brachyscome_smithwhitei         Brachyscome_smithwhitei                                       7                             0.0862758

  Brachyscome_multifida           Brachyscome_multifida                                         7                             0.0863625

  Brachyscome_decipiens           Brachyscome_decipiens                                         7                             0.0864810

  Brachyscome_curvicarpa          Brachyscome_curvicarpa                                        7                             0.0865544

  Brachyscome_dentata             Brachyscome_dentata                                           7                             0.0865845

  Chondropyxis_halophila          Chondropyxis_halophila                                        8                             0.0871172

  Brachyscome_trachycarpa         Brachyscome_trachycarpa                                       7                             0.0873506

  Brachyscome_radicans            Brachyscome_radicans                                          7                             0.0875452

  Brachyscome_sieberi             Brachyscome_sieberi                                           7                             0.0880030

  Brachyscome_stuartii            Brachyscome_stuartii                                          7                             0.0880690

  Brachyscome_kaputarensis        Brachyscome_kaputarensis                                      7                             0.0883164

  Brachyscome_perpusilla          Brachyscome_perpusilla                                        7                             0.0888728

  Brachyscome_readeri             Brachyscome_readeri                                           7                             0.0888807

  Brachyscome_tadgellii           Brachyscome_tadgellii                                         7                             0.0897006

  Bellis_perennis                 Bellis_perennis                                               6                             0.0898949

  Brachyscome_cuneifolia          Brachyscome_cuneifolia                                        7                             0.0899630

  Brachyscome_rigidula            Brachyscome_rigidula                                          7                             0.0906384

  Brachyscome_dalbyensis          Brachyscome_dalbyensis                                        7                             0.0913452

  Brachyscome_casstiana           Brachyscome_casstiana                                         7                             0.0916356

  Brachyscome_eriogona            Brachyscome_eriogona                                          7                             0.0925592

  Brachyscome_ciliaris            Brachyscome_ciliaris                                          7                             0.0926969

  Brachyscome_debilis             Brachyscome_debilis                                           7                             0.0928508

  Brachyscome_tatei               Brachyscome_tatei                                             7                             0.0931240

  Brachyscome_iberidifolia        Brachyscome_iberidifolia                                      7                             0.0942351

  Celmisia_asteliifolia           Celmisia_asteliifolia                                         9                             0.0944826

  Brachyscome_exilis              Brachyscome_exilis                                            7                             0.0945725

  Celmisia_lyallii                Celmisia_lyallii                                              9                             0.0948329

  Brachyscome_parvula             Brachyscome_parvula                                           7                             0.0951245

  Brachyscome_ciliocarpa          Brachyscome_ciliocarpa                                        7                             0.0954911

  Brachyscome_cheilocarpa         Brachyscome_cheilocarpa                                       7                             0.0960027

  Olearia_speciosa                Olearia_speciosa                                              9                             0.0970135

  Brachyscome_billabongensis      Brachyscome_billabongensis                                    7                             0.0974063

  Olearia_avicenniaefolia         Olearia_avicenniaefolia                                       9                             0.0984200

  Olearia_lanceolata              Olearia_lanceolata                                            9                             0.0987004

  Allittia_uliginosa              Allittia_uliginosa                                            7                             0.0987396

  Brachyscome_gilesii             Brachyscome_gilesii                                           7                             0.0989206

  Celmisia_latifolia              Celmisia_latifolia                                            9                             0.1004190

  Olearia_tasmanica               Olearia_tasmanica                                             9                             0.1007450

  Celmisia_saxifraga              Celmisia_saxifraga                                            9                             0.1008971

  Olearia_viscosa                 Olearia_viscosa                                               9                             0.1012949

  Olearia_arborescens             Olearia_arborescens                                           9                             0.1013846

  Celmisia_discolor               Celmisia_discolor                                             9                             0.1015063

  Pleurophyllum_hookeri           Pleurophyllum_hookeri                                         9                             0.1021389

  Celmisia_bellidioides           Celmisia_bellidioides                                         9                             0.1025304

  Olearia_archeri                 Olearia_archeri                                               9                             0.1027637

  Celmisia_costiniana             Celmisia_costiniana                                           9                             0.1031362

  Celmisia_pugioniformis          Celmisia_pugioniformis                                        9                             0.1033010

  Celmisia_spectabilis            Celmisia_spectabilis                                          9                             0.1033617

  Olearia_nummularifolia          Olearia_nummularifolia                                        9                             0.1033921

  Olearia_rosmarinifolia          Olearia_rosmarinifolia                                        9                             0.1035887

  Celmisia_gracilenta             Celmisia_gracilenta                                           9                             0.1038710

  Olearia_argophylla              Olearia_argophylla                                            9                             0.1039109

  Celmisia_longifolia             Celmisia_longifolia                                           9                             0.1046570

  Olearia_cymbifolia              Olearia_cymbifolia                                            9                             0.1049913

  Olearia_pinifolia               Olearia_pinifolia                                             9                             0.1051978

  Olearia_oppositifolia           Olearia_oppositifolia                                         9                             0.1053397

  Olearia_obcordata               Olearia_obcordata                                             9                             0.1053607

  Olearia_ledifolia               Olearia_ledifolia                                             9                             0.1055517

  Olearia_cydoniifolia            Olearia_cydoniifolia                                          9                             0.1056628

  Olearia_grandiflora             Olearia_grandiflora                                           9                             0.1065992

  Olearia_alpicola                Olearia_alpicola                                              9                             0.1066653

  Olearia_pannosa                 Olearia_pannosa                                               9                             0.1069035

  Olearia_solandri                Olearia_solandri                                              9                             0.1069250

  Olearia_covenyi                 Olearia_covenyi                                               9                             0.1069403

  Celmisia_tomentella             Celmisia_tomentella                                           9                             0.1069634

  Olearia_stilwelliae             Olearia_stilwelliae                                           9                             0.1084522

  Olearia_erubescens              Olearia_erubescens                                            9                             0.1089169

  Olearia_myrsinoides             Olearia_myrsinoides                                           9                             0.1095844

  Olearia_megalophylla            Olearia_megalophylla                                          9                             0.1104366

  Olearia_fimbriata               Olearia_fimbriata                                             9                             0.1108652

  Brachyscome_bellidioides        Brachyscome_bellidioides                                      7                             0.1133181

  Ewartia_nubigena                Ewartia_nubigena                                             10                             0.1741899

  Dimorphotheca_pluvialis         Dimorphotheca_pluvialis                                      10                             0.1999102

  Cotula_coronopifolia            Cotula_coronopifolia                                         10                             0.2314641

  Abrotanella_nivigena            Abrotanella_nivigena                                         10                             0.2347103
  -------------------------------------------------------------------------------------------------------------------------------------
:::
:::
