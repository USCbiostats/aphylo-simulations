# On the automatic annotation of gene functions usingobservational data and phylogenetic trees

This repository contains all the materials that are necesary to reproduce the
figures and tables of the paper, including a set of novel predictions.

There are two main sections: the simulations folder and the parameter-estimates
folder. The first deals with the simulation study generating random annotated
phylogenetic trees using the PANTHER database, the latter fits the pooled-data
model using 138 different sets of annotations and makes the predictions, all
using a combination between the GO dataset and the PANTHER phylogenetic trees.

## Overview of the folders:

**data-raw** Contains raw data used in the paper. The main dataset here is the
set of experimental annotations from GOA.

**data** Contains the scripts used to process the data, including, reading
the panther trees, GOA annotations, and combining them into `aphylo_tree`
objects. This also contains the resulting data.

**fig** Most figures of the paper, including the code used to generate the two trees featured
in the paper (the low and high MAE).

**parameter-estimates** All the code used to fit the pooled models, including the
obtained paramter estimates.

**proposed-annotations** Code used to generate annotation proposals on genes
with no experimental annotations. The folder includes the actual table
with proposed annotations.

**sifter** Code used to analyze SIFTER data.

**simulations** Code to generate the large simulation study in which we analyze
the properties of MCMC and MLE estimates.

## Data

For this paper we used the data available in PANTHER version 14.1, which can
be downloaded [here](ftp://ftp.pantherdb.org/panther_library/14.1/).

The experimental annotations were retrieved from the Gene Ontology release
2019-02-01, which can be downloaded [here](https://zenodo.org/record/2555603).

To ease reproducibility, we have included all experimental annotations used
for this paper, including the panther trees, in the [data](data) folder.

The [data-raw](data-raw) folder includes information 


## Simulation study

The simulation study consisted on generating annotations on the full set of the
trees available in the PANTHER project. Files can be found [here](simulations).

## Figures

Most of the figures of the paper can be found [here](fig)

## Application on GO annotations w/ PANTHER

The folder [parameter-estimates](parameter-estimates) contains files to fit the model
both as a pooled-data model
([parameter_estimates.r](parameter-estimates/parameter_estimates.r))
and one-at-a-time
([parameter_estimates_disjoint.r](parameter-estimates/parameter_estimates_disjoint.r)). 

For the paper, we took the parameter estimates based on the pooled-data model
with a uniform prior and calculated confidence intervals for the predictions
for all 141 analized trees. Furthermore, we compared the obtained annotations
with highest certainty (smallest CI closests to either 0 or 1) with the
available set of annotations on GO via a query with Quick-GO. Ultimately, this
concluded on generating a list of 295 proposed annotations (including absent/present).

The file used to generate this new list is in 
[proposed-annotations.r](proposed-annotations/proposed-annotations.r) and the
resulting set of annotations in [proposed-annotations.csv](proposed-annotations/proposed-annotations.csv)


