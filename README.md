# On the automatic annotation of gene functions usingobservational data and phylogenetic trees

This repository contains all the materials that are necesary to reproduce the
figures and tables of the paper, including a set of novel predictions.

There are two main sections: the simulations folder and the parameter-estimates
folder. The first deals with the simulation study generating random annotated
phylogenetic trees using the PANTHER database, the latter fits the pooled-data
model using 141 different sets of annotations and makes the predictions, all
using a combination between the GO dataset and the PANTHER phylogenetic trees.

## Data

For this paper we used the data available in PANTHER version 15.0.


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
[proposed-annotations.r](parameter-estimates/proposed-annotations.r) and the
resulting set of annotations in [proposed-annotations.csv](parameter-estimates/proposed-annotations.csv)


