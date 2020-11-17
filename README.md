# On the automatic annotation of gene functions usingobservational data and phylogenetic trees

This repository contains all the materials that are necesary to reproduce the
figures and tables of the paper, including a set of novel predictions.

There are two main sections: the simulations folder and the parameter-estimates
folder. The first deals with the simulation study generating random annotated
phylogenetic trees using the PANTHER database, the latter fits the pooled-data
model using 138 different sets of annotations and makes the predictions, all
using a combination between the GO dataset and the PANTHER phylogenetic trees.

## The aphylo package

All of the methods presented in this paper are available in the R package
[`aphylo`](https://github.com/USCbiostats/aphylo). To install the `aphylo`
package, you should use the following command:

```r
devtools::install_github("USCbiostats/aphylo")
```

## Overview of the repository

In general, each dataset, figure, or table has its own R script used to be generated.
Furthermore, resulting files have the same name of the R script that created it, only
changing in the extension, for example, the R script `candidate_trees.r` generates the file
`candidate_trees.rds`, which has the 138 phylogenetic trees (including annotations) used
throughtout the paper.

- [**data-raw**](data-raw) Contains raw data used in the paper. The main dataset here is the
set of experimental annotations from GOA.

- [**data**](data) Contains the scripts used to process the data, including, reading
the panther trees, GOA annotations, and combining them into `aphylo_tree`
objects. This also contains the resulting data.

- [**fig**](fig) Most figures of the paper, including the code used to generate the two trees featured
in the paper (the low and high MAE).

- [**parameter-estimates**](parameter-estimates) All the code used to fit the pooled models, including the
obtained paramter estimates.

- [**proposed-annotations**](proposed-annotations) Code used to generate annotation proposals on genes
with no experimental annotations. The folder includes the actual table
with proposed annotations.

- [**sifter**](sifter) Code used to analyze SIFTER data.

- [**simulations**](simulations) Code to generate the large simulation study in which we analyze
the properties of MCMC and MLE estimates.

