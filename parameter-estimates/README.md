The folder [parameter-estimates](parameter-estimates) contains files to fit the model
both as a pooled-data model
([parameter_estimates.r](parameter-estimates/parameter_estimates.r))
and one-at-a-time
([parameter_estimates_disjoint.r](parameter-estimates/parameter_estimates_disjoint.r)). 

For the paper, we took the parameter estimates based on the pooled-data model
with a uniform prior and calculated confidence intervals for the predictions
for all 138 analized trees. Furthermore, we compared the obtained annotations
with highest certainty (smallest CI closests to either 0 or 1) with the
available set of annotations on GO via a query with Quick-GO. Ultimately, this
concluded on generating a list of 295 proposed annotations (including absent/present).

The file used to generate this new list is in 
[proposed-annotations.r](proposed-annotations/proposed-annotations.r) and the
resulting set of annotations in [proposed-annotations.csv](proposed-annotations/proposed-annotations.csv)