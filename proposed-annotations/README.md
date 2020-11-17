# Proposed annotations

For the paper, we took the parameter estimates based on the pooled-data model
with a uniform prior and calculated confidence intervals for the predictions
for all 138 analized trees. Furthermore, we compared the obtained annotations
with highest certainty (smallest CI closests to either 0 or 1) with the
available set of annotations on GO via a query with Quick-GO. Ultimately, this
concluded on generating a list of 225 proposed annotations (including absent/present).

The main file is [proposed-annotations.csv](proposed-annotations.csv), which has the following columns:

- UniProtKB

- go_id The id of the GO term.

- p0_025, p0_500, p0_975 the 2.5%, 50%, and 97.5% centiles of the posterior probabilities.

- qualifier_binary The observed Qualifier in the GO database

- go_evidence_code The available evidence codes for this annotation (if any)

- qualifier_proposed The proposed qualifier

- opposing when `TRUE` there are opposing annotations in GOA

