# Trees

This data excludes the trees. Phylogenetic trees can be downloaded from PantherDB
website. We used the 15.0 version [here](ftp://ftp.pantherdb.org/panther_library/15.0/)

# Annotations

`true_annotations.gz` contains all the functional annotations through manual curation. It is the one you should use in your model.
`inferred_annotations.gz` contains all the inferred annotations from the phylogenetic curation project. It is something you can use to compare with your inferences.

Both files have 4 columns:
Col 1  Node id, for example “PTHR10005:AN136”. 
Col 2  Leaf sequence ID, for example "HUMAN|HGNC=10896|UniProtKB=P12755”.
Col 3  GO id, for example "GO:0000122”. 
Col 4  Annotation qualifier (NOT). 

