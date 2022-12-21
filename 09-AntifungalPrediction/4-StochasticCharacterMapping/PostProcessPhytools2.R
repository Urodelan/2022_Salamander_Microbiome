# Get tip states from trees processed with phytools.

# Load trees.
load("trees3.RData")

# Describe trees.
library(phytools)
description<-describe.simmap(trees3)

# Get tip states.
BdTaxaStates<-description$tips

# Save tip states as a csv file.
write.csv(BdTaxaStates,file="BdTaxaStates2.csv",row.names=T)

# Done!
