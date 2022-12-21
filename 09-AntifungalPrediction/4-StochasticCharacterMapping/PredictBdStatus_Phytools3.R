#################################################
### Predict Bd Status of Salamander Sequences ###
#################################################

# Using phytools in parallel. Assume using 10 cores.

# For all taxa observed in non-control samples.

# Report start of script.
print(paste0("Starting script... (",Sys.time()," ",Sys.timezone(),")"))

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# Set working directory.
setwd("/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/AntifungalPrediction/Take2")

# Load microbe helper functions.
source("/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/HelperFunctions/MicrobeFunctions.R")

# Load phylogenetic tree.
library(phytools)
print(paste0("- Reading in tree: ",elapsed.time(ptm)))
tree<-read.tree("PhyloTree2")

# Set seed.
set.seed(1234)

# Root the tree.
print(paste0("- Rooting the tree: ",elapsed.time(ptm)))
tree2<-midpoint.root(tree)

# Resolve polytomies.
print(paste0("- Resolving polytomies: ",elapsed.time(ptm)))
tree3<-multi2di(tree2)

# Copy the tree and set edges of length zero to be 10^-6 times the total tree length.
print(paste0("- Resetting edges of zero length: ",elapsed.time(ptm)))
tree4<-tree3
tree4$edge.length[tree4$edge.length==0]<-max(nodeHeights(tree4))*1e-6

# Load the sequence metadaa.
print(paste0("- Reading in metadata: ",elapsed.time(ptm)))
meta<-read.csv(file="WoodhamsBdInhibitionStatus5.csv")

# Create Bd inhibition data.
print(paste0("- Creating Bd inhibition matrix: ",elapsed.time(ptm)))
BdInhib<-setNames(object=meta$Bd_Inhibition,nm=meta$TaxaID)
BdInhib<-to.matrix(x=BdInhib,seq=c("Inhibitory","Non-Inhibitory"))
BdInhib[meta$Dataset=="Salamanders",]<-rep(0.5,2)
colnames(BdInhib)[2]<-"NonInhibitory"

# From: http://blog.phytools.org/2017/11/running-makesimmap-in-parallel.html

# First fit our model.
print(paste0("- Fitting first model: ",elapsed.time(ptm)))
fit<-fitMk(tree=tree4,x=BdInhib,model="ER")

# Extract the fitted transition matrix.
print(paste0("- Extracting the transition matrix: ",elapsed.time(ptm)))
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states

# Run the analysis.
print(paste0("- Running the analysis: ",elapsed.time(ptm)))
library(parallel)
trees3<-mclapply(1:10,function(n,tree,x,fixedQ) make.simmap(tree,x,Q=fixedQ,nsim=500),
                tree=tree4,x=BdInhib,fixedQ=fittedQ,
                mc.cores=if(.Platform$OS.type=="windows") 1L else 10L)

# Combine trees.
print(paste0("- Combining trees: ",elapsed.time(ptm)))
trees3<-do.call(c,trees3)
if(!("multiSimmap"%in%class(trees3))) class(trees3)<-c("multiSimmap",class(trees3))

# Saving output.
print(paste0("- Saving output: ",elapsed.time(ptm)))
save(trees3,file="trees3.RData")

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
