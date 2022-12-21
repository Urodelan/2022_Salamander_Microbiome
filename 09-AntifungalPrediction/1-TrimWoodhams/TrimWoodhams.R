####################################################
### Trim Woodhams Sequences to Our Region of 16S ###
####################################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in helper functions.
source("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/HelperFunctions/MicrobeFunctions.R")

# Read in seqence data.
seq<-read.fasta("../AllNonCntrlTaxa/WoodhamsAndSalamandersUniqueSequences3.fasta")

# Read in sequence metadata.
meta<-read.csv("../AllNonCntrlTaxa/WoodhamsBdInhibitionStatus2.csv",stringsAsFactors=F)

# Append field for dataset to sequences.
seq$Dataset<-meta$Dataset[match(seq$Name,meta$TaxaID)]

# Get just Woodhams sequences.
w<-subset(seq,Dataset=="Woodhams")

# Check histogram of Woodhams sequence lengths.
# tiff("WoodhamsLengthHist.tiff",width=7,height=5,units="in",res=300)
hist(unname(sapply(X=w$Sequence,FUN=nchar)),xlab="Number of Basepairs",main="Histogram of Woodhams Sequence Lengths")
# dev.off()

# Define 16S amplification primers.
Forward_amp_primer16S<-"GTGYCAGCMGCCGCGGTAA"
library(Biostrings)
## Getting reverse compliment of reverse amplification primer.
Reverse_amp_primer16S<-as.character(reverseComplement(DNAString("GGACTACHVGGGTWTCTAAT")))

# Define IUPAC wildcard characters.
IUPAC_wildcards<-list(R=c("A","G"),Y=c("C","T"),S=c("G","C"),W=c("A","T"),K=c("G","T"),M=c("A","C"),B=c("C","G","T"),D=c("A","G","T"),H=c("A","C","T"),V=c("A","C","G"),N=c("A","T","G","C"))

# Expand amplification primers.
## Expand forward primer.
## Split up characters in string.
splitCharacters<-strsplit(Forward_amp_primer16S,split="")[[1]]
## Get just wildcard characters from string.
wildCharacters<-splitCharacters[splitCharacters %in% names(IUPAC_wildcards)]
## Create empty list to store wildcard character options.
options<-vector(mode="list",length=length(wildCharacters))
names(options)<-wildCharacters
## Fill in list with wildcard character options.
for(i in 1:length(options)){
  options[[i]]<-IUPAC_wildcards[[names(options)[i]]]
}
## Get all combinations of wildcard character options.
Combn<-expand.grid(options)
## Create a vector of possible forward primer strings.
primerOptionsForward<-c()
for(i in 1:nrow(Combn)){
  splitCharacters_copy<-as.character(splitCharacters)
  splitCharacters_copy[splitCharacters_copy %in% names(options)]<-as.character(unname(unlist(Combn[i,])))
  primerOptionsForward<-c(primerOptionsForward,paste0(splitCharacters_copy,collapse=""))
}
## Expand reverse primer.
## Split up characters in string.
splitCharacters<-strsplit(Reverse_amp_primer16S,split="")[[1]]
## Get just wildcard characters from string.
wildCharacters<-splitCharacters[splitCharacters %in% names(IUPAC_wildcards)]
## Create empty list to store wildcard character options.
options<-vector(mode="list",length=length(wildCharacters))
names(options)<-wildCharacters
## Fill in list with wildcard character options.
for(i in 1:length(options)){
  options[[i]]<-IUPAC_wildcards[[names(options)[i]]]
}
## Get all combinations of wildcard character options.
Combn<-expand.grid(options)
## Create a vector of possible reverse primer strings.
primerOptionsReverse<-c()
for(i in 1:nrow(Combn)){
  splitCharacters_copy<-as.character(splitCharacters)
  splitCharacters_copy[splitCharacters_copy %in% names(options)]<-as.character(unname(unlist(Combn[i,])))
  primerOptionsReverse<-c(primerOptionsReverse,paste0(splitCharacters_copy,collapse=""))
}
# Get all combinations of possible forward and reverse amplification primers.
primerOptions<-vector(mode="list",length=2)
names(primerOptions)<-c("Forward","Reverse")
primerOptions[["Forward"]]<-primerOptionsForward
primerOptions[["Reverse"]]<-primerOptionsReverse
primerOptions<-expand.grid(primerOptions)
primerOptions$Expression<-paste0(primerOptions$Forward,"")

# Check how many foward and reverse matches there are.
NumFowardMatches<-c()
NumReverseMatches<-c()
for(i in 1:nrow(w)){
  w_seq<-w$Sequence[i]
  ForwardMatch<-c()
  for(j in primerOptionsForward){
    ForwardMatch<-c(ForwardMatch,ifelse(grepl(j,w_seq),T,F))
  }
  ReverseMatch<-c()
  for(j in primerOptionsReverse){
    ReverseMatch<-c(ReverseMatch,ifelse(grepl(j,w_seq),T,F))
  }
  NumFowardMatches<-c(NumFowardMatches,sum(ForwardMatch))
  NumReverseMatches<-c(NumReverseMatches,sum(ReverseMatch))
}

# There is only one forward and reverse match per sequence.
max(NumFowardMatches)
max(NumReverseMatches)

# Get just Woodhams sequences which have a forward and reverse match.
w_sub<-w[NumFowardMatches==1 & NumReverseMatches==1,]
nrow(w_sub) # There are 877 Woodhams sequences with forward and reverse matches.

# Trim sequences by forward and reverse primers.
TrimmedSeq<-c()
for(i in 1:nrow(w_sub)){
  w_seq<-w_sub$Sequence[i]
  for(j in primerOptionsForward){
    if(grepl(j,w_seq)){
      ForwardMatch<-j
    }
  }
  for(j in primerOptionsReverse){
    if(grepl(j,w_seq)){
      ReverseMatch<-j
    }
  }
  pattern<-paste0(ForwardMatch,"(.*?)",ReverseMatch)
  TrimmedSeq<-c(TrimmedSeq,regmatches(w_seq,regexec(pattern,w_seq))[[1]][2])
}

# Add trimmed sequences to the dataframe.
w_sub$TrimmedSequence<-TrimmedSeq

# How many sequences are not with in the 250 to 256 bp length?
sum(!(unname(sapply(w_sub$TrimmedSequence,FUN=nchar)) %in% 250:256)) # 5 sequences.

# Subset to just sequences within the 250 to 256 bp length range.
w_sub2<-w_sub[unname(sapply(w_sub$TrimmedSequence,FUN=nchar)) %in% 250:256,]

# Get summary of lengths.
summary(unname(sapply(w_sub2$TrimmedSequence,FUN=nchar)))

# Get just salamander sequences.
s<-subset(seq,Dataset=="Salamanders")

# Get just desired fields from Woodhams and salamander data.
## Woodhams.
w_sub3<-w_sub2[,c("Name","TrimmedSequence")]
colnames(w_sub3)<-c("Name","Sequence")
## Salamander.
s_2<-s[,c("Name","Sequence")]

# Combine Woodhams and salamander data.
ws<-rbind(w_sub3,s_2)

# Subset metadata to just select sequences.
meta2<-meta[meta$TaxaID %in% ws$Name,]

# Check that taxa IDs match across sequences and metadata.
identical(sort(meta2$TaxaID),sort(ws$Name))

# Write out sequences.
write.fasta(names=ws$Name,sequences=ws$Sequence,file="WoodhamsAndSalamandersUniqueSequences5.fasta")

# Write out metadata.
write.csv(x=meta2,file="WoodhamsBdInhibitionStatus5.csv",row.names=F)

# Check how many trimmed Woodhams sequences are unique.
length(unique(w_sub3$Sequence)) # Only 361 unique sequences (out of 872).

# Append field for antifungal status to Woodhams sequences.
w_sub3$Bd_Inhibition<-meta2$Bd_Inhibition[match(w_sub3$Name,meta2$TaxaID)]

# Check each unique Woodhams sequence for consistency in antifungal-ness.
## Loop through each unique sequence and get the number which are inhibitory and non-inhibitory.
num_seq_antifungal<-c()
num_seq_nonantifungal<-c()
for(i in unique(w_sub3$Sequence)){
  sub<-subset(w_sub3,Sequence==i)
  num_seq_antifungal<-c(num_seq_antifungal,sum(sub$Bd_Inhibition=="Inhibitory"))
  num_seq_nonantifungal<-c(num_seq_nonantifungal,sum(sub$Bd_Inhibition=="Non-Inhibitory"))
}
## Combine data into a data frame.
consistency<-data.frame(Sequence=unique(w_sub3$Sequence),NumberAntifungal=num_seq_antifungal,NumberNonAntifungal=num_seq_nonantifungal)
## Get just sequences which have more than one observation in the Woodhams data.
consistency<-consistency[rowSums(consistency[,c("NumberAntifungal","NumberNonAntifungal")])>1,]
## Calculate the proportion of antifungal traits in matching sequences.
consistency$ProportionAntifungal<-consistency$NumberAntifungal/(consistency$NumberAntifungal+consistency$NumberNonAntifungal)

# Create a histogram of the proportion of antifungal traits in matching Woodhams sequences.
#tiff("WoodhamsTraitConsistency.tiff",width=7,height=5,units="in",res=300)
hist(consistency$ProportionAntifungal,xlab="Proportion of Antifungal Traits in Matching Sequences",main="Consistency of Woodhams Antifungal Traits")
#dev.off()

# Summary of Woodhams sequence preparation:
## Removed sequences with unknown Bd inhibition status.
## Classified all sequences as inhibitory or non-inhibitory.
## Subsetted sequences to just Bacteria (using Woodhams Uclust taxonomy).
## Trimmed sequences to the region between our 16S amplification primers.
## Sequences without both forward and reverse amplification primers were discarded.
## Sequences that were outside the range of 250 to 256 bp were discarded.
## This leaves 872 Woodhams sequences available for alignment.
##
## There are 15,690 salamander sequences to be aligned.
## There are only 361 unique sequences in the Woodhams data (out of the 872).
## We cannot be certain that salamander sequences are antifungal with exact matching to the
## Woodhams data due to variable consistency in antifungal traits across the same Woodhams sequences.

# Done!
print("Done!")
