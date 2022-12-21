# Remake custom reference libraries for use with NextSeq data.

# Clear environment.
rm(list=ls())

# Read in helper functions.
source("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/HelperFunctions/MicrobeFunctions.R")

# Define 16S amplification primers.
Forward_amp_primer16S<-"GTGYCAGCMGCCGCGGTAA"
Reverse_amp_primer16S<-"GGACTACHVGGGTWTCTAAT"

# Read in 16S reference sequences.
df<-read.fasta("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/ReferenceLibraries/Custom_16S_Library.fasta")

# Get just bacterial reference sequences.
df<-df[grepl("^Bacteria;",df$Name),]

# Read in 16S coligos.
coligos_16S<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Lab/LibraryPrep/Coligos16S.csv")
coligos_16S$Name<-paste0("Coligo16S_",coligos_16S$Position)
length(unique(coligos_16S$Name))
coligos_16S<-coligos_16S[,c("Name","SequenceNoAmpPrimers")]
colnames(coligos_16S)[2]<-"Sequence"

# Combine reference and coligo sequences.
df<-rbind(df,coligos_16S)

# Get 16S synthgene sequence.
synthgene16S<-"GCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT"

# Add in 16S synthgene.
df<-rbind(df,data.frame(Name="Synthgene16S",Sequence=synthgene16S))

# All NextSeq reads will be trimmed to 140 bp length. Then amplification primers will be removed,
# except for the first base, which was already removed during demultiplexing.

# So truncate reads to 140 - (length of 5' amplification primer - 1).
library(Biostrings)
## Truncate forward sequences.
df$ForwardSequence<-substr(df$Sequence,start=1,stop=140-(nchar(Forward_amp_primer16S)-1))
## Truncate reverse sequences. Derive the reverse compliment from the sequence such that
## the reverse read appears as it would in the NextSeq data.
df$ReverseSequence<-NA
for(i in 1:nrow(df)){
  df$ReverseSequence[i]<-base::substr(reverseComplement(DNAString(df$Sequence[i])),start=1,
                                      stop=140-(nchar(Reverse_amp_primer16S)-1))
}

# Check if there are any duplicates after truncation.
## Concatenate forward and reverse reads to get a unique ID of sorts.
df$ForwardSequenceAndReverseSequence<-paste0(df$ForwardSequence,"+",df$ReverseSequence)
paste0(length(unique(df$ForwardSequenceAndReverseSequence))," of ",nrow(df)," sequences are still unique.")
sum(duplicated(df$ForwardSequenceAndReverseSequence)) # 52 duplicate sequences after truncation.
sum(duplicated(df$Name)) # 15390 duplicate names.

# Get duplicated read combos (from the concatenated forward and reverse reads).
duplicated_read_combos<-unique(df$ForwardSequenceAndReverseSequence[duplicated(df$ForwardSequenceAndReverseSequence)])

# Separate duplicated and not duplicated sequences into different data frames.
df_duplicated<-df[df$ForwardSequenceAndReverseSequence %in% duplicated_read_combos,]
df_notduplicated<-df[!(df$ForwardSequenceAndReverseSequence %in% duplicated_read_combos),]

# Define IUPAC wildcard characters.
IUPAC_wildcards<-list(R=c("A","G"),Y=c("C","T"),S=c("G","C"),W=c("A","T"),K=c("G","T"),M=c("A","C"),B=c("C","G","T"),D=c("A","G","T"),H=c("A","C","T"),V=c("A","C","G"),N=c("A","T","G","C"))

# Condense duplicated sequences (truncated reads) to just consensus names and consensus full sequences.
## Create storage data frame.
df_duplicated_condensed<-data.frame(NULL)
## Create storage vectors for monitoring sequence set lengths and number of mismatches.
are_same_lengths<-c()
num_mismatches<-c()
## Loop through each unique combination of truncated forward and reverse reads.
for(i in unique(df_duplicated$ForwardSequenceAndReverseSequence)){
  ## Subset to just a single combination of truncated forward and reverse reads.
  sub<-df_duplicated[df_duplicated$ForwardSequenceAndReverseSequence==i,]
  ## Split up taxonomies into individual levels of classification.
  sub_names<-base::strsplit(x=sub$Name,split=";")
  ## What is the lowest level of classification of any of the taxonomies?
  min_len_sub_names<-min(sapply(X=sub_names,FUN=length))
  
  ## Create new consensus names.
  ### Create storage vector for whether each level of classification is consistent
  ### across sequences in the set.
  all_sub_names2_equal<-c()
  ### Loop through all taxonomic levels (up to the lowest level of any sequence in the set).
  for(j in 1:min_len_sub_names){
    ### Get the names of a single taxonomic level.
    sub_names2<-sapply(X=sub_names,FUN="[[",j)
    ### Check to see whether all names in this single taxonomic level are the same.
    all_sub_names2_equal<-c(all_sub_names2_equal,all(sub_names2==sub_names2[1]))
  }
  ### Get the index of the highest agreed-upon taxonomic level.
  max_of_same<-max(which(all_sub_names2_equal==T))
  ### Extract values of the highest agreed-upon taxonomy from the first sequence in the set.
  vector_of_name_elements<-sub_names[[1]][1:max_of_same]
  ### Collapse the highest agreed-upon taxonomy from a vector into a full taxonomic names
  ### with levels separated by a semi-colon.
  full_name<-paste0(paste0(vector_of_name_elements,collapse=";"),";")
  
  # For 16S, let's try and get a good idea of a consensus sequence.
  # This can be used in the phylogenetic tree and character mapping.
  
  # Monitor whether all sequences in the set are the same length.
  are_same_lengths<-c(are_same_lengths,length(unique(nchar(sub$Sequence)))==1)
  
  ## Create new consensus sequences.
  ### Split up all sequences in the set into individual nucleotides.
  splitCharacters<-base::strsplit(x=sub$Sequence,split="")
  ### Create storage vector checking to see which nucleotides match across the sequences in the set.
  match_across_seq<-c()
  ### Loop through each nucleotide position.
  for(k in 1:length(splitCharacters[[1]])){
    ### Get the set of nucleotide characters at the position.
    characters_at_position<-sapply(X=splitCharacters,FUN="[[",k)
    ### Check whether all nucleotide characters at the position are the same.
    match_across_seq<-c(match_across_seq,all(characters_at_position==characters_at_position[1]))
  }
  ### Get the positions at which nucleotides do not match across the sequences in the set.
  mismatch_location<-which(match_across_seq==F)
  ### Store the number of mismatches within each set of sequences.
  num_mismatches<-c(num_mismatches,length(mismatch_location))
  ### Create an empty storage vector for storing which wildcard characters map to the mismatches.
  wildcards<-c()
  ### Loop through each mismatch position.
  for(k in mismatch_location){
    ### Get the set of mismatching nucleotides.
    mismatching_nts<-sapply(X=splitCharacters,FUN="[[",k)
    ### Create storage vector indicating whether a given IUPAC wildcard includes all the
    ### nucleotides in the mismatching set.
    all_nts_included<-c()
    ### Loop through each vector of nucleotides representing each IUPAC wildcard.
    for(h in IUPAC_wildcards){
      ### Check to see if all nucleotides in the mismatching set are represented in an IUPAC wildcard.
      all_nts_included<-c(all_nts_included,all(mismatching_nts %in% h))
    }
    ### Get the index of the simplest wildcard which encompasses all observed nucleotides
    ### in the mismatching set.
    simplest_wildcard_index<-min(which(all_nts_included==T))
    ### Get the simplest wildcard character which encompasses all observed nucleotides in the
    ### mismatching set.
    simplest_wildcard_character<-names(IUPAC_wildcards)[simplest_wildcard_index]
    ### Store the wildcard characters for the sequence in a storage vector.
    wildcards<-c(wildcards,simplest_wildcard_character)
  }
  ### Copy the first sequence in the set as the new sequence.
  newSequence<-splitCharacters[[1]]
  ### At the mismatch locations, replace the nucleotides with their appropriate wildcard characters.
  newSequence[mismatch_location]<-wildcards
  ### Collapse the new sequence from a vector of individual characters to a full names string.
  newSequence<-paste0(newSequence,collapse="")
  
  # Create data frame for the sequence to bind to the storage data frame.
  df_tobind<-data.frame(Name=full_name, # Use consensus taxonomy.
                        Sequence=newSequence, # Use consensus sequence.
                        ForwardSequence=sub$ForwardSequence[1], # Get truncated forward sequence.
                        ReverseSequence=sub$ReverseSequence[1], # Get truncated reverse sequence.
                        # Get concatenated truncated forward and reverse sequence.
                        ForwardSequenceAndReverseSequence=sub$ForwardSequenceAndReverseSequence[1],
                        stringsAsFactors=F)
  # Append sequence set information to storage data frame.
  df_duplicated_condensed<-rbind(df_duplicated_condensed,df_tobind)
}

# There is at least one wildcard character in each consensus sequence.
checking_for_wildcards<-base::strsplit(x=df_duplicated_condensed$Sequence,split="")
all(sapply(X=checking_for_wildcards,FUN=function(x) sum(names(IUPAC_wildcards) %in% x)>0))

# All duplicate NextSeq-observable forward and reverse reads belonged to sequences of the same length.
all(are_same_lengths)

# How many mismatches were observed in each set of sequences?
table(num_mismatches)
## 1 mismatch: 49 sequences
## 2 mismatches: 3 sequences

# Combine non-duplicated and condensed duplicated data frames.
df_full16S<-rbind(df_notduplicated,df_duplicated_condensed)
sum(duplicated(df_full16S$ForwardSequenceAndReverseSequence)) # All sequences unique.

# Pull out coligos and synthgene.
df_full16S_colAndSyn<-df_full16S[grepl("^Coligo16S_",df_full16S$Name) |
                                   df_full16S$Name=="Synthgene16S",]
df_full16S<-df_full16S[!(grepl("^Coligo16S_",df_full16S$Name) |
                           df_full16S$Name=="Synthgene16S"),]

# Add a number to each of our sequence names so that each is unique.
## Create an empty column.
df_full16S$NewName<-NA
## Loop through each unique sequence name.
for(i in unique(df_full16S$Name)){
  ## Get all records with this particular unique sequence name.
  sub<-subset(df_full16S,Name==i)
  ## Append a number such that all sequence names are unique.
  df_full16S$NewName[df_full16S$Name==i]<-paste0(sub$Name,"_",1:nrow(sub))
}

# Use original names for coligos and synthgene.
df_full16S_colAndSyn$NewName<-df_full16S_colAndSyn$Name

# Add back in coligos and synthgene.
df_full16S<-rbind(df_full16S,df_full16S_colAndSyn)

# Subset columns.
colnames(df_full16S)
ForUseWithWoodhams<-df_full16S[,c("NewName","Sequence")]
colnames(ForUseWithWoodhams)[1]<-"Name"
df_full16S<-df_full16S[,c("NewName","ForwardSequence","ReverseSequence")]
colnames(df_full16S)<-c("Name","ForwardRead","ReverseRead")

# Remove coligos and synthgene from ForUseWithWoodhams.
ForUseWithWoodhams<-ForUseWithWoodhams[!(grepl("^Coligo16S_",ForUseWithWoodhams$Name)
                                         | ForUseWithWoodhams$Name=="Synthgene16S"),]

# What are the read length counts?
table(nchar(df_full16S$ForwardRead))
table(nchar(df_full16S$ReverseRead))
table(nchar(ForUseWithWoodhams$Sequence))

# No duplicate sequences.
sum(duplicated(paste0(df_full16S$ForwardRead,"...",df_full16S$ReverseRead)))
sum(duplicated(ForUseWithWoodhams$Sequence))

# No duplicate names.
sum(duplicated(df_full16S$Name))
sum(duplicated(ForUseWithWoodhams$Name))

# Write out truncated reference file (16S).
write.csv(df_full16S,file="TruncatedRefLibs/NextSeqRefLib16S.csv",row.names=F)

# Write out file for use with Woodhams.
write.csv(ForUseWithWoodhams,file="Reference16SForUseWithWoodhams.csv",row.names=F)

# Clear environment.
rm(list=ls())

# Do same for ITS reference library.

# Read in helper functions.
source("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/HelperFunctions/MicrobeFunctions.R")

# Define amplification primers.
Forward_amp_primerITS<-"CTTGGTCATTTAGAGGAAGTAA"
Reverse_amp_primerITS<-"GCTGCGTTCTTCATCGATGC"

# Read in ITS reference sequences.
df<-read.fasta("/Users/kenengoodwin/Desktop/Microbiome/DataProcessing/ClusterScripts/ReferenceLibraries/Custom_ITS_Library.fasta")

# Get just fungi reference sequences.
df<-df[grepl("^Fungi;",df$Name),]

# Read in ITS coligos.
coligos_ITS<-read.csv("/Users/kenengoodwin/Desktop/Microbiome/Lab/LibraryPrep/ColigosITS.csv")
coligos_ITS$Name<-paste0("ColigoITS_",coligos_ITS$Position)
length(unique(coligos_ITS$Name))
coligos_ITS<-coligos_ITS[,c("Name","SequenceNoAmpPrimer")]
colnames(coligos_ITS)[2]<-"Sequence"

# Combine reference and coligo sequences.
df<-rbind(df,coligos_ITS)

# Get ITS synthgene sequence.
synthgeneITS<-"TGCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT"

# Add in ITS synthgene.
df<-rbind(df,data.frame(Name="SynthgeneITS",Sequence=synthgeneITS))

# All NextSeq reads will be trimmed to 140 bp length - (length of 5' amplification primer - 1).
library(Biostrings)
## Truncate forward sequences.
df$ForwardSequence<-base::substr(
  sapply(base::strsplit(x=df$Sequence,split="NNNNNNNNNN"),"[[",1),start=1,
  stop=140-(nchar(Forward_amp_primerITS)-1))
## Truncate reverse sequences.
df$ReverseSequence<-NA
for(i in 1:nrow(df)){
  if(grepl("NNNNNNNNNN",df$Sequence[i])){
    df$ReverseSequence[i]<-base::substr(
      reverseComplement(DNAString(
        base::strsplit(x=df$Sequence[i],split="NNNNNNNNNN")[[1]][2])),start=1,
      stop=140-(nchar(Reverse_amp_primerITS)-1))
    } else {
      df$ReverseSequence[i]<-base::substr(reverseComplement(DNAString(df$Sequence[i])),start=1,
                                          stop=140-(nchar(Reverse_amp_primerITS)-1))
  }
}

# Are there any Ns remaining?
sum(grepl("N",df$ForwardSequence)) # 0
sum(grepl("N",df$ReverseSequence)) # 0

# Check if there are any duplicates after truncation.
df$ForwardSequenceAndReverseSequence<-paste0(df$ForwardSequence,"+",df$ReverseSequence)
paste0(length(unique(df$ForwardSequenceAndReverseSequence))," of ",nrow(df)," sequences are still unique.")
sum(duplicated(df$ForwardSequenceAndReverseSequence)) # 46 duplicate sequences after truncation.
sum(duplicated(df$Name)) # 3438 duplicate names.

# Get duplicated read combos.
duplicated_read_combos<-unique(df$ForwardSequenceAndReverseSequence[duplicated(df$ForwardSequenceAndReverseSequence)])

# How many Bd variants total?
sum(grepl(";Batrachochytrium;dendrobatidis;",df$Name)) # 25!!!

# Separate duplicated and not duplicated sequences into different data frames.
df_duplicated<-df[df$ForwardSequenceAndReverseSequence %in% duplicated_read_combos,]
df_notduplicated<-df[!(df$ForwardSequenceAndReverseSequence %in% duplicated_read_combos),]

# Condense duplicated sequences to just consensus names.
df_duplicated_condensed<-data.frame(NULL)
for(i in unique(df_duplicated$ForwardSequenceAndReverseSequence)){
  sub<-df_duplicated[df_duplicated$ForwardSequenceAndReverseSequence==i,]
  sub_names<-base::strsplit(x=sub$Name,split=";")
  min_len_sub_names<-min(sapply(X=sub_names,FUN=length))
  all_sub_names2_equal<-c()
  for(j in 1:min_len_sub_names){
    sub_names2<-sapply(X=sub_names,FUN="[[",j)
    all_sub_names2_equal<-c(all_sub_names2_equal,all(sub_names2==sub_names2[1]))
  }
  max_of_same<-max(which(all_sub_names2_equal==T))
  vector_of_name_elements<-sub_names[[1]][1:max_of_same]
  full_name<-paste0(paste0(vector_of_name_elements,collapse=";"),";")
  df_tobind<-data.frame(Name=full_name,Sequence=sub$Sequence[1],
                        ForwardSequence=sub$ForwardSequence[1],
                        ReverseSequence=sub$ReverseSequence[1],
                        ForwardSequenceAndReverseSequence=sub$ForwardSequenceAndReverseSequence[1])
  df_duplicated_condensed<-rbind(df_duplicated_condensed,df_tobind)
}

# Combine non-duplicated and condensed duplicated data frames.
df_fullITS<-rbind(df_notduplicated,df_duplicated_condensed)
sum(duplicated(df_fullITS$ForwardSequenceAndReverseSequence)) # All sequences unique.

# Pull out coligos and synthgene.
df_fullITS_colAndSyn<-df_fullITS[grepl("^ColigoITS_",df_fullITS$Name) |
                                   df_fullITS$Name=="SynthgeneITS",]
df_fullITS<-df_fullITS[!(grepl("^ColigoITS_",df_fullITS$Name) |
                           df_fullITS$Name=="SynthgeneITS"),]

# Add a number to each of our sequence names so that each is unique.
## Create an empty column.
df_fullITS$NewName<-NA
## Loop through each unique sequence name.
for(i in unique(df_fullITS$Name)){
  ## Get all records with this particular unique sequence name.
  sub<-subset(df_fullITS,Name==i)
  ## Append a number such that all sequence names are unique.
  df_fullITS$NewName[df_fullITS$Name==i]<-paste0(sub$Name,"_",1:nrow(sub))
}

# Use original names for coligos and synthgene.
df_fullITS_colAndSyn$NewName<-df_fullITS_colAndSyn$Name

# Add back in coligos and synthgene.
df_fullITS<-rbind(df_fullITS,df_fullITS_colAndSyn)

# Subset columns.
colnames(df_fullITS)
df_fullITS<-df_fullITS[,c("NewName","ForwardSequence","ReverseSequence")]
colnames(df_fullITS)<-c("Name","ForwardRead","ReverseRead")

# What are the read length counts?
table(nchar(df_fullITS$ForwardRead))
table(nchar(df_fullITS$ReverseRead))

# No duplicate sequences.
sum(duplicated(paste0(df_fullITS$ForwardRead,"...",df_fullITS$ReverseRead)))

# No duplicate names.
sum(duplicated(df_fullITS$Name))

# Write out truncated reference file (ITS).
write.csv(df_fullITS,file="TruncatedRefLibs/NextSeqRefLibITS.csv",row.names=F)

# Done!
print("Done!")
