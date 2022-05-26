########################
### Create ASV Table ###
########################

# Provide paths to input and output files.
## Path to microbe functions script.
MicrobeFunctionsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/HelperFunctions/MicrobeFunctions.R"
## Path to forward reads.
ForwardReadsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/RegionITS/NextSeq/NoPrimers/NextSeq_ITS_NoPrimers-fwd.fastq"
## Path to reverse reads.
ReverseReadsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/RegionITS/NextSeq/NoPrimers/NextSeq_ITS_NoPrimers-rev.fastq"
## Path to reference taxonomies.
ReferenceTaxonomiesPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/RegionITS/NextSeq/ReferenceTaxonomies/NextSeqRefLibITS.csv"
## Path to output ASV table.
ASVTablePath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/RegionITS/NextSeq/ASVTable/ASV_Table_ITS.csv"

# Load microbe functions.
source(MicrobeFunctionsPath)

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# Report beginning of computations.
print("Starting script...")

# Read in NextSeq fastq files.
## Forward reads.
print(paste0("- Reading in forward fastq file: ",elapsed.time(ptm)))
fwd<-read.fastq(file=ForwardReadsPath)
## Reverse reads.
print(paste0("- Reading in reverse fastq file: ",elapsed.time(ptm)))
rev<-read.fastq(file=ReverseReadsPath)

# Create data frame with sample information and forward and reverse reads.
## Get sample information from the forward read data.
print(paste0("- Getting sample information: ",elapsed.time(ptm)))
samples<-data.frame(Sample=sapply(strsplit(x=fwd$Name,split="-"),"[[",1))
## Add in forward sequences.
print(paste0("- Adding foward sequences to full DF: ",elapsed.time(ptm)))
samples$ForwardRead<-fwd$Sequence
## Add in reverse sequences.
print(paste0("- Adding reverse sequences to full DF: ",elapsed.time(ptm)))
samples$ReverseRead<-rev$Sequence

# Clear forward and reverse sequence data frames from memory.
print(paste0("- Clearing foward and reverse DFs from memory: ",elapsed.time(ptm)))
rm(fwd); rm(rev)

# Concatenate forward and reverse NextSeq reads with a plus sign to get unique full length sequences.
print(paste0("- Concatenating forward and reverse NextSeq reads: ",elapsed.time(ptm)))
samples$FullRead<-paste0(samples$ForwardRead,"+",samples$ReverseRead)

# Read in taxonomic reference dataset.
print(paste0("- Reading in reference taxonomies: ",elapsed.time(ptm)))
reference_df<-read.csv(file=ReferenceTaxonomiesPath)

# Concatenate forward and reverse reference reads with a plus sign to get full length sequences.
print(paste0("- Concatenating forward and reverse reference reads: ",elapsed.time(ptm)))
reference_df$FullRead<-paste0(reference_df$ForwardRead,"+",reference_df$ReverseRead)

# Get number of total unique NextSeq sequences.
print(paste0("- Getting number of total unique NextSeq sequences: ",elapsed.time(ptm)))
number_of_total_unique_NextSeq_sequences<-length(unique(samples$FullRead))

# Get number of NextSeq reads.
print(paste0("- Getting number of NextSeq reads: ",elapsed.time(ptm)))
number_of_NextSeq_reads<-nrow(samples)

# Subset NextSeq reads to just those sequences which match the reference data.
print(paste0("- Subsetting to NextSeq reads with reference matches: ",elapsed.time(ptm)))
samples<-samples[samples$FullRead %in% reference_df$FullRead,]

# Create an ASV table from the sample information and full reads.
print(paste0("- Creating ASV table: ",elapsed.time(ptm)))
ASV_table<-as.data.frame.matrix(table(samples[,c("Sample","FullRead")]))

# Clear the samples data frame from memory.
print(paste0("- Clearing concatenated sequences DF from memory: ",elapsed.time(ptm)))
rm(samples)

# Extract full sequences from ASV table column names.
print(paste0("- Getting full sequences from ASV table: ",elapsed.time(ptm)))
sequences<-colnames(ASV_table)

# Get forward reads from full sequences.
print(paste0("- Extracting ASV table forward reads: ",elapsed.time(ptm)))
ForwardReads<-sapply(strsplit(x=sequences,split="\\+"),"[[",1)

# Get reverse reads from full sequences.
print(paste0("- Extracting ASV table reverse reads: ",elapsed.time(ptm)))
ReverseReads<-sapply(strsplit(x=sequences,split="\\+"),"[[",2)

# Clear full sequences from memory.
print(paste0("- Removing ASV table full sequences from memory: ",elapsed.time(ptm)))
rm(sequences)

# Assign taxonomy to ASV table columns.
print(paste0("- Assigning taxonomies to ASV table columns: ",elapsed.time(ptm)))
colnames(ASV_table)<-assign.taxa(forward_reads=ForwardReads,
                                 reverse_reads=ReverseReads,
                                 reference_df=reference_df)

# Clear reference taxonomies and forward and reverse reads from memory.
print(paste0("- Clearing reference taxonomies and ASV table forward and reverse reads from memory: ",elapsed.time(ptm)))
rm(reference_df); rm(ForwardReads); rm(ReverseReads)

# How many unique identified sequences are there?
print(paste0("- Getting number of unique identified sequences: ",elapsed.time(ptm)))
HowManyUniqueIdentifiedSequences<-ncol(ASV_table)

# What is the total number of reads which are identified?
print(paste0("- Getting total number of identified reads: ",elapsed.time(ptm)))
NumberOfIdentifiedReads<-sum(ASV_table)

# How many unique unidentified sequences are there?
print(paste0("- Getting number of unique unidentified sequences: ",elapsed.time(ptm)))
HowManyUniqueUnidentifiedSequences<-number_of_total_unique_NextSeq_sequences-HowManyUniqueIdentifiedSequences

# What is the total number of reads which are unidentified?
print(paste0("- Getting total number of unidentified reads: ",elapsed.time(ptm)))
NumberOfUnidentifiedReads<-number_of_NextSeq_reads-NumberOfIdentifiedReads

# Order ASV table columns.
print(paste0("- Reordering ASV columns: ",elapsed.time(ptm)))
ASV_table<-ASV_table[,order(colnames(ASV_table))]

# Order ASV table rows.
print(paste0("- Reordering ASV rows: ",elapsed.time(ptm)))
ASV_table<-ASV_table[order(row.names(ASV_table)),]

# Write out ASV table.
print(paste0("- Writing out ASV table: ",elapsed.time(ptm)))
write.csv(x=ASV_table,file=ASVTablePath,row.names=T)

# Report number of unidentied and identified sequences and reads.
print(paste0("- Reporting identification summary: ",elapsed.time(ptm)))
print(paste0("-- There were ",HowManyUniqueUnidentifiedSequences," unique unidentified sequences comprising ",NumberOfUnidentifiedReads," reads."))
print(paste0("-- There were ",HowManyUniqueIdentifiedSequences," unique identified sequences comprising ",NumberOfIdentifiedReads," reads."))

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
