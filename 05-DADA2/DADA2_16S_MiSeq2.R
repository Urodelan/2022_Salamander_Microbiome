#######################################
### Process MiSeq 16S Data in DADA2 ###
#######################################

# Clear lists and graphics.
rm(list=ls())
graphics.off()

# Set working directory.
setwd("/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/Region16S/MiSeq")

# Load packages.
#BiocManager::install("dada2",version="3.10")
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")

# Specify path to sequence data files.
path<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/Region16S/MiSeq/OrigSeqData/OrigNameRepair/NamesRepaired"
list.files(path) # Check that files are in path.

# Specify full paths to forward and reverse files. Sorting ensures matching when there are multiple files.
fnFs<-sort(list.files(path,pattern="ForwardNew.fastq",full.names=T))
fnRs<-sort(list.files(path,pattern="ReverseNew.fastq",full.names=T))

# Trim all reads to 290 bp with cutadapt.

## RUN CUTADAPT THROUGH R OR RSTUDIO ON THE TERMINAL PROPER (NOT RSTUDIO SERVER).

# Specify path to cutadapt.
cutadapt<-"/uufs/chpc.utah.edu/common/home/u0980895/.local/bin/cutadapt"

# Check that cutadapt works.
system2("/uufs/chpc.utah.edu/common/home/u0980895/.local/bin/cutadapt",args="--version")

# Run the following to perform cutadapt-ing!
path.SameLength <- file.path(path, "SameLength")
if(!dir.exists(path.SameLength)) dir.create(path.SameLength)
fnFs.SameLength <- file.path(path.SameLength, basename(fnFs))
fnRs.SameLength <- file.path(path.SameLength, basename(fnRs))

# Use cutadapt to trim all reads to 290 bp in both forward and reverse files.
# system2(cutadapt,args=c("-l 290", # New fixed read length.
#                         "-o",fnFs.SameLength[1],"-p",fnRs.SameLength[1], # Output files.
#                         fnFs[1],fnRs[1])) # Input files.

# Remove amplification primers with cutadapt.

# Specify amplication primers.
FWD<-"GTGYCAGCMGCCGCGGTAA" # Forward primer sequence.
REV<-"GGACTACHVGGGTWTCTAAT" # Reverse primer sequence.

# Define function to find all primer orientations.
allOrients<-function(primer){
  # Create all orientations of the input sequence.
  require(Biostrings)
  dna<-DNAString(primer) # The Biostrings works with DNAString objects rather than character vectors.
  orients<-c(Forward=dna,Complement=complement(dna),Reverse=reverse(dna),
             RevComp=reverseComplement(dna))
  return(sapply(orients,toString)) # Convert back to character vector.
}

# Get all orientations of forward and reverse amplification primers.
(FWD.orients<-allOrients(FWD))
(REV.orients<-allOrients(REV))

# Remove reads with Ns using cutadapt.

# Put N-filterd files in filtN subdirectory.
path.filtN <- file.path(path, "filtN")
if(!dir.exists(path.filtN)) dir.create(path.filtN)
fnFs.filtN<-file.path(path,"filtN",basename(fnFs))
fnRs.filtN<-file.path(path,"filtN",basename(fnRs))

# Use cutadapt.
# system2(cutadapt,args=c("--max-n 0", # Remove read pairs which contain Ns.
#                         "-o",fnFs.filtN[1],"-p",fnRs.filtN, # Output files.
#                         fnFs.SameLength[1],fnRs.SameLength[1])) # Input files.

# Create function to identify primers in sequence reads.
primerHits<-function(primer,fn){
  # Counts number of reads in which the primer is found.
  nhits<-vcountPattern(primer,sread(readFastq(fn)),fixed=F)
  return(sum(nhits>0))
}

# Check first sample for primer hits.
# rbind(FWD.ForwardReads=sapply(FWD.orients,primerHits,fn=fnFs.filtN[[1]]),
#       FWD.ReverseReads=sapply(FWD.orients,primerHits,fn=fnRs.filtN[[1]]),
#       REV.ForwardReads=sapply(REV.orients,primerHits,fn=fnFs.filtN[[1]]),
#       REV.ReverseReads=sapply(REV.orients,primerHits,fn=fnRs.filtN[[1]]))

# Run the following to perform cutadapt-ing!
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Forward reads should have the forward amplification primer minus the first character, and possibly
# a reverse complement of the full reverse amplification primer.
(R1.flags<-paste0('-a "^',gsub('^.','',FWD),'...',reverseComplement(DNAString(REV)),';required"'))

# Reverse reads should have the reverse amplification primer minus the first character, and possibly
# a reverse complement of the full forward amplification primer.
(R2.flags<-paste0('-A "^',gsub('^.','',REV),'...',reverseComplement(DNAString(FWD)),';required"'))

# Run cutadapt
# for(i in seq_along(fnFs)){
#   system2(cutadapt,args=c(R1.flags,R2.flags,"-e 0.15","--discard-untrimmed",
#                           "--minimum-length 1", # Keep only reads at least 1 bp long.
#                           "-o",fnFs.cut[i],"-p",fnRs.cut[i], # output files
#                           fnFs.filtN[i],fnRs.filtN[i])) # input files
# }

# Check that primers have been removed.
# rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
#       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
#       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
#       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "ForwardNew.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "ReverseNew.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
## Modified function from online tutorial.
get.sample.name <- function(fname) substr(basename(fname),start=1,stop=8)
# strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# Inspect read quality profiles.
# plotQualityProfile(cutFs[1]) # Forward.
# plotQualityProfile(cutRs[1]) # Reverse.

# Filter and trim.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
# out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
#                      truncQ = 2, minLen = 10, rm.phix = TRUE, compress = FALSE, multithread = TRUE)
# save(out,file="out.RData") # Save out object.
# head(out)
# 
# # Learn the error rates.
# errF <- learnErrors(filtFs, multithread = TRUE)
# save(errF,file="errF.RData") # Save forward error rate object.
# errR <- learnErrors(filtRs, multithread = TRUE)
# save(errR,file="errR.RData") # Save reverse error rate object.
# 
# # Visualize error rates.
# # plotErrors(errF, nominalQ = TRUE)
# 
# # Sample inference algorithm.
# dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
# save(dadaFs,file="dadaFs.RData") # Save forward dada2 object.
# dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
# save(dadaRs,file="dadaRs.RData") # Save reverse dada2 object.

###

# Load out object.
load("out.RData")

# Load error rate objects.
load("errF.RData")
load("errR.RData")

# Load dada2 objects,
load("dadaFs.RData")
load("dadaRs.RData")

# Inspect dada2 object.
head(dadaFs[[1]])

# # Merge pairs.
# merger<-mergePairs(dadaFs,filtFs,dadaRs,filtRs,minOverlap=12,maxMismatch=0,verbose=T)
# 
# # Format merger as a proper data frame.
# merger2<-as.data.frame(do.call(cbind,merger))
# 
# # Format merger2 data frame.
# merger2$sequence<-as.character(merger2$sequence)
# merger2$abundance<-as.integer(merger2$abundance)
# merger2$forward<-as.integer(merger2$forward)
# merger2$reverse<-as.integer(merger2$reverse)
# merger2$nmatch<-as.integer(merger2$nmatch)
# merger2$nmismatch<-as.integer(merger2$nmismatch)
# merger2$nindel<-as.integer(merger2$nindel)
# merger2$prefer<-as.numeric(merger2$prefer)
# merger2$accept<-as.logical(merger2$accept)
# 
# # Add merger2 to the first element of a list in merger3.
# merger3<-list(merger2)
# 
# # Inspect the merger3 data.frame from the first sample.
# head(merger3[[1]])
# 
# # Save merger3 list.
# save(merger3,file="merger3.RData")

# Load merger3 list.
load("merger3.RData")

# # Construct sequence table.
# seqtab<-makeSequenceTable(merger3)
# dim(seqtab)
# 
# # Inspect distribution of sequence lengths.
# table(nchar(getSequences(seqtab)))
# 
# ### Synthgenes from the protocol -
# ## ITS: TGCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT
# ## 16S: GCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT
# 
# # Create vector of non-biological sequences.
# synthgene<-"GCCACAGATACGTACCGCTCATAACGCGAACCGAAGCGCAGTAGAAGTACTCCGTATCCTACCTCGGTCGTGGTTTAGGCTATCGACATCTTGCATGGGCTTCCCTAGTGAACTCTTGGGATGT"
# coligos<-read.csv("/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/Region16S/Coligos/Coligos16S.csv")
# nonbio_seq<-c(synthgene,coligos$SequenceNoAmpPrimer)
# 
# # Get proportion of non-biological sequences which are included in the denoised sequences.
# sum(nonbio_seq %in% colnames(seqtab))/length(nonbio_seq) # 94.8%
# 
# # Which non-biological sequences are not included in the denoised sequences.
# ## Synthgene is included.
# synthgene %in% colnames(seqtab)
# ## 5 coligos are missing from the denoised sequences.
# sum(!(coligos$SequenceNoAmpPrimer %in% colnames(seqtab)))
# ## Which coligos are missing from the denoised sequences.
# ## 5 coligos are missing from the denoised sequences.
# coligos$SequenceNoAmpPrimer[!(coligos$SequenceNoAmpPrimer %in% colnames(seqtab))]
# 
# # Pull out non-biological sequences (synthgene & coligos).
# seqtab.nonbio_which_keep<-colnames(seqtab) %in% nonbio_seq
# ## Ackwardly format a matrix for seqtab.nonbio.
# seqtab.nonbio_seq<-colnames(seqtab)[seqtab.nonbio_which_keep]
# seqtab.nonbio_counts<-seqtab[1,seqtab.nonbio_which_keep]
# seqtab.nonbio<-matrix(seqtab.nonbio_counts,nrow=1)
# colnames(seqtab.nonbio)<-seqtab.nonbio_seq
# 
# # Inspect dimensions of seqtab.nonbio.
# dim(seqtab.nonbio)
# 
# # Inspect distribution of non-biological sequence lengths.
# table(nchar(getSequences(seqtab.nonbio)))
# 
# # Pull out biological sequences (removing synthgene & coligos) which are between 250 and 256 bp.
# seqtab.bio_which_keep<-(!(colnames(seqtab) %in% nonbio_seq)) & (nchar(colnames(seqtab)) %in% 250:256)
# ## Ackwardly format a matrix for seqtab.bio.
# seqtab.bio_seq<-colnames(seqtab)[seqtab.bio_which_keep]
# seqtab.bio_counts<-seqtab[1,seqtab.bio_which_keep]
# seqtab.bio<-matrix(seqtab.bio_counts,nrow=1)
# colnames(seqtab.bio)<-seqtab.bio_seq
# 
# # Inspect dimensions of seqtab.bio.
# dim(seqtab.bio)
# 
# # Inspect distribution of biological sequence lengths.
# table(nchar(getSequences(seqtab.bio)))
# 
# # Combine biological and non-biological sequences back together.
# seqtab.LegitimateSeq<-cbind(seqtab.nonbio,seqtab.bio)
# dim(seqtab.LegitimateSeq)
# table(nchar(getSequences(seqtab.LegitimateSeq)))
# 
# # Remove chimeras.
# seqtab.nochim<-removeBimeraDenovo(seqtab.LegitimateSeq,method="consensus",multithread=T,verbose=T)
# dim(seqtab.nochim)
# 
# # Check proportion of reads not removed by chimera removal.
# sum(seqtab.nochim)/sum(seqtab)
# 
# # Re-inspect distribution of sequence lengths.
# table(nchar(getSequences(seqtab.nochim)))
# 
# # Remove non-biological sequences again.
# seqtab.nochimbio_which_keep<-(!(colnames(seqtab.nochim) %in% nonbio_seq)) & (nchar(colnames(seqtab.nochim)) %in% 250:256)
# ## Ackwardly format a matrix for seqtab.nochimbio.
# seqtab.nochimbio_seq<-colnames(seqtab.nochim)[seqtab.nochimbio_which_keep]
# seqtab.nochimbio_counts<-seqtab.nochim[1,seqtab.nochimbio_which_keep]
# seqtab.nochimbio<-matrix(seqtab.nochimbio_counts,nrow=1)
# colnames(seqtab.nochimbio)<-seqtab.nochimbio_seq
# 
# # Inspect dimensions of seqtab.nochimbio.
# dim(seqtab.nochimbio)
# 
# # Inspect distribution of biological sequence lengths.
# table(nchar(getSequences(seqtab.nochimbio)))
# 
# # Add back in non-biological sequences from before chimera removal, since many coligos were removed
# # during chimera removal.
# seqtab.nochimfull<-cbind(seqtab.nonbio,seqtab.nochimbio)
# dim(seqtab.nochimfull)
# table(nchar(getSequences(seqtab.nochimfull)))
# 
# # Check how many duplicate sequences are in the full ASV/non-bio table. This should be zero.
# sum(duplicated(colnames(seqtab.nochimfull)))
# 
# # Check that the synthgene is included in the full ASV/non-bio table.
# synthgene %in% colnames(seqtab.nochimfull)
# 
# # Check how many sequence reads are of the synthgene. (395,646 reads)
# sum(seqtab.nochimfull[,colnames(seqtab.nochimfull)==synthgene])
# 
# # Check what proportion of sequence reads are of the synthgene. (7.4%)
# paste0(round(sum(seqtab.nochimfull[,colnames(seqtab.nochimfull)==synthgene])/sum(seqtab.nochimfull),3)*100,"%")
# 
# # Save no-chimera full sequence table.
# save(seqtab.nochimfull,file="seqtab.nochimfull.RData")
# 
# # Save no-chimera biological sequence table.
# save(seqtab.nochimbio,file="seqtab.nochimbio.RData")
# 
# # Save non-biological sequence table.
# save(seqtab.nonbio,file="seqtab.nonbio.RData")

# Load no-chimera full sequence table.
load("seqtab.nochimfull.RData")

# Load no-chimera  biological sequence table.
load("seqtab.nochimbio.RData")

# Load non-biological sequence table.
load("seqtab.nonbio.RData")

# Track reads through the pipeline.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(merger3[[1]]), rowSums(seqtab.nochimfull), rowSums(seqtab.nochimbio), rowSums(seqtab.nonbio))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimfull", "nochimbio", "nonbio")
rownames(track) <- sample.names
head(track)

# Save track.
save(track,file="track.RData")

# Assign taxonomy.
taxa<-assignTaxonomy(seqtab.nochimbio,"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/Region16S/MiSeq/ReferenceDBs/Silva_Full_DB.fa",multithread=T)

# Add species.
taxa2<-addSpecies(taxa,"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/DADA2/Region16S/MiSeq/ReferenceDBs/Silva_AddSpecies_DB.fa")

# Print taxa.
taxa.print <- taxa2 # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print,n=20)

# Save taxonomy.
save(taxa2,file="taxa.RData")

# Finish formatting custom reference data base.

# Copy taxa.
taxa_reformat<-as.data.frame(taxa2)

# Reformat taxa names.
formatted_taxa_names<-c()
for(i in 1:nrow(taxa_reformat)){
  taxa_row<-unname(unlist(taxa_reformat[i,]))
  taxa_row<-taxa_row[!is.na(taxa_row)]
  taxa_row_reformat<-base::paste0(base::paste(taxa_row,collapse=";"),";")
  formatted_taxa_names<-c(formatted_taxa_names,taxa_row_reformat)
}

# Get taxa sequences.
taxa_sequences<-row.names(taxa2)

# Write out custom 16S library.
library(seqinr)
write.fasta(sequences=as.list(taxa_sequences),
            names=formatted_taxa_names,
            file.out="Generated_Library/Custom_16S_Library.fasta")

# Done!
