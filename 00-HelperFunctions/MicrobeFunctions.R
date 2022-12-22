### These are functions to assist in working with the amphibian microbe data.
### Use the command 'source("MicrobeFunctions.R")' to load these functions.

# Define the simplex function.
## Arguments:
### x: A vector of values.
## Returns: A vector of proportions (i.e., x/sum(x)).
simplex<-function(x) x/sum(x)

# Define the softmax function.
## Arguments:
### x: A vector of values.
## Returns: A vector containing the output of the softmax function.
softmax<-function(x) exp(x)/sum(exp(x))

# Define function for reading in taxa tables.
## Arguments:
### file: Pathname to a taxa table file to read.
## Returns: A data frame with sample ID as row names and taxa names as column names.
read.taxa.table<-function(file){
  # Read in taxa table.
  reads<-read.csv(file,row.names=1,check.names=F)
  # Return the data frame.
  return(reads)
}

# Define function for reading in fasta files.
## Arguments:
### file: Pathname to a fasta file to read.
## Returns: A data frame with fields for sequence names and sequences.
read.fasta<-function(file){
  # Read in fasta file.
  fasta<-read.delim(file=file,header=F,stringsAsFactors=F)
  # Get sequence names.
  names<-gsub(pattern="^>",replacement="",x=fasta$V1[seq(from=1,to=nrow(fasta),by=2)])
  # Get sequences.
  sequences<-fasta$V1[seq(from=2,to=nrow(fasta),by=2)]
  # Create data frame with names and sequences.
  df<-data.frame(Name=names,Sequence=sequences,stringsAsFactors=F)
  # Return the data frame.
  return(df)
}

# Define function for writing fasta files.
## Arguments:
### names: A vector of sequence names.
### sequences: A vector of sequences.
### file: Pathname to a fasta file to write.
write.fasta<-function(names,sequences,file){
  # Check that names and sequences are of the same length.
  if(length(names)!=length(sequences)) stop("Names and sequences are not of the same length.")
  # Append '>' to the start of sequence names.
  names<-paste0(">",names)
  # Create an empty data frame to store what will be written to a text file.
  df<-as.data.frame(matrix(NA,nrow=2*length(names),ncol=1))
  # Add names to the text file data frame.
  df$V1[seq(from=1,to=nrow(df),by=2)]<-names
  # Add sequences to the text file data frame.
  df$V1[seq(from=2,to=nrow(df),by=2)]<-as.character(sequences)
  # Write out fasta file.
  write.table(x=df,file=file,row.names=F,col.names=F,quote=F)
}

# Define function for reading in fastq files.
## Arguments:
### file: Pathname to a fastq file to read.
## Returns: A data frame with fields for sequence names, sequences, and quality scores.
read.fastq<-function(file){
  # Read in fastq file.
  fastq<-read.delim(file=file,header=F,stringsAsFactors=F)
  # Get sequence names.
  names<-gsub(pattern="^@",replacement="",x=fastq$V1[seq(from=1,to=nrow(fastq),by=4)])
  # Get sequences.
  sequences<-fastq$V1[seq(from=2,to=nrow(fastq),by=4)]
  # Get quality scores.
  Qscores<-fastq$V1[seq(from=4,to=nrow(fastq),by=4)]
  # Create data frame with names, sequences, and quality scores.
  df<-data.frame(Name=names,Sequence=sequences,Qscore=Qscores,stringsAsFactors=F)
  # Return the data frame.
  return(df)
}

# Define function for writing fastq files.
## Arguments:
### names: A vector of sequence names.
### sequences: A vector of sequences.
### Qscores: A vector of quality scores.
### file: Pathname to a fastq file to write.
write.fastq<-function(names,sequences,Qscores,file){
  # Check that names, sequences, and quality scores are of the same length.
  if(!all(sapply(list(length(names),length(sequences),length(Qscores)),
                 function(x) x==length(Qscores)))){
    stop("Names, sequences, and quality scores are not of the same length.")
  }
  # Append '@' to the start of sequence names.
  names<-paste0("@",names)
  # Create an empty data frame to store what will be written to a text file.
  df<-as.data.frame(matrix(NA,nrow=4*length(names),ncol=1))
  # Add names to the text file data frame.
  df$V1[seq(from=1,to=nrow(df),by=4)]<-names
  # Add sequences to the text file data frame.
  df$V1[seq(from=2,to=nrow(df),by=4)]<-as.character(sequences)
  # Add '+' separator between sequences and quality scores in the text file data frame.
  df$V1[seq(from=3,to=nrow(df),by=4)]<-"+"
  # Add quality scores to the text file data frame.
  df$V1[seq(from=4,to=nrow(df),by=4)]<-as.character(Qscores)
  # Write out fastq file.
  write.table(x=df,file=file,row.names=F,col.names=F,quote=F)
}

# Define function for displaying elasped time.
## Arguments:
### result_of_initial_proc.time: Result of an earlier call of the proc.time function.
## Returns: A character string of formatted elasped time.
elapsed.time<-function(result_of_initial_proc.time){
  # Get number of seconds since the initial proc.time result.
  tot_num_sec<-c(proc.time()-result_of_initial_proc.time)[3]
  # Get number of full hours since initial proc.time result.
  num_hrs<-floor(tot_num_sec/3600)
  # Get number of full minutes minus full hours since initial proc.time result.
  num_min<-floor((tot_num_sec-3600*num_hrs)/60)
  # Get number of full seconds minus full hours and full minutes since initial proc.time result.
  num_sec<-floor(tot_num_sec-3600*num_hrs-60*num_min)
  # Get formatted elapsed time.
  elapsed_time<-paste0(num_hrs," hrs, ",num_min," min, ",num_sec," sec")
  # Return character string of formatted elasped time.
  return(elapsed_time)
}

# Define function for assigning taxa to forward and reverse reads by exact matching with
# a reference dataset.
## Arguments:
### forward_reads: A vector of forward reads.
### reverse_reads: A vector of reverse reads.
### reference_df: A data frame containing the reference data. The following fields must be included
###   in the data frame: "Name", "ForwardRead", "ReverseRead".
## Returns: A vector of taxonomic assignments. Taxa names are taken from the reference dataset. If
##  no taxonomic assignment is made for a sequence, an NA is returned in the vector.
assign.taxa<-function(forward_reads,reverse_reads,reference_df){
  # Check that reference_df is a data frame.
  if(!is.data.frame(reference_df)){
    stop("reference_df is not a data frame. reference_df must be a data frame.")
  }
  # Check that the reference data frame contains all required fields.
  if(!all(c("Name","ForwardRead","ReverseRead") %in% colnames(reference_df))){
    stop('The reference data frame does not contain all required fields. The fields "Name", "ForwardRead", and "ReverseRead" are required.')
  }
  # Check that forward and reverse reads are of the same length.
  if(!length(forward_reads)==length(reverse_reads)){
    stop("Forward and reverse reads are not of the same length.")
  }
  # Create a data frame of forward and reverse reads.
  read.df<-data.frame(ForwardRead=forward_reads,ReverseRead=reverse_reads)
  # Ensure that the names field in the reference data frame is formatted as character.
  reference_df$Name<-as.character(reference_df$Name)
  # Define function to map forward and reverse reads to reference data.
  map.to.reference<-function(read.df_row,reference_df){
    # Get a vector of TRUE or FALSE's that states whether forward and reverse reads match those
    # of the reference data.
    matches<-read.df_row["ForwardRead"]==reference_df$ForwardRead &
      read.df_row["ReverseRead"]==reference_df$ReverseRead
    # If there is a match between forward and reverse reads and the reference data...
    if(sum(matches)==1){
      # Get the name of the matching taxon from the reference data.
      matching_name<-reference_df$Name[matches]
    }else{ # Otherwise...
      # Set the matching taxon name to NA.
      matching_name<-NA
    }
    # Return the matching taxon name.
    return(matching_name)
  }
  # Apply the reference mapping function to all provided forward and reverse reads.
  matching_names<-apply(X=read.df,MARGIN=1,FUN=map.to.reference,reference_df=reference_df)
  # Return the matching taxon names.
  return(matching_names)
}
