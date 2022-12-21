#############################################################################################
### Process Proportion Predictions From 16S Dirichlet-Multinomial Regression Model Output ###
#############################################################################################

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Output/Region16S"

# Set barcoding region (16S or ITS).
barcode_region<-"16S"

# Path to microbe functions script.
MicrobeFunctionsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/HelperFunctions/MicrobeFunctions.R"

####################
### Begin Script ###
####################

# Report starting script.
print("Starting script...")

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# Set working directory.
setwd(working_directory)

# Load microbe functions.
source(MicrobeFunctionsPath)

# Define Hill's diversity index (with alpha=2) function.
hills2<-function(x) 1/sum(x^2)

# Ensure barcode region variable is set to either 16S or ITS.
if(!(barcode_region %in% c("16S","ITS"))) stop("Barcode region must be 16S or ITS.")

# Load proportions.
props<-read.csv(paste0("props_",barcode_region,".csv"),check.names=F)

# If the barcode region is 16S.
if(barcode_region=="16S"){
  
  # Read in antifungal status metadata.
  antifungal<-read.csv("/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/BdStatusMetadata.csv")
  
  # Assign antifungal status based on very high or low confidence for Bd inhibition.
  antifungal$Status<-ifelse(antifungal$ProbabilityOfBdInhibition <= 0.1,"NonInhibitory",ifelse(antifungal$ProbabilityOfBdInhibition >= 0.9,"Inhibitory","Uncertain"))
  
  # Get options for antifungal status, plus a category of 'other' taxa.
  status_options<-c(unique(antifungal$Status),"Other")
  
  # Create empty storage data frame for antifungal proportional abundance summaries.
  antifungal_summary<-data.frame(NULL)
  
}

### Summarize proportions with credible intervals.

# Create empty storage data frames.
## For proportional abundance summaries.
props_summary<-data.frame(NULL)
## For Hills diversity index summaries (for salamanders only).
hills2_summary<-data.frame(NULL)

# Loop through each sample type.
for(i in c("Salamander","Water","Substrate")){
  
  # Report beginning to process predictions.
  print(paste("- Beginning to process",tolower(i),"predictions."))
  
  # Subset to the sample type.
  sub<-props[props$Type==i,]
  
  # Get unique combinations of predictor levels.
  unique_combinations<-unique(sub[,c("Type","Date","Site","Age","LM")])
  row.names(unique_combinations)<-1:nrow(unique_combinations)
  
  # Loop through each unique predictor combination.
  for(j in 1:nrow(unique_combinations)){
    
    # Get unique predictor combination.
    unique_combination<-unique_combinations[j,]
    
    # If the sample type is salamander.
    if(i=="Salamander"){
      
      # Subset HMC samples to unique predictor combination.
      sub2<-sub[sub$Date==unique_combination$Date &
                  sub$Age==unique_combination$Age &
                  sub$LM==unique_combination$LM,]
      
    } else { # If the sample type is water or substrate.
      
      # Subset HMC samples to unique predictor combination.
      sub2<-sub[sub$Date==unique_combination$Date,]
      
    }
    
    # Get just HMC samples in the data frame.
    sub2<-sub2[,-which(colnames(sub2) %in% c("Type","Date","Site","Age","LM","HMC_sample"))]
    
    # Get quantiles of microbe proportional abundances.
    sub3<-as.data.frame(t(apply(X=sub2,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
    
    # Rename columns of the proportion summary data frame.
    colnames(sub3)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=colnames(sub3)))/100)
    
    # Repeat predictor information as many times as there are rows in the proportion summary data frame.
    unique_combination2<-unique_combination[rep(1,nrow(sub3)),]
    
    # Add a field to the predictor information data frame containing taxa names.
    unique_combination2$Taxon<-row.names(sub3)
    
    # Combine predictor information and proportional abundance summary data frames.
    sub3<-cbind(unique_combination2,sub3)
    
    # Rename the rows of the proportional abundance summary data frame.
    row.names(sub3)<-1:nrow(sub3)
    
    # Add proportional abundance summary to the storage data frame.
    props_summary<-rbind(props_summary,sub3)
    
    # If the barcode region is 16S.
    if(barcode_region=="16S"){
      
      # Create an intermediate empty storage list for calculate antifungal status proportions.
      sub5<-vector(mode="list",length=length(status_options))
      
      # Loop through each antifungal status.
      for(k in 1:length(sub5)){
        
        # If the antifungal status category is other.
        if(k==4){
          
          # Set the taxa of interest to be the 'other' taxa category.
          taxa_of_interest<-"Other"
          
        } else { # If the antifungal status category is not for the 'other' taxa.
          
          # Subset the antifungal metadata data frame to the current antifungal status.
          taxa_of_interest<-antifungal[antifungal$Status==status_options[k],]
          
          # Get just the taxa of interest.
          taxa_of_interest<-taxa_of_interest$Taxon
          
        }
        
        # Subset the proportions to just the taxa of interest.
        sub4<-as.data.frame(sub2[,colnames(sub2) %in% taxa_of_interest])
        
        # If there are taxa of interest in the proportions data frame.
        if(ncol(sub4) > 0){
          
          # Sum the proportions of the taxa of interest.
          sub4<-data.frame(Props=rowSums(sub4))
          
        } else { # If there are no taxa of interest in the proportions data frame.
          
          # Set the summed proportion of the taxa of interest to zero.
          sub4<-as.data.frame(matrix(data=0,nrow=nrow(sub2),ncol=1))
          
        }
        
        # Rename the data frame to the current antifungal status.
        colnames(sub4)<-status_options[k]
        
        # Store the current data frame in the storage list.
        sub5[[k]]<-sub4
        
      }
      
      # Create a data frame containing the proportions of taxa by antifungal status.
      sub5<-do.call(cbind,sub5)
      
      # Get quantiles of antifungal proportional abundances.
      sub6<-as.data.frame(t(apply(X=sub5,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
      
      # Rename columns of the proportion summary data frame.
      colnames(sub6)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=colnames(sub6)))/100)
      
      # Repeat predictor information as many times as there are rows in the antifungal
      # proportion summary data frame.
      unique_combination3<-unique_combination[rep(1,nrow(sub6)),]
      
      # Add a field to the predictor information data frame containing antifungal status.
      unique_combination3$Taxon<-row.names(sub6)
      
      # Combine predictor information and antifungal proportional abundance summary data frames.
      sub6<-cbind(unique_combination3,sub6)
      
      # Rename the rows of the antifungal proportional abundance summary data frame.
      row.names(sub6)<-1:nrow(sub6)
      
      # Add antifungal proportional abundance summary to the storage data frame.
      antifungal_summary<-rbind(antifungal_summary,sub6)
      
    }
    
    # If the sample type is salamander.
    if(i=="Salamander"){
      
      # Calculate Hills diversity index for each HMC sample.
      hills2_vector<-apply(X=sub2,MARGIN=1,FUN=hills2)
      
      # Summarize Hills diversity index across HMC samples.
      hills2_df<-as.data.frame(t(quantile(hills2_vector,probs=c(0.025,0.25,0.5,0.75,0.975))))
      
      # Rename columns of the Hills diversity index summary data frame.
      colnames(hills2_df)<-
        paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=colnames(hills2_df)))/100)
      
      # Add predictor information to Hills diversity index summary data frame.
      hills2_df<-cbind(unique_combination,hills2_df)
      
      # Remove the sample type field from the Hills diversity index summary data frame,
      # since Hills diversity index is only calculated for salamanders.
      hills2_df<-hills2_df[,-which(colnames(hills2_df)=="Type")]
      
      # Add Hills diversity index summary to the storage data frame.
      hills2_summary<-rbind(hills2_summary,hills2_df)
      
    }
    
    # Report processing progress.
    print(paste0("-- ",i," processing progress: ",j," of ",nrow(unique_combinations),"."))
    
  }
  
}

# Write out predicted proportions as a csv file.
print("- Writing out processed proportions.")
write.csv(props_summary,file=paste0("props_summary_",barcode_region,".csv"),row.names=F)

# If the barcode region is 16S, write out predicted antifungal proportions as a csv file.
if(barcode_region=="16S"){
  print("- Writing out antifungal proportions.")
  write.csv(antifungal_summary,file=paste0("antifungal_summary_",barcode_region,".csv"),row.names=F)
}

# Write out Hills diversity as a csv file.
print("- Writing out Hills diversity.")
write.csv(hills2_summary,file=paste0("hills2_summary_",barcode_region,".csv"),row.names=F)

# Report beginning calculate Bray-Curtis dissimilarity for salamanders.
print("- Calculating Bray-Curtis dissimilarity for salamanders.")

# Subset to salamaders.
props<-props[props$Type=="Salamander",]

# Get unique age and LM combinations.
unique_combinations<-unique(props[,c("Age","LM")])

# Sort unique age and LM combinations by age and LM.
unique_combinations<-unique_combinations[order(unique_combinations$Age,unique_combinations$LM),]
row.names(unique_combinations)<-1:nrow(unique_combinations)

# Load the vegan package.
library(vegan)

# Create a storage list for Bray-Curtis dissimilarity summaries.
bc_summary<-vector(mode="list",length=nrow(unique_combinations))

# Loop through each unique age and LM combination.
for(i in 1:nrow(unique_combinations)){
  
  # Get unique age and LM combination.
  unique_combination<-unique_combinations[i,]
  
  # Subset HMC samples to unique age and LM combination.
  sub<-props[props$Age==unique_combination$Age & props$LM==unique_combination$LM,]
  
  # Format date field.
  sub$Date<-as.Date(sub$Date)
  
  # Only proceed if there are at least two dates to compare.
  if(length(unique(sub$Date))!=1){
    
    # Order HMC samples by sample and date.
    sub<-sub[order(sub$HMC_sample,sub$Date),]
    
    # Create the first intermediate storage list for calculating Bray-Curtis dissimilarity.
    bc_int1<-vector(mode="list",length=length(unique(sub$HMC_sample)))
    
    # Loop through each HMC sample.
    for(j in 1:length(unique(sub$HMC_sample))){
      
      # Get predictions for sorted dates for each HMC sample
      # and store in the first intermediate list.
      bc_int1[[j]]<-subset(sub,HMC_sample==j)[7:ncol(sub)]
      row.names(bc_int1[[j]])<-1:nrow(bc_int1[[j]])
      
    }
    
    # Compute Bray-Curtis dissimilarity indices.
    bc_int2<-lapply(bc_int1,vegdist,method="bray")
    
    # Get the 2.5th percentile of Bray-Curtis dissimilarity indices.
    bc_int3_0.025<-bc_int2[[1]]
    for(j in 1:length(bc_int3_0.025)){
      bc_int3_0.025[j]<-quantile(unlist(lapply(bc_int2,"[",j)),probs=0.025)
    }
    
    # Get the 25th percentile of Bray-Curtis dissimilarity indices.
    bc_int3_0.25<-bc_int2[[1]]
    for(j in 1:length(bc_int3_0.25)){
      bc_int3_0.25[j]<-quantile(unlist(lapply(bc_int2,"[",j)),probs=0.25)
    }
    
    # Get the 50th percentile of Bray-Curtis dissimilarity indices.
    bc_int3_0.5<-bc_int2[[1]]
    for(j in 1:length(bc_int3_0.5)){
      bc_int3_0.5[j]<-quantile(unlist(lapply(bc_int2,"[",j)),probs=0.5)
    }
    
    # Get the 75th percentile of Bray-Curtis dissimilarity indices.
    bc_int3_0.75<-bc_int2[[1]]
    for(j in 1:length(bc_int3_0.75)){
      bc_int3_0.75[j]<-quantile(unlist(lapply(bc_int2,"[",j)),probs=0.75)
    }
    
    # Get the 97.5th percentile of Bray-Curtis dissimilarity indices.
    bc_int3_0.975<-bc_int2[[1]]
    for(j in 1:length(bc_int3_0.975)){
      bc_int3_0.975[j]<-quantile(unlist(lapply(bc_int2,"[",j)),probs=0.975)
    }
    
    # Create a list of Bray-Curtis dissimilarity index quantiles.
    bc<-list("Quantile_0.025"=bc_int3_0.025,
             "Quantile_0.25"=bc_int3_0.25,
             "Quantile_0.5"=bc_int3_0.5,
             "Quantile_0.75"=bc_int3_0.75,
             "Quantile_0.975"=bc_int3_0.975)
    
    # Loop through each age and LM combination Bray-Curtis summaries.
    for(j in 1:length(bc)){
      
      # Get element of Bray-Curtis summary storage list.
      bc_sub<-as.matrix(bc[[j]])
      
      # Loop through each column of the list's element.
      for(k in 1:ncol(bc_sub)){
        
        # Remove repeated values.
        bc_sub[1:which(bc_sub[,k]==0),k]<-NA
        
      }
      
      # Get dates from Gibson Lakes.
      sub.GibsonDates<-unique(subset(sub,Site=="Gibson Lakes")$Date)
      
      # Get dates from both lakes.
      unique.dates<-unique(sub$Date)
      
      # Create a vector denoting which dates are from which lake.
      Gibson.v.Ponds_Dates<-ifelse(unique.dates %in% sub.GibsonDates,"Gibson","Ponds")
      
      # Combine lake names and their associated dates for new row and column names.
      new.names<-paste0(Gibson.v.Ponds_Dates,"_",unique.dates)
      
      # Rename the Bray-Curtis matrix columns and rows.
      colnames(bc_sub)<-new.names
      row.names(bc_sub)<-new.names
      
      # Give element back to the Bray-Curtis summary storage list.
      bc[[j]]<-bc_sub
      
    }
    
    # Get age and LM status of the Bray-Curtis list.
    element_age<-unique_combination$Age
    element_LM<-unique_combination$LM
    
    # Decide on larvae or neotene based on salamander age.
    if(element_LM=="Larvae/Neotene"){
      element_LM<-ifelse(element_age=="Age-2+","Neotene","Larvae")
    }
    
    # Name the element of the Bray-Curtis list with salamander age and LM.
    names(bc_summary)[i]<-paste(element_age,element_LM)
    
    # Add element to the Bray-Curtis list.
    bc_summary[[i]]<-bc
    
  }
  
}

# If there are any NA names in the Bray-Curtis list.
if(any(is.na(names(bc_summary)))){
  
  # Identify which list elements have NA names.
  elements.to.remove<-which(is.na(names(bc_summary)))
  
  # Remove elements with NA names from the Bray-Curtis list.
  bc_summary<-bc_summary[-elements.to.remove]
  
}

# Save Bray-Curtis dissimilarity summaries as an RData file.
print("- Writing out Bray-Curtis dissimilarity.")
save(bc_summary,file=paste0("bc_summary_",barcode_region,".RData"))

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
