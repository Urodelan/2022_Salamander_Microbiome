##############################################################################################
### Process Proportion Predictions From Antifungal Dirichlet-Multinomial Regression Output ###
##############################################################################################

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Output/Antifungal"

# Set barcoding region (Antifungal).
barcode_region<-"Antifungal"

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

# Ensure barcode region variable is set to Antifungal.
if(!(barcode_region=="Antifungal")) stop("Barcode region must be Antifungal.")

# Load proportions.
props<-read.csv(paste0("props_",barcode_region,".csv"),check.names=F)

### Summarize proportions with credible intervals.

# Create empty storage data frames.
## For proportional abundance summaries.
props_summary<-data.frame(NULL)

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
    
    # Report processing progress.
    print(paste0("-- ",i," processing progress: ",j," of ",nrow(unique_combinations),"."))
    
  }
  
}

# Write out predicted proportions as a csv file.
print("- Writing out processed proportions.")
write.csv(props_summary,file=paste0("props_summary_",barcode_region,".csv"),row.names=F)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
