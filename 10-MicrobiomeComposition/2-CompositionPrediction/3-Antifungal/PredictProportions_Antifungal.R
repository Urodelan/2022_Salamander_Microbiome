##################################################################################
### Predict Proportions From 16S Dirichlet-Multinomial Regression Model Output ###
##################################################################################

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Output/Antifungal"

# Set barcoding region (Antifungal).
barcode_region<-"Antifungal"

# Declare directory to composition models.
composition_model_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition"

# Declare path to sample metadata in-file (with .csv extension).
metadata_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/SampleMetadata2.csv"

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

# Load metadata.
metadata<-read.csv(metadata_in_file)

# Create storage data frame for predicted proportions.
props<-data.frame(NULL)

# Loop through each sample type.
for(i in c("Salamander","Water","Substrate")){
  
  # Report preparing for predictions.
  print(paste("- Preparing for",tolower(i),"predictions."))
  
  # Load read data.
  print(paste("-- Loading",tolower(i),"read data."))
  reads<-read.csv(
    paste0(composition_model_directory,"/",i,"/Antifungal/Take3/",i,"_Antifungal_Table.csv"),
    row.names=1,check.names=F)
  
  # Load scaling data.
  print(paste("-- Loading",tolower(i),"scaling data."))
  stratum_scaling<-read.csv(
    paste0(composition_model_directory,"/",i,"/Antifungal/Take3/stratum_scaling.csv"),
    row.names=1)
  pred_scaling<-read.csv(
    paste0(composition_model_directory,"/",i,"/Antifungal/Take3/pred_scaling.csv"))
  
  # Subset and order metadata to match the samples included in the reads data.
  metadata_sub<-metadata[match(row.names(reads),metadata$SampleID),]
  
  # Check that metadata and read data samples match.
  if(!identical(as.character(metadata_sub$SampleID),row.names(reads))) stop("The samples in the read and metadata files are not the same.")
  
  # Report preparing predictors.
  print(paste("-- Preparing",tolower(i),"predictors."))
  
  # If the current sample type is salamander.
  if(i=="Salamander"){
    
    # Store salamander taxa names.
    sal_taxa_names<-colnames(reads)
    
    # Get unique combinations of predictors.
    meta<-unique(metadata_sub[,c("Date","Site","Age","LM")])
    row.names(meta)<-1:nrow(meta)
    
    # Format metadata fields.
    ## Format Date as a date.
    meta$Date<-as.Date(meta$Date,format="%Y-%m-%d")
    ## Ensure that Site is a factor.
    meta$Site<-as.factor(meta$Site)
    ## Format Age as a factor.
    meta$Age<-as.factor(meta$Age)
    ## Format LM as a factor.
    meta$LM<-as.factor(meta$LM)
    
    # Prepare non-stratum predictors.
    
    # Convert date to number of weeks since June 9th, 2018.
    meta$Week<-as.vector(meta$Date-as.Date("2018-06-09"))/7
    
    # Get site as a binary variable.
    meta$Site_binary<-as.numeric(meta$Site)-1
    ## Gibson: 0
    ## Ponds: 1
    
    # Get age as a continuous variable.
    meta$Age_numeric<-as.numeric(meta$Age)-1
    ## Age-0: 0
    ## Age-1: 1
    ## Age-2+: 2
    
    # Get LM as a binary variable.
    meta$LM_binary<-as.numeric(meta$LM)-1
    ## Larvae/Neotene: 0
    ## Metamorphosed: 1
    
    # Preparing stratum predictors.
    
    # Create stratum predictor matrix.
    stratum_matrix<-as.data.frame(matrix(data=0,nrow=nrow(meta),ncol=9))
    colnames(stratum_matrix)<-paste0("Stratum",1:9)
    
    # Copy stratum matrix for scaling.
    stratum_matrix_scaled<-stratum_matrix
    
    # Center stratum predictors (no scaling by SD).
    for(j in 1:ncol(stratum_matrix)){
      stratum_matrix_scaled[,j]<-stratum_matrix[,j]-stratum_scaling$Mean[j]
    }
    
    # Get just age, LM, site, and week.
    pred_matrix<-meta[,c("Age_numeric","LM_binary","Site_binary","Week")]
    
    # Get mean and SD of age.
    age_mean<-pred_scaling$Value[pred_scaling$Predictor=="Age" & pred_scaling$Parameter=="Mean"]
    age_sd<-pred_scaling$Value[pred_scaling$Predictor=="Age" & pred_scaling$Parameter=="SD"]
    
    # Get mean of binary LM.
    LM_mean<-pred_scaling$Value[pred_scaling$Predictor=="LM"]
    
    # Get mean of binary site.
    site_mean<-pred_scaling$Value[pred_scaling$Predictor=="Site"]
    
    # Get mean and SD of week.
    week_mean<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="Mean"]
    week_sd<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="SD"]
    
    # Scale the predictor matrix.
    pred_matrix_scaled<-data.frame(Age=(pred_matrix$Age_numeric-age_mean)/age_sd,
                                   LM=pred_matrix$LM_binary-LM_mean,
                                   Site=pred_matrix$Site_binary-site_mean,
                                   Week=(pred_matrix$Week-week_mean)/week_sd)
    
    # Add interaction terms to the scaled predictor matrix.
    pred_matrix_scaled<-as.data.frame(model.matrix(~Age*LM*Site*poly(Week,degree=2,raw=T)-1,
                                                   data=pred_matrix_scaled))
    
    # Call the polynomial for week 'Week1' and 'Week2' in the scaled predictor matrix.
    colnames(pred_matrix_scaled)<-gsub(pattern="poly\\(Week,(.*)\\)",
                                       replacement="Week",
                                       x=colnames(pred_matrix_scaled))
    
  } else { # If the current sample type is water or substrate.
    
    # Get unique combinations of predictors.
    meta<-unique(metadata_sub[,c("Date","Site")])
    row.names(meta)<-1:nrow(meta)
    
    # Remove predictors for Sept 29th.
    meta<-meta[meta$Date!="2018-09-29",]
    
    # If the samples are water, then add a record for a day at Gibson Lakes
    # for which water samples were discarded.
    if(i=="Water"){
      meta$Date<-as.Date(meta$Date)
      meta<-rbind(meta,data.frame(Date=as.Date("2018-07-14"),Site="Gibson Lakes"))
      meta$Date<-as.factor(meta$Date)
    }
    
    # Format metadata fields.
    ## Format Date as a date.
    meta$Date<-as.Date(meta$Date,format="%Y-%m-%d")
    ## Ensure that Site is a factor.
    meta$Site<-as.factor(meta$Site)
    
    # Prepare non-stratum predictors.
    
    # Convert date to number of weeks since June 9th, 2018.
    meta$Week<-as.vector(meta$Date-as.Date("2018-06-09"))/7
    
    # Get site as a binary variable.
    meta$Site_binary<-as.numeric(meta$Site)-1
    ## Gibson: 0
    ## Ponds: 1
    
    # Preparing stratum predictors.
    
    # Create stratum predictor matrix.
    stratum_matrix<-as.data.frame(matrix(data=0,nrow=nrow(meta),ncol=9))
    colnames(stratum_matrix)<-paste0("Stratum",1:9)
    
    # Copy stratum matrix for scaling.
    stratum_matrix_scaled<-stratum_matrix
    
    # Center stratum predictors (no scaling by SD).
    for(j in 1:ncol(stratum_matrix)){
      stratum_matrix_scaled[,j]<-stratum_matrix[,j]-stratum_scaling$Mean[j]
    }
    
    # Get just site and week.
    pred_matrix<-meta[,c("Site_binary","Week")]
    
    # Get mean of binary site.
    site_mean<-pred_scaling$Value[pred_scaling$Predictor=="Site"]
    
    # Get mean and SD of week.
    week_mean<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="Mean"]
    week_sd<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="SD"]
    
    # Scale the predictor matrix.
    pred_matrix_scaled<-data.frame(Site=pred_matrix$Site_binary-site_mean,
                                   Week=(pred_matrix$Week-week_mean)/week_sd)
    
    # Add interaction terms to the scaled predictor matrix.
    pred_matrix_scaled<-as.data.frame(model.matrix(~Site*poly(Week,degree=2,raw=T)-1,
                                                   data=pred_matrix_scaled))
    
    # Call the polynomial for week 'Week1' and 'Week2' in the scaled predictor matrix.
    colnames(pred_matrix_scaled)<-gsub(pattern="poly\\(Week,(.*)\\)",
                                       replacement="Week",
                                       x=colnames(pred_matrix_scaled))
    
  }
  
  # Load model selection information.
  print(paste("-- Loading",tolower(i),"model selection information."))
  model_selection<-read.csv(
    paste0(composition_model_directory,"/",i,"/Antifungal/Take3/model_selection.csv"),
    check.names=F)
  
  # Get the best fitting model from the model selection information.
  best_fit_model<-model_selection[model_selection$Model=="Best",]
  
  # Get the predictors of best-fitting model.
  best_fit_model_predictors<-
    best_fit_model[,3:(3+ncol(pred_matrix_scaled))]
  
  # If the best-fitting model includes any predictors.
  if(sum(best_fit_model_predictors) > 0){
    # Get the names of predictors of the best-fitting model.
    best_fit_model_predictors<-
      colnames(best_fit_model_predictors)[best_fit_model_predictors==1]
  } else {
    # If the best-fitting model has no predictors, then set the new predictor combination to NA.
    # This will result in an intercept-only model.
    best_fit_model_predictors<-NA
  }
  
  # Check if stratum is included in the predictor combination.
  stratum_included<-ifelse("Stratum" %in% best_fit_model_predictors,"Yes","No")
  
  # Check if any non-stratum predictors are included in the predictor combination.
  nonstratum_included<-ifelse(any(best_fit_model_predictors!="Stratum",na.rm=T),"Yes","No")
  
  # If non-stratum predictors are included, get new non-stratum predictor matrix
  # with the predictor combination.
  if(nonstratum_included=="Yes"){
    pred_matrix_scaled_reduced<-
      pred_matrix_scaled[,colnames(pred_matrix_scaled) %in% best_fit_model_predictors]
  }
  
  # Load model output.
  print(paste("-- Loading",tolower(i),"model output."))
  load(paste0(composition_model_directory,"/",i,"/Antifungal/Take3/mod_fit.RData"))
  
  # Get HMC output as a data frame.
  HMC<-as.data.frame(mod_fit)
  
  # Get intercept terms.
  intercept_betas<-HMC[,grepl("^beta_0",colnames(HMC))]
  
  # Get stratum betas if stratum is included as a predictor.
  if(stratum_included=="Yes") stratum_betas<-HMC[,grepl("^beta_stratum",colnames(HMC))]
  
  # Get non-stratum betas if included as predictors.
  if(nonstratum_included=="Yes") pred_betas<-HMC[,grepl("^beta_pred",colnames(HMC))]
  
  # Report beginning predictions.
  print(paste("-- Beginning",tolower(i),"predictions."))
  
  # Loop through each observation.
  for(n in 1:nrow(meta)){
    
    # Get stratum predictor values for the observation if stratum is included as a predictor.
    if(stratum_included=="Yes") stratum_values<-stratum_matrix_scaled[n,]
    
    # Get non-stratum predictor values for the observation if included as predictors.
    if(nonstratum_included=="Yes") pred_values<-pred_matrix_scaled_reduced[n,]
    
    # Create empty storage data frame for linear combination of HMC estimates.
    eta<-as.data.frame(matrix(data=NA,nrow=nrow(HMC),ncol=ncol(reads)))
    
    # Loop through each taxon.
    for(k in 1:ncol(reads)){
      
      # Get stratum betas for the taxon if stratum is included as a predictor.
      if(stratum_included=="Yes"){
        stratum_taxon_betas<-stratum_betas[,grepl(paste0("\\[",k,","),colnames(stratum_betas))]
      }
      
      # Get non-stratum betas for the taxon if included as predictors.
      if(nonstratum_included=="Yes"){
        pred_taxon_betas<-pred_betas[,grepl(paste0("\\[",k,","),colnames(pred_betas))]
      }
      
      # Calculate linear combination of HMC estimates.
      ## Initiate with eta as the intercept.
      eta[,k]<-intercept_betas[,k]
      ## If stratum is included as a predictor,
      ## add on the stratum effects.
      if(stratum_included=="Yes"){
        eta[,k]<-eta[,k]+as.matrix(stratum_taxon_betas) %*% t(stratum_values)
      }
      if(nonstratum_included=="Yes"){
        # If non-stratum predictors are included,
        # add on the effects.
        eta[,k]<-eta[,k]+as.matrix(pred_taxon_betas) %*% t(pred_values)
      }
      
    }
    
    # Get expected microbial proportions.
    pi<-as.data.frame(t(apply(X=eta,MARGIN=1,FUN=softmax)))
    
    # Provide taxon names to the expected microbial proportions.
    colnames(pi)<-colnames(reads)
    
    # If the current sample type is salamander.
    if(i=="Salamander"){
      
      # Get metadata for the prediction.
      prediction_metadata<-cbind(data.frame(Type=i,stringsAsFactors=F),meta[n,c("Date","Site","Age","LM")])
      
    } else { # If the current sample type is water or substrate.
      
      # Get metadata for the prediction.
      prediction_metadata<-cbind(data.frame(Type=i,stringsAsFactors=F),meta[n,c("Date","Site")])
      
      # Add a field for age and LM with NAs.
      prediction_metadata$Age<-NA
      prediction_metadata$LM<-NA
      
      # Get taxa which are missing from the non-salamander dataset.
      missing_taxa<-sal_taxa_names[!(sal_taxa_names %in% colnames(pi))]
      
      # If there are missing taxa.
      if(length(missing_taxa) > 0){
        # Create data frame containing missing taxa.
        missing_taxa_df<-as.data.frame(matrix(data=0,nrow=nrow(pi),ncol=length(missing_taxa)))
        colnames(missing_taxa_df)<-missing_taxa
        # Add the missing taxa to the data and order the same as the salamander data.
        pi<-cbind(pi,missing_taxa_df)
        pi<-pi[,sal_taxa_names]
      }
      
    }
    
    # Repeat prediction metadata as many times as there are HMC samples.
    prediction_metadata<-prediction_metadata[rep(1,nrow(pi)),]
    
    # Add field in prediction metadata for HMC sample.
    prediction_metadata$HMC_sample<-1:nrow(prediction_metadata)
    
    # Add prediction metadata to the predicted proportions.
    pi<-cbind(prediction_metadata,pi)
    
    # Rename pi rows.
    row.names(pi)<-1:nrow(pi)
    
    # Add predicted proportions to the storage data frame.
    props<-rbind(props,pi)
    
    # Report prediction progress.
    print(paste0("--- ",i," prediction progress: ",n," of ",nrow(meta),"."))
    
  }
  
}

# Order the rows of the predicted proportions data frame.
props<-props[order(props$Type,props$Date,props$Age,props$LM,props$HMC_sample),]

# Write out predicted proportions as a csv file.
print("- Writing out predicted proportions.")
write.csv(props,file=paste0("props_",barcode_region,".csv"),row.names=F)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
