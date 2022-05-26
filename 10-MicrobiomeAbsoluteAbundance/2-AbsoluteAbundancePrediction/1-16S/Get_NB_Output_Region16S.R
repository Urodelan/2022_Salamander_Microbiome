########################################################################################
### Work Up Output from Negative Binomial Regression Model with Salamander 16S Reads ###
########################################################################################

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Output/Absolute/Region16S"

# Declare path to reads in-file (with .csv extension).
reads_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Absolute/Salamander/Region16S/Take3/Salamander_ASVTable_16S_Rate_2.csv"

# Declare path to directory with model outputs.
path_to_model_output<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Absolute/Salamander/Region16S/Take3/Parallel"

# Declare barcode region ("16S" or "ITS").
barcode_region<-"16S"

# Declare path to sample metadata in-file (with .csv extension).
metadata_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/SampleMetadata2.csv"

# Declare path to antifungal status metadata in-file (with .csv extension).
path_to_antifungal_status_metadata<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/BdStatusMetadata.csv"

# Path to microbe functions script.
MicrobeFunctionsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/HelperFunctions/MicrobeFunctions.R"

####################
### Begin Script ###
####################

# Report starting script.
print("Starting script...")

# Set working directory.
setwd(working_directory)

# Load microbe functions.
source(MicrobeFunctionsPath)

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# Check that barcode region is valid.
if(!(barcode_region %in% c("16S","ITS"))) stop("User-defined barcode region is not valid.")

# Load read data.
print(paste0("- Loading read data: ",elapsed.time(ptm)))
reads<-read.csv(reads_in_file,row.names=1,check.names=F)

# Create a data frame with a field for the taxon name and its numeric values.
taxa<-data.frame(Taxon_number=1:101,Taxon=colnames(reads)[1:101],stringsAsFactors=F)

# Locate convergence files.
print(paste0("- Locating convergence files: ",elapsed.time(ptm)))
convergence_files<-list.files(path=paste0(path_to_model_output,"/Convergence"),pattern=".csv$",full.names=T)

# Load convergence information.
print(paste0("- Loading convergence information: ",elapsed.time(ptm)))

# Create an empty storage data frame.
convergence<-data.frame(NULL)

# Loop through each convergence file.
for(i in convergence_files){
  # Read in the convergence file.
  df<-read.csv(file=i,row.names=1,check.names=F)
  # Turn the row names into a parameters field.
  df$Parameter<-row.names(df)
  # Rename the rows.
  row.names(df)<-1:nrow(df)
  # Create a field for the taxon number
  df$Taxon_number<-as.numeric(sub(pattern="^Taxon",replacement="",x=strsplit(x=basename(i),split="_")[[1]][2]))
  # Add the read-in convergence file to the storage data frame.
  convergence<-rbind(convergence,df)
}

# Format convergence information.
print(paste0("- Formatting convergence information: ",elapsed.time(ptm)))

# Add taxa names to the convergence information.
convergence$Taxon<-taxa$Taxon[match(convergence$Taxon_number,taxa$Taxon_number)]

# Sort convergence information by taxon number.
convergence<-convergence[order(convergence$Taxon_number),]

# Re-order fields of convergence information.
convergence<-convergence[,c("Taxon","Parameter","Upper C.I.")]

# Rename last field to friendlier characters.
colnames(convergence)[3]<-"Upper_CI"

# Reshape convergence information to wide format.
convergence<-reshape(data=convergence,idvar="Taxon",timevar="Parameter",direction="wide")

# Turn taxon into the convergence row names.
row.names(convergence)<-convergence$Taxon

# Remove the taxon field from the convergence information.
convergence<-convergence[,-which(colnames(convergence)=="Taxon")]

# Remove the upper CI prefixes from the convergence field names.
colnames(convergence)<-sub(pattern="^Upper_CI.",replacement="",x=colnames(convergence))

# Get convergence information for just the regression coefficients.
convergence_betas<-convergence[,grepl(pattern="^beta_",x=colnames(convergence))]

# Get convergence information for just the inclusion variables.
convergence_IVs<-convergence[,grepl(pattern="^IV_",x=colnames(convergence))]

# Summarize convergence information for the regression coefficients.
## Get number of convergence diagnostics greater than 1.05.
beta_conv_gt_1.05<-as.data.frame(apply(X=convergence_betas,MARGIN=1,FUN=function(x) sum(x > 1.05,na.rm=T)))
colnames(beta_conv_gt_1.05)<-"Num_conv_gt_1.05"
## Get number of convergence diagnostics greater than 1.1.
beta_conv_gt_1.1<-as.data.frame(apply(X=convergence_betas,MARGIN=1,FUN=function(x) sum(x > 1.1,na.rm=T)))
colnames(beta_conv_gt_1.1)<-"Num_conv_gt_1.1"
## Get number of convergence diagnostics greater than 1.2.
beta_conv_gt_1.2<-as.data.frame(apply(X=convergence_betas,MARGIN=1,FUN=function(x) sum(x > 1.2,na.rm=T)))
colnames(beta_conv_gt_1.2)<-"Num_conv_gt_1.2"
## Get number of convergence diagnostics greater than 1.5.
beta_conv_gt_1.5<-as.data.frame(apply(X=convergence_betas,MARGIN=1,FUN=function(x) sum(x > 1.5,na.rm=T)))
colnames(beta_conv_gt_1.5)<-"Num_conv_gt_1.5"
## Get number of convergence diagnostics greater than 1.2.
beta_conv_gt_2<-as.data.frame(apply(X=convergence_betas,MARGIN=1,FUN=function(x) sum(x > 2,na.rm=T)))
colnames(beta_conv_gt_2)<-"Num_conv_gt_2"
## Combine this information together.
beta_conv_gt<-cbind(beta_conv_gt_1.05,beta_conv_gt_1.1,beta_conv_gt_1.2,beta_conv_gt_1.5,beta_conv_gt_2)
## Add a field for the total number of parameters being considered.
beta_conv_gt$Num_params<-ncol(convergence_betas)

# Summarize convergence information for the inclusion variables.
## Get number of convergence diagnostics greater than 1.05.
IV_conv_gt_1.05<-as.data.frame(apply(X=convergence_IVs,MARGIN=1,FUN=function(x) sum(x > 1.05,na.rm=T)))
colnames(IV_conv_gt_1.05)<-"Num_conv_gt_1.05"
## Get number of convergence diagnostics greater than 1.1.
IV_conv_gt_1.1<-as.data.frame(apply(X=convergence_IVs,MARGIN=1,FUN=function(x) sum(x > 1.1,na.rm=T)))
colnames(IV_conv_gt_1.1)<-"Num_conv_gt_1.1"
## Get number of convergence diagnostics greater than 1.2.
IV_conv_gt_1.2<-as.data.frame(apply(X=convergence_IVs,MARGIN=1,FUN=function(x) sum(x > 1.2,na.rm=T)))
colnames(IV_conv_gt_1.2)<-"Num_conv_gt_1.2"
## Get number of convergence diagnostics greater than 1.5.
IV_conv_gt_1.5<-as.data.frame(apply(X=convergence_IVs,MARGIN=1,FUN=function(x) sum(x > 1.5,na.rm=T)))
colnames(IV_conv_gt_1.5)<-"Num_conv_gt_1.5"
## Get number of convergence diagnostics greater than 1.2.
IV_conv_gt_2<-as.data.frame(apply(X=convergence_IVs,MARGIN=1,FUN=function(x) sum(x > 2,na.rm=T)))
colnames(IV_conv_gt_2)<-"Num_conv_gt_2"
## Combine this information together.
IV_conv_gt<-cbind(IV_conv_gt_1.05,IV_conv_gt_1.1,IV_conv_gt_1.2,IV_conv_gt_1.5,IV_conv_gt_2)
## Add a field for the total number of parameters being considered.
IV_conv_gt$Num_params<-ncol(convergence_IVs)

# Write out convergence information.
print(paste0("- Writing out convergence information: ",elapsed.time(ptm)))
## For all parameters.
write.csv(x=convergence,file=paste0("Region",barcode_region,"_NB_Reg_Convergence.csv"),row.names=T)
## For regression coefficients.
write.csv(x=beta_conv_gt,file=paste0("Region",barcode_region,"_NB_Reg_Beta_Convergence.csv"),row.names=T)
## For inclusion variables.
write.csv(x=IV_conv_gt,file=paste0("Region",barcode_region,"_NB_Reg_IV_Convergence.csv"),row.names=T)

# Locate PIP files.
print(paste0("- Locating PIP files: ",elapsed.time(ptm)))
PIP_files<-list.files(path=paste0(path_to_model_output,"/PIP"),pattern=".csv$",full.names=T)

# Load PIP information.
print(paste0("- Loading PIP information: ",elapsed.time(ptm)))

# Create an empty storage data frame.
PIP<-data.frame(NULL)

# Loop through each PIP file.
for(i in PIP_files){
  # Read in the PIP file.
  df<-read.csv(file=i,row.names=1,check.names=F)
  # Turn the row names into a predictors field.
  df$Predictor<-row.names(df)
  # Rename the rows.
  row.names(df)<-1:nrow(df)
  # Create a field for the taxon number
  df$Taxon_number<-as.numeric(sub(pattern="^Taxon",replacement="",x=strsplit(x=basename(i),split="_")[[1]][2]))
  # Add the read-in PIP file to the storage data frame.
  PIP<-rbind(PIP,df)
}

# Format PIP information.
print(paste0("- Formatting PIP information: ",elapsed.time(ptm)))

# Add taxa names to the PIP information.
PIP$Taxon<-taxa$Taxon[match(PIP$Taxon_number,taxa$Taxon_number)]

# Sort PIP information by taxon number.
PIP<-PIP[order(PIP$Taxon_number),]

# Re-order fields of PIP information.
PIP<-PIP[,c("Taxon","Predictor","PIP")]

# Reshape PIP information to wide format.
PIP<-reshape(data=PIP,idvar="Taxon",timevar="Predictor",direction="wide")

# Turn taxon into the PIP row names.
row.names(PIP)<-PIP$Taxon

# Remove the taxon field from the PIP information.
PIP<-PIP[,-which(colnames(PIP)=="Taxon")]

# Remove the upper PIP prefixes from the PIP field names.
colnames(PIP)<-sub(pattern="^PIP.",replacement="",x=colnames(PIP))

# Get means of PIPs.
PIP_mean<-as.data.frame(apply(X=PIP,MARGIN=2,FUN=mean))
colnames(PIP_mean)<-"PIP.Mean"

# Get SDs of PIPs.
PIP_SD<-as.data.frame(apply(X=PIP,MARGIN=2,FUN=sd))
colnames(PIP_SD)<-"PIP.SD"

# Combine PIP mean and SD information.
PIP_summary<-cbind(PIP_mean,PIP_SD)

# Write out PIP information.
print(paste0("- Writing out PIP information: ",elapsed.time(ptm)))
## All PIP information.
write.csv(x=PIP,file=paste0("Region",barcode_region,"_NB_Reg_PIP.csv"),row.names=T)
## PIP summaries.
write.csv(x=PIP_summary,file=paste0("Region",barcode_region,"_NB_Reg_PIP_Summary.csv"),row.names=T)

# Locate MCMC quantiles files.
print(paste0("- Locating MCMC quantiles files: ",elapsed.time(ptm)))
mcmc_quantiles_files<-list.files(path=paste0(path_to_model_output,"/MCMC_Quantiles"),pattern=".csv$",full.names=T)

# Load density predictions.
print(paste0("- Loading MCMC quantiles files: ",elapsed.time(ptm)))

# Create an empty storage data frame.
mcmc_quantiles<-data.frame(NULL)

# Loop through each MCMC quantiles file.
for(i in mcmc_quantiles_files){
  # Read in the MCMC quantiles file.
  df<-read.csv(file=i,check.names=F,row.names=1)
  # Turn the row names into a parameters field.
  df$Parameter<-row.names(df)
  # Rename the rows.
  row.names(df)<-1:nrow(df)
  # Create a field for the taxon number.
  df$Taxon_number<-as.numeric(sub(pattern="^Taxon",replacement="",x=strsplit(x=basename(i),split="_")[[1]][2]))
  # Add the read-in MCMC quantiles file to the storage data frame.
  mcmc_quantiles<-rbind(mcmc_quantiles,df)
}

# Add taxa names to the MCMC quantiles.
mcmc_quantiles$Taxon<-taxa$Taxon[match(mcmc_quantiles$Taxon_number,taxa$Taxon_number)]

# Sort MCMC quantiles by taxon number.
mcmc_quantiles<-mcmc_quantiles[order(mcmc_quantiles$Taxon_number),]

# Re-order fields of MCMC quantiles.
mcmc_quantiles<-mcmc_quantiles[,c("Taxon","Parameter","2.5%","25%","50%","75%","97.5%")]

# Write out MCMC quantiles.
print(paste0("- Writing out MCMC quantiles: ",elapsed.time(ptm)))
write.csv(x=mcmc_quantiles,file=paste0("Region",barcode_region,"_NB_Reg_MCMC_Quantiles.csv"),row.names=F)

# Locate density prediction files.
print(paste0("- Locating density prediction files: ",elapsed.time(ptm)))
density_files<-list.files(path=paste0(path_to_model_output,"/Density_Predictions"),pattern=".csv$",full.names=T)

# Load density predictions.
print(paste0("- Loading density predictions: ",elapsed.time(ptm)))

# Create an empty storage data frame.
density<-data.frame(NULL)

# Loop through each density prediction file.
for(i in density_files){
  # Read in the density prediction file.
  df<-read.csv(file=i,check.names=F,stringsAsFactors=F)
  # Create a field for the taxon number
  df$Taxon_number<-as.numeric(sub(pattern="^Taxon",replacement="",x=strsplit(x=basename(i),split="_")[[1]][2]))
  # Add the read-in density prediction file to the storage data frame.
  density<-rbind(density,df)
}

# Sort records.
## Format the date field.
density$Date<-as.Date(density$Date)
## Sort density predictions by date, age, life stage, and taxon number.
density<-density[order(density$Date,density$Age,density$LM,density$Taxon_number),]

# Remove the taxon number field from the density predictions.
density<-density[,-which(colnames(density)=="Taxon_number")]

# Write out density predictions.
print(paste0("- Writing out density predictions: ",elapsed.time(ptm)))
write.csv(x=density,file=paste0("Region",barcode_region,"_NB_Reg_Density_Predictions.csv"),row.names=F)

# If the barcode region is 16S.
if(barcode_region=="16S"){
  
  # Report preparing metadata for prediction.
  print(paste0("- Preparing metadata for antifungal prediction: ",elapsed.time(ptm)))
  
  # Load metadata.
  metadata<-read.csv(metadata_in_file)
  
  # Format metadata fields.
  ## Format Date as a date.
  metadata$Date<-as.Date(metadata$Date,format="%Y-%m-%d")
  ## Ensure that Site is a factor.
  metadata$Site<-as.factor(metadata$Site)
  ## Format Stratum as a factor.
  metadata$Stratum<-as.factor(metadata$Stratum)
  ## Format Age as a factor.
  metadata$Age<-as.factor(metadata$Age)
  ## Format LM as a factor.
  metadata$LM<-as.factor(metadata$LM)
  
  # Subset metadata to just samples included in the reads data.
  metadata<-metadata[metadata$SampleID %in% row.names(reads),]
  
  # Order metadata to match read data.
  metadata<-metadata[match(metadata$SampleID,row.names(reads)),]
  
  # Check that metadata and read data samples match.
  if(!identical(as.character(metadata$SampleID),row.names(reads))) stop("The samples in the read and metadata files are not the same.")
  
  # Read in stratum scaling parameters.
  stratum_scaling<-read.csv(file=paste0(path_to_model_output,"/Stratum_Scaling/Region16S_Taxon1_StratumScaling.csv"),row.names=1)
  
  # Read in non-stratum scaling parameters.
  pred_scaling<-read.csv(file=paste0(path_to_model_output,"/Predictor_Scaling/Region16S_Taxon1_PredictorScaling.csv"))
  
  # Get mean and SD of age.
  age_mean<-pred_scaling$Value[pred_scaling$Predictor=="Age" & pred_scaling$Parameter=="Mean"]
  age_sd<-pred_scaling$Value[pred_scaling$Predictor=="Age" & pred_scaling$Parameter=="SD"]
  
  # Get mean of binary LM.
  LM_mean<-pred_scaling$Value[pred_scaling$Predictor=="LM" & pred_scaling$Parameter=="Mean"]
  
  # Get mean of binary site.
  site_mean<-pred_scaling$Value[pred_scaling$Predictor=="Site" & pred_scaling$Parameter=="Mean"]
  
  # Get mean and SD of week.
  week_mean<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="Mean"]
  week_sd<-pred_scaling$Value[pred_scaling$Predictor=="Week" & pred_scaling$Parameter=="SD"]
  
  # Get unique combinations of predictors.
  meta<-unique(metadata[,c("Date","Site","Age","LM")])
  row.names(meta)<-1:nrow(meta)
  
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
  
  # Scale the predictor matrix.
  # Continuous covariates are scaled by two standard deviations.
  pred_matrix_scaled<-data.frame(Age=(pred_matrix$Age_numeric-age_mean)/(2*age_sd),
                                 LM=pred_matrix$LM_binary-LM_mean,
                                 Site=pred_matrix$Site_binary-site_mean,
                                 Week=(pred_matrix$Week-week_mean)/(2*week_sd))
  
  # Add interaction terms to the scaled predictor matrix.
  pred_matrix_scaled<-as.data.frame(model.matrix(~poly(Age,degree=2,raw=T)*LM*Site*poly(Week,degree=4,raw=T)-1,data=pred_matrix_scaled))
  
  # Split predictor matrix field names by colon.
  col.names.split<-strsplit(x=colnames(pred_matrix_scaled),split=":")
  
  # Loop through each predictor matrix field.
  for(i in 1:ncol(pred_matrix_scaled)){
    # Get the split column name.
    col.name<-col.names.split[i][[1]]
    # Simplify the polynomial name for age.
    col.name<-gsub(pattern="poly\\(Age,(.*)\\)",replacement="Age",x=col.name)
    # Simplify the polynomial name for week.
    col.name<-gsub(pattern="poly\\(Week,(.*)\\)",replacement="Week",x=col.name)
    # Collapse the split column name by colon.
    col.name<-paste0(col.name,collapse=":")
    # Assign the simplified column name to the predictor matrix field.
    colnames(pred_matrix_scaled)[i]<-col.name
  }
  
  # Read in antifungal status metadata.
  antifungal<-read.csv(path_to_antifungal_status_metadata,stringsAsFactors=F)
  
  # Assign antifungal status based on very high or low confidence for Bd inhibition.
  antifungal$Status<-ifelse(antifungal$ProbabilityOfBdInhibition <= 0.1,"NonInhibitory",ifelse(antifungal$ProbabilityOfBdInhibition >= 0.9,"Inhibitory","Uncertain"))
  
  # Create an empty storage vector for antifungal predictions.
  antifungal_predictions<-vector(mode="list",length=4)
  names(antifungal_predictions)<-c("NonInhibitory","Inhibitory","Uncertain","Other")
  
  # Locate MCMC samples files.
  print(paste0("- Locating MCMC samples files: ",elapsed.time(ptm)))
  mcmc_samples_files<-list.files(path=paste0(path_to_model_output,"/MCMC_Samples"),pattern=".csv$",full.names=T)
  
  # Begin antifungal predictions.
  print(paste0("- Beginning antifungal predictions: ",elapsed.time(ptm)))
  
  # Loop through each MCMC samples file.
  for(i in 1:length(mcmc_samples_files)){
    
    # Read in the MCMC samples file.
    NB_Reg_MCMC<-read.csv(file=mcmc_samples_files[i],check.names=F)
    
    # Get intercept terms.
    intercept_betas<-NB_Reg_MCMC[,grepl("^beta_0",colnames(NB_Reg_MCMC))]
    
    # Get stratum betas.
    stratum_betas<-NB_Reg_MCMC[,grepl("^beta_stratum",colnames(NB_Reg_MCMC))]
    
    # Get non-stratum betas.
    pred_betas<-NB_Reg_MCMC[,grepl("^beta_pred",colnames(NB_Reg_MCMC))]
    
    # Get the taxon number.
    taxon_numeric<-as.numeric(sub(pattern="^Taxon",replacement="",x=strsplit(x=basename(mcmc_samples_files[i]),split="_")[[1]][2]))
    
    # Get the taxon name.
    taxon_name<-taxa$Taxon[taxa$Taxon_number==taxon_numeric]
    
    # Get the taxon's antifungal status.
    if(taxon_name=="Other"){
      antifungal_status<-"Other"
    } else {
      antifungal_status<-antifungal$Status[antifungal$Taxon==taxon_name]
    }
    
    # Create an empty storage data frame for density predictions.
    density_predictions_storage<-as.data.frame(matrix(NA,nrow=nrow(NB_Reg_MCMC),ncol=nrow(meta)))
    colnames(density_predictions_storage)<-paste0("meta_record_",1:ncol(density_predictions_storage))
    
    # Loop through each observation.
    for(n in 1:nrow(meta)){
      
      # Get stratum predictor values for the observation.
      stratum_values<-stratum_matrix_scaled[n,]
      
      # Get non-stratum predictor values for the observation.
      pred_values<-pred_matrix_scaled[n,]
      
      # Predict density for the observation with the MCMC samples.
      density_predictions<-exp(intercept_betas+
                                 as.matrix(stratum_betas) %*% t(stratum_values)+
                                 as.matrix(pred_betas) %*% t(pred_values))
      
      # Add density predictions to the storage data frame.
      density_predictions_storage[,n]<-density_predictions
      
    }
    
    # If this is the first time an antifungal category is reached.
    if(is.null(antifungal_predictions[[antifungal_status]])){
      # Add the first data frame to the antifungal category.
      antifungal_predictions[[antifungal_status]]<-density_predictions_storage
    } else {
      # If this is not the first time an antifungal category is reached.
      # Add the data frame to the antifungal category.
      ## Get the existing antifungal category density predictions.
      density_predictions_old<-antifungal_predictions[[antifungal_status]]
      ## Add the current taxon's density predictions to the existing data frame.
      density_predictions_new<-density_predictions_old+density_predictions_storage
      ## Fix the column names of the density predictions.
      colnames(density_predictions_new)<-paste0("meta_record_",1:ncol(density_predictions_new))
      ## Add the updated density predictions to the storage list.
      antifungal_predictions[[antifungal_status]]<-density_predictions_new
    }
    
    # Report antifungal prediction progress.
    print(paste0("-- Finished antifungal predictions for taxon ",i,": ",elapsed.time(ptm)))
    
  }
  
  # Summarize antifungal predictions.
  print(paste0("- Summarizing antifungal predictions: ",elapsed.time(ptm)))
  
  # Create an emtpy storage data frame for antifungal prediction summaries.
  antifungal_summary<-data.frame(NULL)
  
  # Loop through each antifungal status.
  for(i in 1:length(antifungal_predictions)){
    # Get the antifungal status predictions.
    antifungal<-antifungal_predictions[[i]]
    # Get quantiles of the antifungal status predictions.
    antifungal_quantiles<-as.data.frame(t(apply(X=antifungal,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))
    # Format the column names.
    colnames(antifungal_quantiles)<-paste0("Quantile_",as.numeric(sub(pattern="%$",replacement="",x=colnames(antifungal_quantiles)))/100)
    # Rename the rows.
    row.names(antifungal_quantiles)<-1:nrow(antifungal_quantiles)
    # Get the antifungal status name.
    antifungal_status<-names(antifungal_predictions)[i]
    # Gather information into a data frame.
    antifungal_summary_df<-cbind(data.frame(Type="Salamander",stringsAsFactors=F),meta[,c("Date","Site","Age","LM")],data.frame(Taxon=antifungal_status,stringsAsFactors=F),antifungal_quantiles)
    # Store information in the storage data frame.
    antifungal_summary<-rbind(antifungal_summary,antifungal_summary_df)
  }
  
  # Create a field representing the numeric ordering of the antifungal status.
  antifungal_summary$Taxon_order<-ifelse(antifungal_summary$Taxon=="Uncertain",1,ifelse(antifungal_summary$Taxon=="NonInhibitory",2,ifelse(antifungal_summary$Taxon=="Inhibitory",3,4)))
  
  # Sort antifungal predictions by date, age, life stage, and taxon number.
  antifungal_summary<-antifungal_summary[order(antifungal_summary$Date,antifungal_summary$Age,antifungal_summary$LM,antifungal_summary$Taxon_order),]
  
  # Remove the taxon order field from the antifungal predictions.
  antifungal_summary<-antifungal_summary[,-which(colnames(antifungal_summary)=="Taxon_order")]
  
  # Write out summaries of antifungal predictions.
  print(paste0("- Writing out summaries of antifungal predictions: ",elapsed.time(ptm)))
  write.csv(x=antifungal_summary,file="Region16S_NB_Reg_Antifungal_Predictions.csv",row.names=F)
  
}

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
