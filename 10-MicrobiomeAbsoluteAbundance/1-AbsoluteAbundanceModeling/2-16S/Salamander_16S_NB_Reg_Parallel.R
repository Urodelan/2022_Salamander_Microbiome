########################################################################
### Fit Negative Binomial Regression Model with Salamander 16S Reads ###
########################################################################

# Use the taxon denoted by this column number in the reads matrix as the response taxon.
# (Value must be in the range of 1 to 101.)
taxon_number<-as.numeric(commandArgs(trailingOnly=T))

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Absolute/Salamander/Region16S/Take3/Parallel"

# Declare path to reads in-file (with .csv extension).
reads_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Absolute/Salamander/Region16S/Take3/Salamander_ASVTable_16S_Rate_2.csv"

# Declare barcode region ("16S" or "ITS").
barcode_region<-"16S"

# Declare path to sample metadata in-file (with .csv extension).
metadata_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/SampleMetadata2.csv"

# Declare path to JAGS model in-file (with .txt extension).
model_script_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Absolute/Model/NB_Reg_Model_JAGS.txt"

# Set model parameters.
## Number of chains.
number_of_chains<-3
## Number of adaptation iterations.
number_adaptation_iterations<-50000
## Number of warmup iterations.
number_warmup_iterations<-50000
## Number of sampling iterations.
number_sampling_iterations<-500000
## Thinning interval.
thinning_interval<-2

# Path to microbe functions script.
MicrobeFunctionsPath<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/HelperFunctions/MicrobeFunctions.R"

# Number of the most likely models to save out in the variable selection output.
# If the full variable selection output is desired, set this variable to NA.
num_most_likely_models<-10000

# Number of decimal places to report memory size of MCMC objects in gigabytes.
memory_size_decimal_places<-2

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

# Check that taxon number is valid.
if(!(taxon_number >= 1 & taxon_number <= 101)) stop("Taxon number must be in the range of 1 to 101.")

# Check that barcode region is valid.
if(!(barcode_region %in% c("16S","ITS"))) stop("User-defined barcode region is not valid.")

# Load read data.
print(paste0("- Loading read data: ",elapsed.time(ptm)))
reads<-read.csv(reads_in_file,row.names=1,check.names=F)

# Check if synthgene is present in the read data.
if(sum(colnames(reads)==paste0("Synthgene",barcode_region))!=1) stop("Synthgene not found in read data.")

# Define function for getting the most precise taxonomic level.
getMostPreciseTaxonomicLevel<-function(name){
  # Split taxon name by semi-colon.
  x<-strsplit(name,split=";")[[1]]
  # If the taxon is bacterial.
  if(x[1]=="Bacteria"){
    # Store the trailing number part.
    num<-gsub(pattern="_",replacement="",x=x[length(x)])
    # Remove the trailing number part from the name.
    x<-x[-length(x)]
    # If species is included.
    if(length(x)==7){
      # Then get both the binomial name.
      x<-paste(x[6],x[7],num)
    } else {
      # Otherwise, get the most precise taxonomic level available.
      x<-paste(x[length(x)],num)
    }
  } else { # If the taxon is fungal.
    # If species is included.
    if(length(x)==7){
      # Then get both the binomial name.
      x<-paste(x[6],x[7])
    } else {
      # Otherwise, get the most precise taxonomic level available.
      x<-x[length(x)]
    }
  }
  # Return the name.
  return(x)
}

# Get the taxon name from the reads data.
taxon_name_full<-colnames(reads)[taxon_number]

# Get the taxon's most precise taxonomic level.
taxon_name<-getMostPreciseTaxonomicLevel(name=taxon_name_full)

# Report getting taxon counts from read data.
print(paste0("- Taxon of interest: ",taxon_name))

# Report getting taxon counts from read data.
print(paste0("- Getting taxon counts from read data: ",elapsed.time(ptm)))

# Get taxon counts from the reads data.
taxon_count<-reads[,taxon_number]

# Report getting synthgene counts from read data.
print(paste0("- Getting synthgene counts from read data: ",elapsed.time(ptm)))

# Get synthgene counts from the reads data.
synthgene_count<-reads[,colnames(reads)==paste0("Synthgene",barcode_region)]

# Check if any synthgene counts are zero.
if(any(synthgene_count==0)) stop("At least one synthgene count is zero.")

# Report preparing predictors.
print(paste0("- Preparing predictors: ",elapsed.time(ptm)))

# Load metadata.
metadata<-read.csv(metadata_in_file)

# Subset metadata to just samples included in the reads data.
metadata<-metadata[metadata$SampleID %in% row.names(reads),]

# Order metadata to match read data.
metadata<-metadata[match(metadata$SampleID,row.names(reads)),]

# Check that metadata and read data samples match.
if(!identical(as.character(metadata$SampleID),row.names(reads))) stop("The samples in the read and metadata files are not the same.")

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

# Report estimating swabbed area.
print(paste0("- Estimating swabbed area: ",elapsed.time(ptm)))

# Create a field in the metadata for estimated swabbed area using
# an equation derived from length-weight regression.
metadata$Area<-metadata$SVL_mm^1.675231

# Report preparing stratum predictors.
print(paste0("- Preparing stratum predictors: ",elapsed.time(ptm)))

# Get stratum predictor matrix.
stratum_matrix<-as.data.frame(model.matrix(~Stratum-1,data=metadata))

# Get stratum scaling parameters.
stratum_scaling<-data.frame(Mean=apply(X=stratum_matrix,MARGIN=2,FUN=mean))

# Copy stratum matrix for scaling.
stratum_matrix_scaled<-stratum_matrix

# Center stratum predictors (no scaling by SD).
for(i in 1:ncol(stratum_matrix)){
  stratum_matrix_scaled[,i]<-stratum_matrix[,i]-stratum_scaling$Mean[i]
}

# Report preparing non-stratum predictors.
print(paste0("- Preparing non-stratum predictors: ",elapsed.time(ptm)))

# Convert date to number of weeks since June 9th, 2018.
metadata$Week<-as.vector(metadata$Date-as.Date("2018-06-09"))/7

# Get site as a binary variable.
metadata$Site_binary<-as.numeric(metadata$Site)-1
## Gibson: 0
## Ponds: 1

# Get age as a continuous variable.
metadata$Age_numeric<-as.numeric(metadata$Age)-1
## Age-0: 0
## Age-1: 1
## Age-2+: 2

# Get LM as a binary variable.
metadata$LM_binary<-as.numeric(metadata$LM)-1
## Larvae/Neotene: 0
## Metamorphosed: 1

# Get just age, LM, site, and week.
pred_matrix<-metadata[,c("Age_numeric","LM_binary","Site_binary","Week")]

# Get mean and SD of age.
age_mean<-mean(pred_matrix$Age_numeric)
age_sd<-sd(pred_matrix$Age_numeric)

# Get mean of binary LM.
LM_mean<-mean(pred_matrix$LM_binary)

# Get mean of binary site.
site_mean<-mean(pred_matrix$Site_binary)

# Get mean and SD of week.
week_mean<-mean(pred_matrix$Week)
week_sd<-sd(pred_matrix$Week)

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

# Collect non-stratum predictor scaling parameters into a data frame.
pred_scaling<-data.frame(Predictor=c("Age","Age","LM","Site","Week","Week"),
                         Parameter=c("Mean","SD","Mean","Mean","Mean","SD"),
                         Value=c(age_mean,age_sd,LM_mean,site_mean,week_mean,week_sd))

# Write out scaling parameters.
print(paste0("- Writing out scaling parameters: ",elapsed.time(ptm)))
write.csv(stratum_scaling,file=paste0("Stratum_Scaling/Region",barcode_region,"_Taxon",taxon_number,"_StratumScaling.csv"),row.names=T)
write.csv(pred_scaling,file=paste0("Predictor_Scaling/Region",barcode_region,"_Taxon",taxon_number,"_PredictorScaling.csv"),row.names=F)

# Report defining model data.
print(paste0("- Defining model data: ",elapsed.time(ptm)))

# Define model data.
NB_Reg_Data<-list(
  NSamples=length(taxon_count), # Number of samples.
  NStrata=ncol(stratum_matrix_scaled), # Number of strata.
  StratumMatrix=stratum_matrix_scaled, # Predictor matrix for strata.
  NPredictors=ncol(pred_matrix_scaled), # Number of non-strata predictors.
  PredictorMatrix=pred_matrix_scaled, # Non-strata predictor matrix.
  TaxonCount=taxon_count, # Taxon count.
  SynthgeneCount=synthgene_count, # Synthgene count.
  Area=metadata$Area # Estimated swabbed area.
)

# Load rjags.
print(paste0("- Loading rjags: ",elapsed.time(ptm)))
library(rjags)

# Initialize the model.
print(paste0("- Initializing model: ",elapsed.time(ptm)))
NB_Reg<-jags.model(file=model_script_in_file,
                   data=NB_Reg_Data,
                   n.adapt=number_adaptation_iterations,
                   n.chains=3)

# Burn-in the model.
print(paste0("- Burning-in model: ",elapsed.time(ptm)))
update(NB_Reg,n.iter=number_warmup_iterations)

# Sample the model.
print(paste0("- Sampling model: ",elapsed.time(ptm)))
NB_Reg_Out<-coda.samples(model=NB_Reg,
                         variable.names=c("beta_0","beta_stratum","beta_pred","r",
                                          "tau_squared_stratum","tau_squared_pred","lambda_squared",
                                          "IV_stratum","IV_pred","IP"),
                         n.iter=number_sampling_iterations,
                         thin=thinning_interval)

# Get convergence summaries.
print(paste0("- Getting convergence summaries: ",elapsed.time(ptm)))
gelman_diag_values<-gelman.diag(NB_Reg_Out,multivariate=F)$psrf

# Check convergence statistics.
check_gelman_diag_values<-sum(gelman_diag_values[,"Upper C.I."] > 1.1,na.rm=T)

# Write out model Gelman convergence diagnostics.
print(paste0("- Writing out convergence summaries: ",elapsed.time(ptm)))
write.csv(gelman_diag_values,file=paste0("Convergence/Region",barcode_region,"_Taxon",taxon_number,"_Convergence.csv"),row.names=T)

# Print information on Gelman convergence diagnostics.
print(paste0("- Gelman convergence diagnostics greater than 1.1: ",check_gelman_diag_values," of ",nrow(gelman_diag_values)))

# # Traceplots.
# par(mar=c(1,1,1,1))
# plot(NB_Reg_Out,ask=T)
# par(mar=c(5.1,4.1,4.1,2.1))

# Print the memory size of MCMC output.
print(paste0("- Memory size of MCMC output: ",
             as.character(format(object.size(NB_Reg_Out),
                                 units="GB",standard="SI",
                                 digits=memory_size_decimal_places))))

# Combine all MCMC estimates.
print(paste0("- Getting MCMC samples: ",elapsed.time(ptm)))
NB_Reg_MCMC<-as.data.frame(do.call("rbind",NB_Reg_Out))

# Get summary statistics of MCMC samples.
print(paste0("- Getting MCMC sample quantiles: ",elapsed.time(ptm)))
NB_Reg_MCMC_quantiles<-as.data.frame(t(apply(NB_Reg_MCMC,MARGIN=2,FUN=quantile,probs=c(0.025,0.25,0.5,0.75,0.975))))

# Write out summary statistics of MCMC samples.
print(paste0("- Writing out MCMC sample quantiles: ",elapsed.time(ptm)))
write.csv(NB_Reg_MCMC_quantiles,file=paste0("MCMC_Quantiles/Region",barcode_region,"_Taxon",taxon_number,"_MCMCQuantiles.csv"),row.names=T)

# Get just the binary inclusion variables from the MCMC output.
print(paste0("- Getting binary inclusion variables from MCMC output: ",elapsed.time(ptm)))
NB_Reg_IVs<-NB_Reg_MCMC[,grepl(pattern="^IV",x=colnames(NB_Reg_MCMC))]

# Calculate posterior inclusion probabilities for the predictors.
print(paste0("- Calculating posterior inclusion probabilities: ",elapsed.time(ptm)))
PIP<-as.data.frame(apply(X=NB_Reg_IVs,MARGIN=2,FUN=mean))
colnames(PIP)<-"PIP"
row.names(PIP)<-c(colnames(pred_matrix_scaled),"Stratum")

# Write out posterior inclusion probabilities.
print(paste0("- Writing out posterior inclusion probabilities: ",elapsed.time(ptm)))
write.csv(PIP,file=paste0("PIP/Region",barcode_region,"_Taxon",taxon_number,"_PIP.csv"),row.names=T)

# Get summary of counts by combination of variables.
NB_Reg_IVs$Count<-1
selections<-aggregate(Count~.,data=NB_Reg_IVs,FUN=sum)
colnames(selections)[1:(ncol(selections)-1)]<-c(colnames(pred_matrix_scaled),"Stratum")

# Add field for proportion of iterations.
selections$Proportion<-selections$Count/sum(selections$Count)

# Sort selection data frame by the counts.
selections<-selections[order(selections$Count,decreasing=T),]

# If the number of most likely models requested is not NA.
if(!is.na(num_most_likely_models)){
  # If the number of records in the selection data frame exceeds the
  # number of most likely models requested, then subset the data frame
  # to the number of most likely models requested.
  if(nrow(selections) > num_most_likely_models){
    selections<-selections[1:num_most_likely_models,]
  }
}

# Write out variable selection summary.
print(paste0("- Writing out variable selection summary: ",elapsed.time(ptm)))
write.csv(selections,file=paste0("Variable_Selection/Region",barcode_region,"_Taxon",taxon_number,"_VariableSelection.csv"),row.names=F)

# Report preparing predictors for prediction.
print(paste0("- Preparing predictors for prediction: ",elapsed.time(ptm)))

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

# Get intercept terms.
intercept_betas<-NB_Reg_MCMC[,grepl("^beta_0",colnames(NB_Reg_MCMC))]

# Get stratum betas.
stratum_betas<-NB_Reg_MCMC[,grepl("^beta_stratum",colnames(NB_Reg_MCMC))]

# Get non-stratum betas.
pred_betas<-NB_Reg_MCMC[,grepl("^beta_pred",colnames(NB_Reg_MCMC))]

# If the barcode region is 16S.
if(barcode_region=="16S"){
  
  # Get all MCMC samples for regression coefficients in a data frame.
  MCMC_samples<-cbind(data.frame(beta_0=intercept_betas),stratum_betas,pred_betas)
  
  # Write out MCMC samples for regression coefficients.
  print(paste0("- Writing out regression coefficient MCMC samples: ",elapsed.time(ptm)))
  write.csv(MCMC_samples,file=paste0("MCMC_Samples/Region16S_Taxon",taxon_number,"_MCMCSamples.csv"),row.names=F)
  
  # Print the memory size of regression coefficient MCMC samples.
  print(paste0("- Memory size of regression coefficient MCMC samples: ",
               as.character(format(object.size(MCMC_samples),
                                   units="GB",standard="SI",
                                   digits=memory_size_decimal_places))))
  
}

# Report beginning predictions.
print(paste0("- Predicting densities: ",elapsed.time(ptm)))

# Create empty storage data frame for density predictions.
density_predictions_df<-data.frame(NULL)

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
  
  # Get quantiles of microbe density.
  density_prediction_quantiles<-quantile(x=density_predictions,probs=c(0.025,0.25,0.5,0.75,0.975))
  
  # Rename elements of the density quantile vector.
  names(density_prediction_quantiles)<-paste0("Quantile_",as.numeric(gsub(pattern="%",replacement="",x=names(density_prediction_quantiles)))/100)
  
  # Create a data frame containing density prediction quantiles and metadata.
  density_prediction_quantiles_df<-cbind(data.frame(Type="Salamander"),
                                         meta[n,c("Date","Site","Age","LM")],
                                         data.frame(Taxon=taxon_name_full,stringsAsFactors=F),
                                         as.data.frame(t(density_prediction_quantiles)))
  
  # Add the density prediction quantiles and metadata to the storage data frame.
  density_predictions_df<-rbind(density_predictions_df,density_prediction_quantiles_df)
  
}

# Write out density predictions.
print(paste0("- Writing out density predictions: ",elapsed.time(ptm)))
write.csv(density_predictions_df,file=paste0("Density_Predictions/Region",barcode_region,"_Taxon",taxon_number,"_DensityPredictions.csv"),row.names=F)

# Done!
print(paste0("Done! (",elapsed.time(ptm),")"))
