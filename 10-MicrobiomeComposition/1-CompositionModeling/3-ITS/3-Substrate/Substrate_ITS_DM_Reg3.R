###################################################################################
### Fit Primary Dirichlet-Multinomial Regression Model with Substrate ITS Reads ###
###################################################################################

##############################
### User-Defined Variables ###
##############################

# Declare working directory.
working_directory<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition/Substrate/RegionITS/Take3"

# Declare path to reads in-file (with .csv extension).
reads_in_file<-"Substrate_ASVTable_ITS_Comp_2.csv"

# Declare path to sample metadata in-file (with .csv extension).
metadata_in_file<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Metadata/SampleMetadata2.csv"

# Declare path to file containing model Stan script in-file (with .stan extension).
# This model script includes both stratum and non-stratum predictors.
model_script_in_file_with_both<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition/Model/PrimaryCompositionModel.stan"

# Declare path to file containing model Stan script in-file (with .stan extension).
# This model script does not include stratum as a predictor.
model_script_in_file_without_stratum<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition/Model/PrimaryCompositionModel_NS.stan"

# Declare path to file containing model Stan script in-file (with .stan extension).
# This model script does not include non-stratum predictors.
model_script_in_file_without_nonstratum<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition/Model/PrimaryCompositionModel_NN.stan"

# Declare path to file containing model Stan script in-file (with .stan extension).
# This is an intercept only model.
model_script_in_file_intercept_only<-"/uufs/chpc.utah.edu/common/home/gompert-group2/data/kg_microbes/Scripts/Models/Composition/Model/PrimaryCompositionModel_IO.stan"

# Set model parameters.
## Number of chains.
number_of_chains<-4
## Number of cores.
number_of_cores<-4
## Number of warmup iterations.
number_warmup_iterations<-500
## Number of sampling iterations.
number_sampling_iterations<-500
## Thinning interval.
thinning_interval<-1
## Adapt delta argument.
adapt_delta_argument<-0.95
## Maximum tree depth argument.
max_treedepth_argument<-20
## Set seed for Stan.
stan_seed<-1234

## Path to microbe functions script.
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

# Load rstan.
print("- Loading packages.")
suppressPackageStartupMessages(library(rstan))

# Load Dirichlet-multinomial functions.
library(extraDistr)

# Load read data.
print("- Loading read data.")
reads<-read.csv(reads_in_file,row.names=1,check.names=F)

# Add one to the read data so that there are no zeros.
reads<-reads+1

# Report preparing predictors.
print("- Preparing predictors.")

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

# Preparing stratum predictors.

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

# Prepare non-stratum predictors.

# Convert date to number of weeks since June 9th, 2018.
metadata$Week<-as.vector(metadata$Date-as.Date("2018-06-09"))/7

# Get site as a binary variable.
metadata$Site_binary<-as.numeric(metadata$Site)-1
## Gibson: 0
## Ponds: 1

# Get just site and week.
pred_matrix<-metadata[,c("Site_binary","Week")]

# Get mean of binary site.
site_mean<-mean(pred_matrix$Site_binary)

# Get mean and SD of week.
week_mean<-mean(pred_matrix$Week)
week_sd<-sd(pred_matrix$Week)

# Scale the predictor matrix.
pred_matrix_scaled<-data.frame(Site=pred_matrix$Site_binary-site_mean,
                               Week=(pred_matrix$Week-week_mean)/week_sd)

# Add interaction terms to the scaled predictor matrix.
pred_matrix_scaled<-as.data.frame(model.matrix(~Site*poly(Week,degree=2,raw=T)-1,data=pred_matrix_scaled))

# Call the polynomial for week 'Week1' and 'Week2' in the scaled predictor matrix.
colnames(pred_matrix_scaled)<-gsub(pattern="poly\\(Week,(.*)\\)",replacement="Week",x=colnames(pred_matrix_scaled))

# Collect non-stratum predictor scaling parameters into a data frame.
pred_scaling<-data.frame(Predictor=c("Site","Week","Week"),
                         Parameter=c("Mean","Mean","SD"),
                         Value=c(site_mean,week_mean,week_sd))

# Write out scaling parameters.
print("- Writing out scaling parameters.")
write.csv(stratum_scaling,file="stratum_scaling.csv",row.names=T)
write.csv(pred_scaling,file="pred_scaling.csv",row.names=F)

# Report beginning initial model.
print("- Beginning initial model.")

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# Define model parameters.
data_Model<-list(
  "NSamples"=nrow(reads), # Number of samples.
  "NTaxa"=ncol(reads), # Number of taxa.
  "NStrata"=ncol(stratum_matrix_scaled), # Number of strata.
  "StratumMatrix"=stratum_matrix_scaled, # Predictor matrix for strata.
  "NPredictors"=ncol(pred_matrix_scaled), # Number of non-strata predictors.
  "PredictorMatrix"=pred_matrix_scaled, # Non-strata predictor matrix.
  "ReadsMatrix"=reads, # Reads response matrix.
  "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
)

# Compile Stan program into C++ code.
DM_Reg<-stan_model(file=model_script_in_file_with_both,model_name="DM_Reg")

# Run the model with HMC sampling.
mod_fit<-sampling(DM_Reg,
                  data=data_Model,
                  chains=number_of_chains,
                  warmup=number_warmup_iterations,
                  iter=number_warmup_iterations+number_sampling_iterations,
                  thin=thinning_interval,
                  seed=stan_seed,
                  algorithm="NUTS",
                  cores=number_of_cores,
                  pars<-c("beta_0","beta_stratum","beta_pred",
                          "sigma_stratum","theta"),
                  verbose=F,
                  control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                  refresh=0
                  )

# Get convergence statistics.
convergence<-as.data.frame(summary(mod_fit)$summary[,c("n_eff","Rhat")])

# Check convergence statistics.
check_n_eff<-sum(convergence$n_eff < 100,na.rm=T)
check_Rhat<-sum(convergence$Rhat > 1.05,na.rm=T)

# Get HMC output.
HMC<-as.data.frame(mod_fit)

# Get intercept terms.
intercept_betas<-HMC[,grepl("^beta_0",colnames(HMC))]

# Get stratum betas.
stratum_betas<-HMC[,grepl("^beta_stratum",colnames(HMC))]

# Get non-stratum betas.
pred_betas<-HMC[,grepl("^beta_pred",colnames(HMC))]

# Get theta estimates.
theta<-HMC[,colnames(HMC)=="theta"]

# Get exponentiated theta estimates.
exptheta<-exp(theta)

# Create storage data frame for probability of observing data given the estimated parameters.
log_p<-data.frame(Field_Sample=rep(metadata$SampleID,each=nrow(HMC)),
                  HMC_Sample=rep(1:nrow(HMC),nrow(metadata)),
                  log_P_obs=rep(NA,nrow(metadata)*nrow(HMC)))

# Loop through each observation.
for(i in 1:nrow(reads)){
  
  # Get stratum predictor values for the observation.
  stratum_values<-stratum_matrix_scaled[i,]
  
  # Get non-stratum predictor values for the observation.
  pred_values<-pred_matrix_scaled[i,]
  
  # Get microbe count data for the observation.
  read_values<-reads[i,]
  
  # Create empty storage data frame for linear combination of HMC estimates.
  eta<-as.data.frame(matrix(data=NA,nrow=nrow(HMC),ncol=ncol(reads)))
  
  # Loop through each taxon.
  for(k in 1:ncol(reads)){
    
    # Get stratum betas for the taxon.
    stratum_taxon_betas<-stratum_betas[,grepl(paste0("\\[",k,","),colnames(stratum_betas))]
    
    # Get non-stratum betas for the taxon.
    pred_taxon_betas<-pred_betas[,grepl(paste0("\\[",k,","),colnames(pred_betas))]
    
    # Calculate linear combination of HMC estimates.
    eta[,k]<-intercept_betas[,k]+
      as.matrix(stratum_taxon_betas) %*% t(stratum_values)+
      as.matrix(pred_taxon_betas) %*% t(pred_values)
    
  }
  
  # Get expected microbial proportions.
  pi<-as.data.frame(t(apply(X=eta,MARGIN=1,FUN=softmax)))
  
  # Loop through each HMC estimate.
  for(j in 1:nrow(pi)){
    
    # Get expected microbial proprotions for the HMC estimate.
    pi_estimates<-unlist(pi[j,])
    
    # Get the probability of observing this particular set of microbe counts
    # given the Dirichlet precision and expected microbial proportions.
    log_prob_obs<-ddirmnom(x=read_values,size=sum(read_values),alpha=pi_estimates*exptheta[j],log=T)
    
    # Store this probability in the storage data frame.
    log_p$log_P_obs[(i-1)*nrow(HMC)+j]<-log_prob_obs
    
  }
  
}

# Exponentiate log(p).
log_p$P_obs<-exp(log_p$log_P_obs)

# Check that none of the probabilities are 0.
if(any(log_p$P_obs==0)) stop("Probabilities of 0 were produced during WAIC calculation.")

# Calculate computed log pointwise predictive density (lppd).
## Take the average of HMC sample probability masses for each field sample.
lppd_calc<-aggregate(P_obs~Field_Sample,data=log_p,FUN=mean)
## Take the log of the mean probability mass for each field sample.
lppd_calc$log_P_obs<-log(lppd_calc$P_obs)
## Calculate lppd by summing the log mean probability mass for each field sample.
lppd<-sum(lppd_calc$log_P_obs)

# Calculate the second version of the WAIC bias correction.
## Take the variance of HMC sample log probability masses for each field sample.
P_WAIC2_calc<-aggregate(log_P_obs~Field_Sample,data=log_p,FUN=var)
## Calculate the WAIC bias correction by summing the variance of log probability
## masses for each field sample.
P_WAIC2<-sum(P_WAIC2_calc$log_P_obs)

# Calculate the expected log pointwise predictive density for a new dataset
# (elppd) by subtracting the bias correction from the lppd.
elppd<-lppd-P_WAIC2

# Multiply the elppd by -2 to obtain WAIC.
WAIC<-elppd*-2

# Create a data frame for monitoring the variable selection process.
## Create data frame.
df<-as.data.frame(matrix(ncol=3+ncol(pred_matrix_scaled)+5))
colnames(df)<-c("Set","Model","Stratum",colnames(pred_matrix_scaled),
                "Num_Rhats_gt_1.05","Num_Eff_SS_ls_100",
                "Num_params","Run_Time","WAIC")
## Add values to the data frame.
### Provide the set number.
df$Set<-1
### Provide a model number.
df$Model<-1
### State wether stratum was included as a predictor.
df$Stratum<-1
### Check to see which other predictors were used in the model.
for(i in 4:(3+ncol(pred_matrix_scaled))){
  df[,i]<-as.numeric(colnames(df)[i] %in% colnames(pred_matrix_scaled))
}
### Provide number of R-hats greater than 1.05.
df$Num_Rhats_gt_1.05<-check_Rhat
### Provide number of effective sample sizes less than 100.
df$Num_Eff_SS_ls_100<-check_n_eff
### Provide the number of parameters in the model.
df$Num_params<-nrow(convergence)
### Provide model run time (including time for WAIC calculations).
df$Run_Time<-elapsed.time(ptm)
### Provide WAIC.
df$WAIC<-WAIC

# Report finished with initial model.
print(paste0("- Finished fitting initial model: Run time of ",elapsed.time(ptm)))

# Write out model selection summary.
write.csv(df,file="model_selection.csv",row.names=F)
print("- Initial model written to model selection file.")

# Create variable denoting wether WAIC has been minimized.
WAIC_minimized<-"No"

# Create a model set number to iterate, beginning at 2.
model_set<-2

# Set full model as the initial model to compare to.
best_fit_model_from_last_set<-df

# Report beginning backwards variable selection.
print("- Beginning backwards variable selection.")

# Generate combinations of predictors to fit the next set of models with, removing a single
# variable from the model at a time to minimize WAIC.
## Get the names of predictors in the initial model.
initial_set_of_predictors<-c("Stratum",colnames(pred_matrix_scaled))
## Get new predictor combinations to fit models to.
pred_comb<-combn(x=initial_set_of_predictors,
                 m=length(initial_set_of_predictors)-1,
                 simplify=F)

# Begin backwards variable selection. Remove variables sequentially
# until WAIC is minimized.
while(WAIC_minimized=="No"){
  
  # Loop through each predictor combination to fit models to.
  for(z in 1:length(pred_comb)){
    
    # Get predictor combination.
    new_pred<-pred_comb[[z]]
    
    # Check if stratum is included in the predictor combination.
    stratum_included<-ifelse("Stratum" %in% new_pred,"Yes","No")
    
    # Check if any non-stratum predictors are included in the predictor combination.
    nonstratum_included<-ifelse(any(new_pred!="Stratum",na.rm=T),"Yes","No")
    
    # If non-stratum predictors are included, get new non-stratum predictor matrix
    # with the predictor combination.
    if(nonstratum_included=="Yes"){
      pred_matrix_scaled_reduced<-pred_matrix_scaled[,colnames(pred_matrix_scaled) %in% new_pred]
    }
    
    # Store results of initial proc.time for timing purposes.
    ptm<-proc.time()
    
    # If included as predictors, use model with both stratum and non-stratum predictors.
    if(stratum_included=="Yes" & nonstratum_included=="Yes"){
      
      # Define model parameters.
      data_Model<-list(
        "NSamples"=nrow(reads), # Number of samples.
        "NTaxa"=ncol(reads), # Number of taxa.
        "NStrata"=ncol(stratum_matrix_scaled), # Number of strata.
        "StratumMatrix"=stratum_matrix_scaled, # Predictor matrix for strata.
        "NPredictors"=ncol(pred_matrix_scaled_reduced), # Number of non-strata predictors.
        "PredictorMatrix"=pred_matrix_scaled_reduced, # Non-strata predictor matrix.
        "ReadsMatrix"=reads, # Reads response matrix.
        "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
      )
      
      # Compile Stan program into C++ code.
      DM_Reg<-stan_model(file=model_script_in_file_with_both,model_name="DM_Reg")
      
      # Run the model with HMC sampling.
      mod_fit<-sampling(DM_Reg,
                        data=data_Model,
                        chains=number_of_chains,
                        warmup=number_warmup_iterations,
                        iter=number_warmup_iterations+number_sampling_iterations,
                        thin=thinning_interval,
                        seed=stan_seed,
                        algorithm="NUTS",
                        cores=number_of_cores,
                        pars<-c("beta_0","beta_stratum","beta_pred",
                                "sigma_stratum","theta"),
                        verbose=F,
                        control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                        refresh=0
      )
      
    } else if(stratum_included=="No" & nonstratum_included=="Yes"){
      
      # If stratum is not included but non-stratum predictors are,
      # use model without stratum and with non-stratum predictors.
      
      # Define model parameters.
      data_Model<-list(
        "NSamples"=nrow(reads), # Number of samples.
        "NTaxa"=ncol(reads), # Number of taxa.
        "NPredictors"=ncol(pred_matrix_scaled_reduced), # Number of non-strata predictors.
        "PredictorMatrix"=pred_matrix_scaled_reduced, # Non-strata predictor matrix.
        "ReadsMatrix"=reads, # Reads response matrix.
        "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
      )
      
      # Compile Stan program into C++ code.
      DM_Reg<-stan_model(file=model_script_in_file_without_stratum,model_name="DM_Reg")
      
      # Run the model with HMC sampling.
      mod_fit<-sampling(DM_Reg,
                        data=data_Model,
                        chains=number_of_chains,
                        warmup=number_warmup_iterations,
                        iter=number_warmup_iterations+number_sampling_iterations,
                        thin=thinning_interval,
                        seed=stan_seed,
                        algorithm="NUTS",
                        cores=number_of_cores,
                        pars<-c("beta_0","beta_pred","theta"),
                        verbose=F,
                        control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                        refresh=0
      )
      
    } else if(stratum_included=="Yes" & nonstratum_included=="No"){
      
      # If stratum is included but non-stratum predictors are not,
      # use model with stratum and without non-stratum predictors.
      
      # Define model parameters.
      data_Model<-list(
        "NSamples"=nrow(reads), # Number of samples.
        "NTaxa"=ncol(reads), # Number of taxa.
        "NStrata"=ncol(stratum_matrix_scaled), # Number of strata.
        "StratumMatrix"=stratum_matrix_scaled, # Predictor matrix for strata.
        "ReadsMatrix"=reads, # Reads response matrix.
        "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
      )
      
      # Compile Stan program into C++ code.
      DM_Reg<-stan_model(file=model_script_in_file_without_nonstratum,model_name="DM_Reg")
      
      # Run the model with HMC sampling.
      mod_fit<-sampling(DM_Reg,
                        data=data_Model,
                        chains=number_of_chains,
                        warmup=number_warmup_iterations,
                        iter=number_warmup_iterations+number_sampling_iterations,
                        thin=thinning_interval,
                        seed=stan_seed,
                        algorithm="NUTS",
                        cores=number_of_cores,
                        pars<-c("beta_0","beta_stratum","sigma_stratum","theta"),
                        verbose=F,
                        control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                        refresh=0
      )
      
    } else {
      
      # If neither stratum nor non-stratum predictors are included,
      # use the intercept only model.
      
      # Define model parameters.
      data_Model<-list(
        "NSamples"=nrow(reads), # Number of samples.
        "NTaxa"=ncol(reads), # Number of taxa.
        "ReadsMatrix"=reads, # Reads response matrix.
        "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
      )
      
      # Compile Stan program into C++ code.
      DM_Reg<-stan_model(file=model_script_in_file_intercept_only,model_name="DM_Reg")
      
      # Run the model with HMC sampling.
      mod_fit<-sampling(DM_Reg,
                        data=data_Model,
                        chains=number_of_chains,
                        warmup=number_warmup_iterations,
                        iter=number_warmup_iterations+number_sampling_iterations,
                        thin=thinning_interval,
                        seed=stan_seed,
                        algorithm="NUTS",
                        cores=number_of_cores,
                        pars<-c("beta_0","theta"),
                        verbose=F,
                        control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                        refresh=0
      )
      
    }
    
    # Get convergence statistics.
    convergence<-as.data.frame(summary(mod_fit)$summary[,c("n_eff","Rhat")])
    
    # Check convergence statistics.
    check_n_eff<-sum(convergence$n_eff < 100,na.rm=T)
    check_Rhat<-sum(convergence$Rhat > 1.05,na.rm=T)
    
    # Get HMC output.
    HMC<-as.data.frame(mod_fit)
    
    # Get intercept terms.
    intercept_betas<-HMC[,grepl("^beta_0",colnames(HMC))]
    
    # Get stratum betas if stratum is included as a predictor.
    if(stratum_included=="Yes") stratum_betas<-HMC[,grepl("^beta_stratum",colnames(HMC))]
    
    # Get non-stratum betas if included as predictors.
    if(nonstratum_included=="Yes") pred_betas<-HMC[,grepl("^beta_pred",colnames(HMC))]
    
    # Get theta estimates.
    theta<-HMC[,colnames(HMC)=="theta"]
    
    # Get exponentiated theta estimates.
    exptheta<-exp(theta)
    
    # Create storage data frame for probability of observing data given the estimated parameters.
    log_p<-data.frame(Field_Sample=rep(metadata$SampleID,each=nrow(HMC)),
                      HMC_Sample=rep(1:nrow(HMC),nrow(metadata)),
                      log_P_obs=rep(NA,nrow(metadata)*nrow(HMC)))
    
    # Loop through each observation.
    for(i in 1:nrow(reads)){
      
      # Get stratum predictor values for the observation if stratum is included as a predictor.
      if(stratum_included=="Yes") stratum_values<-stratum_matrix_scaled[i,]
      
      # Get non-stratum predictor values for the observation if included as predictors.
      if(nonstratum_included=="Yes") pred_values<-pred_matrix_scaled_reduced[i,]
      
      # Get microbe count data for the observation.
      read_values<-reads[i,]
      
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
      
      # Loop through each HMC estimate.
      for(j in 1:nrow(pi)){
        
        # Get expected microbial proprotions for the HMC estimate.
        pi_estimates<-unlist(pi[j,])
        
        # Get the probability of observing this particular set of microbe counts
        # given the Dirichlet precision and expected microbial proportions.
        log_prob_obs<-ddirmnom(x=read_values,size=sum(read_values),alpha=pi_estimates*exptheta[j],log=T)
        
        # Store this probability in the storage data frame.
        log_p$log_P_obs[(i-1)*nrow(HMC)+j]<-log_prob_obs
        
      }
      
    }
    
    # Exponentiate log(p).
    log_p$P_obs<-exp(log_p$log_P_obs)
    
    # Check that none of the probabilities are 0.
    if(any(log_p$P_obs==0)) stop("Probabilities of 0 were produced during WAIC calculation.")
    
    # Calculate computed log pointwise predictive density (lppd).
    ## Take the average of HMC sample probability masses for each field sample.
    lppd_calc<-aggregate(P_obs~Field_Sample,data=log_p,FUN=mean)
    ## Take the log of the mean probability mass for each field sample.
    lppd_calc$log_P_obs<-log(lppd_calc$P_obs)
    ## Calculate lppd by summing the log mean probability mass for each field sample.
    lppd<-sum(lppd_calc$log_P_obs)
    
    # Calculate the second version of the WAIC bias correction.
    ## Take the variance of HMC sample log probability masses for each field sample.
    P_WAIC2_calc<-aggregate(log_P_obs~Field_Sample,data=log_p,FUN=var)
    ## Calculate the WAIC bias correction by summing the variance of log probability
    ## masses for each field sample.
    P_WAIC2<-sum(P_WAIC2_calc$log_P_obs)
    
    # Calculate the expected log pointwise predictive density for a new dataset
    # (elppd) by subtracting the bias correction from the lppd.
    elppd<-lppd-P_WAIC2
    
    # Multiply the elppd by -2 to obtain WAIC.
    WAIC<-elppd*-2
    
    # Create a data frame for monitoring the variable selection process.
    ## Create data frame.
    df_new_mod<-as.data.frame(matrix(ncol=3+ncol(pred_matrix_scaled)+5))
    colnames(df_new_mod)<-c("Set","Model","Stratum",colnames(pred_matrix_scaled),
                    "Num_Rhats_gt_1.05","Num_Eff_SS_ls_100",
                    "Num_params","Run_Time","WAIC")
    ## Add values to the data frame.
    ### Provide the set number.
    df_new_mod$Set<-model_set
    ### Provide a model number.
    df_new_mod$Model<-z
    ### State wether stratum was included as a predictor.
    df_new_mod$Stratum<-as.numeric(stratum_included=="Yes")
    ### If non-stratum predictors were used in the model, check to see which.
    if(nonstratum_included=="Yes"){
      for(i in 4:(3+ncol(pred_matrix_scaled))){
        df_new_mod[,i]<-as.numeric(colnames(df_new_mod)[i] %in% colnames(pred_matrix_scaled_reduced))
      }
    } else {
      df_new_mod[,4:(3+ncol(pred_matrix_scaled))]<-0
    }
    ### Provide number of R-hats greater than 1.05.
    df_new_mod$Num_Rhats_gt_1.05<-check_Rhat
    ### Provide number of effective sample sizes less than 100.
    df_new_mod$Num_Eff_SS_ls_100<-check_n_eff
    ### Provide the number of parameters in the model.
    df_new_mod$Num_params<-nrow(convergence)
    ### Provide model run time (including time for WAIC calculations).
    df_new_mod$Run_Time<-elapsed.time(ptm)
    ### Provide WAIC.
    df_new_mod$WAIC<-WAIC
    
    # Append new model information to cumulative model information data frame.
    df<-rbind(df,df_new_mod)
    
    # Add model information to file.
    write.csv(df,file="model_selection.csv",row.names=F)
    print(paste0("-- Model added to file: Set ",model_set,
                 ", Model ",z," of ",length(pred_comb),
                 ", Run time of ",elapsed.time(ptm)))
    
  }
  
  # Get WAIC value from last set's best-fit model.
  best_fit_model_from_last_set_WAIC<-best_fit_model_from_last_set$WAIC
  
  # Get current model set.
  current_model_set<-df[df$Set==model_set,]
  
  # Get the best-fitting model from the current model set. If multiple models have the lowest WAIC,
  # then use the first of these models.
  best_fit_model_from_current_set<-
    current_model_set[which(current_model_set$WAIC==min(current_model_set$WAIC))[1],]
  
  # Get WAIC value from the current set's best-fit model.
  best_fit_model_from_current_set_WAIC<-best_fit_model_from_current_set$WAIC
  
  ## Get predictors of best-fitting model from the current model set.
  best_fit_model_from_current_set_predictors<-
    best_fit_model_from_current_set[,3:(3+ncol(pred_matrix_scaled))]
  
  # If the best WAIC from the last model set is lower than the best WAIC from
  # the current model set, or the last model set included no stratum nor non-stratum
  # predictors, then set the WAIC minimized variable to 'Yes' and exit the while loop.
  if((best_fit_model_from_last_set_WAIC < best_fit_model_from_current_set_WAIC) |
     (sum(best_fit_model_from_current_set_predictors)==0)){
    
    # Set the WAIC minimized variable to 'Yes'.
    WAIC_minimized<-"Yes"
    
    # If the best-fit model from the last set fits better than the
    # best-fit model from the current set, then select the model
    # from the last set.
    if(best_fit_model_from_last_set_WAIC < best_fit_model_from_current_set_WAIC){
      # Display selected model at base of model selection data frame.
      best_fit_model_from_last_set$Set<-"End"
      best_fit_model_from_last_set$Model<-"Best"
      df<-rbind(df,best_fit_model_from_last_set)
    } else {
      # If the best-fit model from the current set fits as well or better than
      # the best-fit model from the last set, and the current set uses
      # an intercept-only model, then select the model from the current set and
      # display the selected model at base of model selection data frame.
      # The while loop ends since there are no more predictors to remove from
      # the intercept-only model.
      best_fit_model_from_current_set$Set<-"End"
      best_fit_model_from_current_set$Model<-"Best"
      df<-rbind(df,best_fit_model_from_current_set)
    }
    
    # Add selected model information to file.
    write.csv(df,file="model_selection.csv",row.names=F)
    print("-- A model has been selected and added to the file!")
    
  } else {
    
    # Generate combinations of predictors to fit the next set of models with, removing a single
    # variable from the model with the lowest WAIC from the current model set.
    ## If the best-fitting model from the current model set includes more than one predictor.
    if(sum(best_fit_model_from_current_set_predictors) > 1){
      ## Get the names of predictors of the best-fitting model from the current model set.
      best_fit_model_from_current_set_predictors<-
        colnames(best_fit_model_from_current_set_predictors)[best_fit_model_from_current_set_predictors==1]
      ## Get new predictor combinations to fit models to.
      pred_comb<-combn(x=best_fit_model_from_current_set_predictors,
                       m=length(best_fit_model_from_current_set_predictors)-1,
                       simplify=F)
    } else {
      ## If the best-fitting model from the current model set has a single predictor,
      ## then set the new predictor combination to NA. This will result in an
      ## intercept-only model in the next set.
      pred_comb<-vector(mode="list",length=1)
      pred_comb[[1]]<-NA
    }
    
    # Increment the model set by one.
    model_set<-model_set+1
    
    # Set the current set's best-fitting model as the next model to compare to.
    best_fit_model_from_last_set<-best_fit_model_from_current_set
    
  }
  
}

# Run the best-fit model and save its HMC output and convergence summaries.

# Get best-fit model from the model selection data frame.
best_mod<-df[df$Model=="Best",]

# Get predictors from the best-fit model.
best_mod_pred<-best_mod[,3:(3+ncol(pred_matrix_scaled))]

# If the best-fitting model includes any predictors.
if(sum(best_mod_pred) > 0){
  # Get the names of predictors of the best-fitting model.
  new_pred<-colnames(best_mod_pred)[best_mod_pred==1]
} else {
  # If the best-fitting model has no predictors, then set the new predictor combination to NA.
  # This will result in an intercept-only model.
  new_pred<-NA
}

# Check if stratum is included in the predictor combination.
stratum_included<-ifelse("Stratum" %in% new_pred,"Yes","No")

# Check if any non-stratum predictors are included in the predictor combination.
nonstratum_included<-ifelse(any(new_pred!="Stratum",na.rm=T),"Yes","No")

# If non-stratum predictors are included, get new non-stratum predictor matrix
# with the predictor combination.
if(nonstratum_included=="Yes"){
  pred_matrix_scaled_reduced<-pred_matrix_scaled[,colnames(pred_matrix_scaled) %in% new_pred]
}

# Report beginning sampling of best-fit model.
print("- Beginning sampling of best-fit model.")

# Store results of initial proc.time for timing purposes.
ptm<-proc.time()

# If included as predictors, use model with both stratum and non-stratum predictors.
if(stratum_included=="Yes" & nonstratum_included=="Yes"){
  
  # Define model parameters.
  data_Model<-list(
    "NSamples"=nrow(reads), # Number of samples.
    "NTaxa"=ncol(reads), # Number of taxa.
    "NStrata"=ncol(stratum_matrix_scaled), # Number of strata.
    "StratumMatrix"=stratum_matrix_scaled, # Predictor matrix for strata.
    "NPredictors"=ncol(pred_matrix_scaled_reduced), # Number of non-strata predictors.
    "PredictorMatrix"=pred_matrix_scaled_reduced, # Non-strata predictor matrix.
    "ReadsMatrix"=reads, # Reads response matrix.
    "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
  )
  
  # Compile Stan program into C++ code.
  DM_Reg<-stan_model(file=model_script_in_file_with_both,model_name="DM_Reg")
  
  # Run the model with HMC sampling.
  mod_fit<-sampling(DM_Reg,
                    data=data_Model,
                    chains=number_of_chains,
                    warmup=number_warmup_iterations,
                    iter=number_warmup_iterations+number_sampling_iterations,
                    thin=thinning_interval,
                    seed=stan_seed,
                    algorithm="NUTS",
                    cores=number_of_cores,
                    pars<-c("beta_0","beta_stratum","beta_pred",
                            "sigma_stratum","theta"),
                    verbose=F,
                    control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                    refresh=0
  )
  
} else if(stratum_included=="No" & nonstratum_included=="Yes"){
  
  # If stratum is not included but non-stratum predictors are,
  # use model without stratum and with non-stratum predictors.
  
  # Define model parameters.
  data_Model<-list(
    "NSamples"=nrow(reads), # Number of samples.
    "NTaxa"=ncol(reads), # Number of taxa.
    "NPredictors"=ncol(pred_matrix_scaled_reduced), # Number of non-strata predictors.
    "PredictorMatrix"=pred_matrix_scaled_reduced, # Non-strata predictor matrix.
    "ReadsMatrix"=reads, # Reads response matrix.
    "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
  )
  
  # Compile Stan program into C++ code.
  DM_Reg<-stan_model(file=model_script_in_file_without_stratum,model_name="DM_Reg")
  
  # Run the model with HMC sampling.
  mod_fit<-sampling(DM_Reg,
                    data=data_Model,
                    chains=number_of_chains,
                    warmup=number_warmup_iterations,
                    iter=number_warmup_iterations+number_sampling_iterations,
                    thin=thinning_interval,
                    seed=stan_seed,
                    algorithm="NUTS",
                    cores=number_of_cores,
                    pars<-c("beta_0","beta_pred","theta"),
                    verbose=F,
                    control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                    refresh=0
  )
  
} else if(stratum_included=="Yes" & nonstratum_included=="No"){
  
  # If stratum is included but non-stratum predictors are not,
  # use model with stratum and without non-stratum predictors.
  
  # Define model parameters.
  data_Model<-list(
    "NSamples"=nrow(reads), # Number of samples.
    "NTaxa"=ncol(reads), # Number of taxa.
    "NStrata"=ncol(stratum_matrix_scaled), # Number of strata.
    "StratumMatrix"=stratum_matrix_scaled, # Predictor matrix for strata.
    "ReadsMatrix"=reads, # Reads response matrix.
    "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
  )
  
  # Compile Stan program into C++ code.
  DM_Reg<-stan_model(file=model_script_in_file_without_nonstratum,model_name="DM_Reg")
  
  # Run the model with HMC sampling.
  mod_fit<-sampling(DM_Reg,
                    data=data_Model,
                    chains=number_of_chains,
                    warmup=number_warmup_iterations,
                    iter=number_warmup_iterations+number_sampling_iterations,
                    thin=thinning_interval,
                    seed=stan_seed,
                    algorithm="NUTS",
                    cores=number_of_cores,
                    pars<-c("beta_0","beta_stratum","sigma_stratum","theta"),
                    verbose=F,
                    control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                    refresh=0
  )
  
} else {
  
  # If neither stratum nor non-stratum predictors are included,
  # use the intercept only model.
  
  # Define model parameters.
  data_Model<-list(
    "NSamples"=nrow(reads), # Number of samples.
    "NTaxa"=ncol(reads), # Number of taxa.
    "ReadsMatrix"=reads, # Reads response matrix.
    "sd_prior"=1 # Regression coefficient and precision standard deviation prior.
  )
  
  # Compile Stan program into C++ code.
  DM_Reg<-stan_model(file=model_script_in_file_intercept_only,model_name="DM_Reg")
  
  # Run the model with HMC sampling.
  mod_fit<-sampling(DM_Reg,
                    data=data_Model,
                    chains=number_of_chains,
                    warmup=number_warmup_iterations,
                    iter=number_warmup_iterations+number_sampling_iterations,
                    thin=thinning_interval,
                    seed=stan_seed,
                    algorithm="NUTS",
                    cores=number_of_cores,
                    pars<-c("beta_0","theta"),
                    verbose=F,
                    control=list(adapt_delta=adapt_delta_argument,max_treedepth=max_treedepth_argument),
                    refresh=0
  )
  
}

# Report finished with sampling the best-fit model.
print(paste0("- Finished sampling best-fit model: Run time of ",elapsed.time(ptm)))

# Save HMC samples of the best-fit model.
print("- Saving HMC samples of best-fit model.")
save(mod_fit,file="mod_fit.RData")

# Get convergence statistics.
convergence<-as.data.frame(summary(mod_fit)$summary[,c("n_eff","Rhat")])

# Write out convergence statistics of the best-fit model.
print("- Writing out convergence statistics of best-fit model.")
write.csv(convergence,file="convergence.csv",row.names=T)

# Create a simplified model selection summary data frame.
## Copy model selection summary data frame.
df2<-df
## Remove the best-fit model record from the data frame.
df2<-df2[df2$Model!="Best",]
## Create an empty storage data frame.
df3<-data.frame(NULL)
## Loop through each model set.
for(i in unique(df2$Set)){
  ## Get the model set.
  df2_sub<-df2[df2$Set==i,]
  ## Get the best-fitting model from the set.
  best_fit_model_from_set<-df2_sub[which(df2_sub$WAIC==min(df2_sub$WAIC))[1],]
  ## Add the best-fitting model from the set to the storage data frame.
  df3<-rbind(df3,best_fit_model_from_set)
}
## Rename the set column to a field which indicates the number of predictors in the model.
colnames(df3)[1]<-"Num_Predictors"
## Clear the contents of the number of predictors field.
df3$Num_Predictors<-NA
## Rename the model column to a field which indicates if the model is the best-fit model.
colnames(df3)[2]<-"Best"
## Clear the contents of the best field.
df3$Best<-NA
## Get the predictor used in the models.
set_best_mod_pred<-df3[,3:(3+ncol(pred_matrix_scaled))]
## Loop through the best-fit model from each set.
for(i in 1:nrow(set_best_mod_pred)){
  ## Get the number of predictors included in the model.
  df3$Num_Predictors[i]<-sum(set_best_mod_pred[i,])
  ## Check to see if the best-fit model from the set is the overall best-fit model.
  df3$Best[i]<-ifelse(identical(unlist(set_best_mod_pred[i,]),unlist(best_mod_pred)),"Yes","No")
}

# Write out the simplified model selection summary.
print("- Writing out simplified model selection summary.")
write.csv(df3,file="model_selection_simple.csv",row.names=F)

# Done!
print("Done!")
