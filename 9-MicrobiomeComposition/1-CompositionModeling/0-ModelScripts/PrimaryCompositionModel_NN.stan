//// Dirichlet-Multinomial Regression Model without Non-Stratum Predictors.

data{

  //// Define variables.
  int<lower=1> NSamples; // Number of samples.
  int<lower=1> NTaxa; // Number of taxa.
  int<lower=1> NStrata; // Number of strata.
  matrix[NSamples,NStrata] StratumMatrix; // Predictor matrix for strata.
  int ReadsMatrix[NSamples,NTaxa]; // Reads response matrix.
  real<lower=0> sd_prior; // Regression coefficient and precision parameter standard deviation prior.

}

parameters{

  //// Specify parameters.
  simplex[NTaxa] p[NSamples]; // Each sample's set of modeled proportions.
  vector[NTaxa-1] beta_0_raw; // Intercept raw regression coefficient vector.
  matrix[NTaxa-1,NStrata-1] beta_stratum_raw; // Stratum raw regression coefficient matrix.
  vector<lower=0>[NTaxa-1] sigma2_stratum; // Stratum common variance vector.
  real theta; // Precision parameter.

}

transformed parameters{

  //// Specify transformed parameters.
  vector[NTaxa] beta_0; // Intercept regression coefficient vector.
  matrix[NTaxa,NStrata] beta_stratum; // Stratum regression coefficient matrix.
  vector<lower=0>[NTaxa-1] sigma_stratum; // Stratum common standard deviation vector.
  real exptheta; // Exponentiated precision parameter.

  //// Define transformed parameters.

  // Intercept regression coefficients.
  beta_0[NTaxa]=0; // Set last taxon's intercept to zero.
  for(k in 1:(NTaxa-1)){
    beta_0[k]=beta_0_raw[k]; // Set intercepts of other taxa to the raw betas.
  }

  // Stratum regression coefficients.
  for(j in 1:NStrata){
    beta_stratum[NTaxa,j]=0; // Set last taxon's stratum regression coefficients to zero.
  }
  for(k in 1:(NTaxa-1)){
    for(j in 1:(NStrata-1)){
      beta_stratum[k,j]=beta_stratum_raw[k,j]; // Set stratum regression coefficients of other taxa to the raw betas (execpt for the last stratum).
    }
    beta_stratum[k,NStrata]=-sum(beta_stratum_raw[k,1:(NStrata-1)]); // Apply a sum to zero constraint for the last stratum.
    sigma_stratum[k]=sqrt(sigma2_stratum[k]); // Calculate stratum common standard deviation as the square root of the common variance.
  }

  // Create exponentiated precision parameter.
  exptheta=exp(theta);

}

model{

  //// Provide priors.

  // Intercept regression coefficients.
  for(k in 1:(NTaxa-1)){
    beta_0_raw[k]~normal(0,sd_prior);
  }

  // Stratum regression coefficients.
  for(k in 1:(NTaxa-1)){
    for(j in 1:(NStrata-1)){
      beta_stratum_raw[k,j]~normal(0,sigma_stratum[k]);
    }
  }

  // Stratum common variance.
  for(k in 1:(NTaxa-1)){
    sigma2_stratum[k]~inv_gamma(0.01,0.01);
  }

  // Precision parameter.
  theta~normal(0,sd_prior);

  //// State likelihood.

  // Loop through each sample.
  for(i in 1:NSamples){

    // Define local variable.
    vector[NTaxa] eta;

    // Loop through each taxon.
    for(j in 1:NTaxa){
      // Linear model for eta.
      eta[j]=beta_0[j]+beta_stratum[j,1:NStrata]*transpose(StratumMatrix[i,1:NStrata]);
    }

    // Relate sample eta's to probabilities with the Dirichlet distribution.
    p[i]~dirichlet(softmax(eta[1:NTaxa])*exptheta);

    // Relate probabilities to counts with the multinomial distribution.
    ReadsMatrix[i,1:NTaxa]~multinomial(p[i]);

  }

}
