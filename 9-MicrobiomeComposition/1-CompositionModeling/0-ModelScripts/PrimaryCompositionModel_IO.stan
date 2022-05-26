//// Dirichlet-Multinomial Regression Model with Intercept Only.

data{

  //// Define variables.
  int<lower=1> NSamples; // Number of samples.
  int<lower=1> NTaxa; // Number of taxa.
  int ReadsMatrix[NSamples,NTaxa]; // Reads response matrix.
  real<lower=0> sd_prior; // Regression coefficient and precision parameter standard deviation prior.

}

parameters{

  //// Specify parameters.
  simplex[NTaxa] p[NSamples]; // Each sample's set of modeled proportions.
  vector[NTaxa-1] beta_0_raw; // Intercept raw regression coefficient vector.
  real theta; // Precision parameter.

}

transformed parameters{

  //// Specify transformed parameters.
  vector[NTaxa] beta_0; // Intercept regression coefficient vector.
  real exptheta; // Exponentiated precision parameter.

  //// Define transformed parameters.

  // Intercept regression coefficients.
  beta_0[NTaxa]=0; // Set last taxon's intercept to zero.
  for(k in 1:(NTaxa-1)){
    beta_0[k]=beta_0_raw[k]; // Set intercepts of other taxa to the raw betas.
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
      eta[j]=beta_0[j];
    }

    // Relate sample eta's to probabilities with the Dirichlet distribution.
    p[i]~dirichlet(softmax(eta[1:NTaxa])*exptheta);

    // Relate probabilities to counts with the multinomial distribution.
    ReadsMatrix[i,1:NTaxa]~multinomial(p[i]);

  }

}
