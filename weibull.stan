///////////////////////// DATA /////////////////////////////////////
data {
  int<lower=0> p;              // Numero di predittori
  int<lower=0> N_obs;          // Numero di osservazioni non censurate
  int<lower=0> N_cens;         // Numero di osservazioni censurate
  vector[N_obs] t_obs;
  vector[N_cens] t_cens;

  matrix[N_obs, p] X_obs;
  matrix[N_cens, p] X_cens;

  vector[p] beta0;
  matrix[p, p] Sigma0;
  
  real<lower=0> a0;     // Prior su alpha (forma della Weibull)
  real<lower=0> b0;     // Prior su alpha (forma della Weibull)
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
  vector[p] beta;          // Coefficienti fissi
  real<lower=0> alpha;     // Parametro di forma Weibull
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters {
  vector[N_obs] lambda_obs = exp(X_obs * beta);  
  vector[N_cens] lambda_cens = exp(X_cens * beta);
}

////////////////// MODEL ////////////////////////
model {
  // Priori
  beta ~ multi_normal(beta0, Sigma0);
  alpha ~gamma(a0,b0);

  // Likelihood per i dati osservati
  t_obs ~ weibull(alpha,lambda_obs);

  // Likelihood per i dati censurati
  target += weibull_lccdf(t_cens | alpha, lambda_cens);
}

generated quantities {
  vector[N_obs] log_lik_obs;
  vector[N_cens] log_lik_cens;
  
  for (i in 1:N_obs) {
    log_lik_obs[i] = weibull_lpdf(t_obs[i] | alpha, lambda_obs[i]);
  }
  
  for (i in 1:N_cens) {
    log_lik_cens[i] = weibull_lccdf(t_cens[i] |alpha,  lambda_cens[i]);
  }
}
