data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> N_obs;
  int<lower=0> N_cens;
  vector[N_obs] t_obs;
  vector[N_cens] t_cens;

  matrix[N_obs, p] X_obs;
  matrix[N_cens, p] X_cens;
  
  vector[p] beta0;
  matrix[p,p] Sigma0;
  real alpha0;
  real lambda0;
}

parameters {
  vector[p] beta;
  real<lower=0> tau;
}

transformed parameters{
  vector[N_obs] mu_obs=X_obs*beta;
  vector[N_cens] mu_cens=X_cens*beta;
}

model {
  tau~gamma(alpha0/2 ,lambda0/2);
  beta~multi_normal(beta0, Sigma0/tau);
  
  t_obs ~ lognormal(mu_obs, pow(1/tau,0.5)); 
  target += lognormal_lccdf(t_cens|mu_cens,pow(1/tau,0.5)); 
}

generated quantities {
  vector[N_obs] log_lik_obs;
  vector[N_cens] log_lik_cens;
  
  for (i in 1:N_obs) {
    log_lik_obs[i] = lognormal_lpdf(t_obs[i] | mu_obs[i], pow(1/tau, 0.5));
  }
  
  for (i in 1:N_cens) {
    log_lik_cens[i] = lognormal_lccdf(t_cens[i] | mu_cens[i], pow(1/tau, 0.5));
  }
}

