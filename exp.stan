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
}

parameters {
  vector[p] beta;
}

transformed parameters{
  vector[N_obs] mu_obs=exp(X_obs*beta);
  vector[N_cens] mu_cens=exp(X_cens*beta);
}

model {
  beta~multi_normal(beta0, Sigma0);
  
  t_obs~ exponential(mu_obs); 
  target += exponential_lccdf(t_cens|mu_cens); 
}

generated quantities {
  vector[N_obs] log_lik_obs;
  vector[N_cens] log_lik_cens;
  
  for (i in 1:N_obs) {
    log_lik_obs[i] = exponential_lpdf(t_obs[i] | mu_obs[i]);
  }
  
  for (i in 1:N_cens) {
    log_lik_cens[i] = exponential_lccdf(t_cens[i] | mu_cens[i]);
  }
}

