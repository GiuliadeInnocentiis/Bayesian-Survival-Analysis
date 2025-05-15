///////////////////////// DATA /////////////////////////////////////
data {
  int<lower=0> p;              // Numero di predittori
  int<lower=0> N_obs;          // Numero di osservazioni non censurate
  int<lower=0> N_cens;         // Numero di osservazioni censurate
  int<lower=1> J;              // Numero di gruppi
  int<lower=1, upper=J> group_obs[N_obs];  // Indicatori di gruppo per dati osservati
  int<lower=1, upper=J> group_cens[N_cens];// Indicatori di gruppo per dati censurati
  vector[N_obs] t_obs;
  vector[N_cens] t_cens;

  matrix[N_obs, p] X_obs;
  matrix[N_cens, p] X_cens;

  vector[p] beta0;
  matrix[p, p] Sigma0;
  
  real<lower=0> tau0;   // Iperparametro per deviazione standard degli effetti casuali
  real<lower=0> a0;     // Prior su alpha (forma della Weibull)
  real<lower=0> b0;     // Prior su alpha (forma della Weibull)
}

//////////////////// PARAMETERS /////////////////////////////////
parameters {
  vector[p] beta;          // Coefficienti fissi
  vector[J] u;             // Effetti casuali del gruppo// Deviazione standard degli effetti casuali
  real<lower=0> alpha;     // Parametro di forma Weibull
}

//////////////////// TRANSFORMED PARAMETERS /////////////////////////////////
transformed parameters {
  vector[N_obs] lambda_obs = exp((X_obs * beta + u[group_obs]));  
  vector[N_cens] lambda_cens = exp((X_cens * beta + u[group_cens]));
}

////////////////// MODEL ////////////////////////
model {
  // Priori
  beta ~ multi_normal(beta0, Sigma0);
   // Prior su deviazione standard degli effetti casuali
  u ~ normal(0, tau0);      // Effetti casuali gerarchici
  alpha ~ gamma(a0, b0);   // Prior su parametro di forma Weibull

  // Likelihood per i dati osservati
  t_obs~ weibull(alpha, lambda_obs);

  // Likelihood per i dati censurati
  target += weibull_lccdf(t_cens | alpha, lambda_cens);
}



generated quantities {
  vector[N_obs] log_lik_obs;
  vector[N_cens] log_lik_cens;
  
  for (i in 1:N_obs) {
    log_lik_obs[i] = weibull_lpdf(t_obs[i]|alpha, lambda_obs[i]);
  }
  
  for (j in 1:N_cens) {
    log_lik_cens[j] = weibull_lccdf(t_cens[j]|alpha, lambda_cens[j]);
  }
}
