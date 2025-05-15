library(loo)
to_one_hot <- function(x) {
  xf<-as.factor(x)
  xn<-as.numeric(xf)
  results <- matrix(0, nrow = length(x), ncol = max(xn))
  colnames(results)<-levels(xf)
  
  for (i in 1:length(xn))
    results[i, xn[i]] <- 1
  results
}
my_waic <- function(out_pred_dens){
  LPPD <- sum(log(colMeans(out_pred_dens)))
  p_waic <- mean(apply(out_pred_dens, 2, function(x) var(log(x))))
  return(- 2 * LPPD + 2 * p_waic)
}

my_lpml <- function(out_pred_dens){
  out <- sum(log(1 / colMeans(1 / out_pred_dens)))
  return(out)
}
library(rstan)
library(bayestestR)
library(bayesplot)
library(rstanarm)
library(ggplot2)
library(coda)
library(BayesSurvival)
library(survival)


dati<-read.csv("gbsg.csv")
str(dati)
hist(dati$pgr)
hist(dati$nodes)
hist(dati$er)
hist(dati$pgr)
#tolgo id e x, pgr, er 
dati<-dati[, -c(1,2)]
str(dati)

#tolgo i dati estremi
ind<-which(dati$pgr>1000)
dati<-dati[-ind,]

#tempi di soppravvivenza
summary(dati$rfstime)
table(dati$status) 

dati_s<- data.frame(scale(dati[,c(1,3,5,6,7) ]), meno=dati$meno, grade=dati$grade, hormon=dati$hormon, status=dati$status, rfstime=dati$rfstime)

#stima frequentista 
library(survival)
library(flexsurv)
model<- survreg(Surv(rfstime, status)~+0+age+meno+grade+nodes+hormon+size+pgr+er, data=dati_s, dist = "exponential" )
exp<-step(model, direction = "backward")
AIC(exp)
he<-flexsurvreg(Surv(rfstime, status)~1, data=dati)
plot(he, type="hazard")
modelw<- survreg(Surv(rfstime, status)~0+age+meno+grade+nodes+hormon+size+pgr+er, data=dati_s, dist = "weibull" )
step(modelw, direction = "backward")

modelog<- survreg(Surv(rfstime, status)~0+age+meno+grade+nodes+hormon+size+pgr+er, data=dati_s, dist = "lognormal" )
step(modelog, direction = "backward")

# Kaplan-Meier estimate
fit_KM <- survfit(Surv(rfstime, status) ~ 1, data=dati)
ggsurvplot(fit_KM,
           data = dati,
           conf.int = TRUE,         # mostra intervalli di confidenza
           pval = TRUE,             # mostra p-value (se c'è confronto tra gruppi)
           risk.table = TRUE,       # tabella con n° a rischio nel tempo
           surv.median.line = "hv", # linea mediana orizzontale + verticale
           ggtheme = theme_minimal(),
           palette = "magenta")


#Data preprocessing

y<- dati$rfstime
y_cens<-dati$rfstime[dati$status==0]
y_obs<-dati$rfstime[dati$status==1]
X<- dati[, -9] #tolgo rfstime
X_s<- data.frame(scale(X[,c(1,3,5,6, 7) ]), meno=X$meno, grade=X$grade, hormon=X$hormon, status=X$status)
X_obs<- X_s[ X_s$status==1, ] 
summary(X_obs)
X_obs<-X_obs[, -9] #tolgo status
X_cens<- X_s[X_s$status==0, ]
summary(X_cens)
X_cens<-X_cens[, -9]

#Exponential model
datablock<-list( n=length(y), 
                 p=ncol(X_obs),
                 N_obs=nrow(X_obs),
                 N_cens=nrow(X_cens),
                 t_obs=y_obs,
                 t_cens=y_cens,
                 X_obs=X_obs,
                 X_cens=X_cens, 
                 beta0=rep(0,ncol(X_obs)),
                 Sigma0=diag(10^4, ncol(X_obs))
)

model_exp<-stan(file="exp.stan", 
                data=datablock, 
                chains=2, 
                iter=5000, 
                warmup=1000,
                seed=123)

rstan::traceplot(model_exp, pars = c("beta"))
rstan::traceplot(model_exp, pars= "beta", inc_warmup = FALSE)
params_STAN_exp <- As.mcmc.list(model_exp, pars = c("beta"))
summary(params_STAN_exp)


effectiveSize(params_STAN_exp)
ci(params_STAN_exp, method = "ETI")
ci(params_STAN_exp, method = "HDI")

#model evaluation
loglikcens<- extract(model_exp, pars = c("log_lik_cens"))[[1]]
loglikoss<-extract(model_exp, pars = c("log_lik_obs"))[[1]]
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
loo(out_pred_dens_STAN)
waic(out_pred_dens_STAN) 

posterior<- as.matrix(params_STAN_exp)
mcmc_areas(posterior)


#Weibull

y<- dati$rfstime
y_cens<-dati$rfstime[dati$status==0]
y_obs<-dati$rfstime[dati$status==1]
X<- dati[, -9] 
X_s<- data.frame(scale(X[,c(1,3,5,6) ]), meno=X$meno, grade=X$grade, hormon=X$hormon, status=X$status) #Without rfstime
X_obs<- X_s[ X_s$status==1, ] 
summary(X_obs)
X_obs<-X_obs[, -8] #without status
X_cens<- X_s[X_s$status==0, ]
summary(X_cens)
X_cens<-X_cens[, -8]


datablock<-list( n=length(y), 
                 p=ncol(X_obs),
                 N_obs=nrow(X_obs),
                 N_cens=nrow(X_cens),
                 t_obs=y_obs,
                 t_cens=y_cens,
                 X_obs=X_obs,
                 X_cens=X_cens, 
                 beta0=rep(0,ncol(X_obs)),
                 Sigma0=diag(10^4, ncol(X_obs)),
                 a0=20,
                 b0=2
)

model_weib<-stan(file="weibull.stan", 
                 data=datablock, 
                 chains=2, 
                 iter=5000, 
                 warmup=1000,
                 seed=123)
library(coda)
rstan::traceplot(model_weib, pars = c("beta", "alpha"))
params_STAN_weib <- As.mcmc.list(model_weib, pars = c("beta", "alpha"))

summary(params_STAN_weib)


effectiveSize(params_STAN_weib)
ci(params_STAN_weib, method = "ETI")
ci(params_STAN_weib, method = "HDI")

#model evaluation
loglikcens<- extract(model_weib, pars = c("log_lik_cens"))[[1]]
loglikoss<-extract(model_weib, pars = c("log_lik_obs"))[[1]]
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
loo(out_pred_dens_STAN) 
waic(out_pred_dens_STAN)


loglikcens<- exp(extract(model_weib, pars = c("log_lik_cens"))[[1]])
loglikoss<-exp(extract(model_weib, pars = c("log_lik_obs"))[[1]])
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
my_lpml(out_pred_dens_STAN)
my_waic(out_pred_dens_STAN) 

###################
#Lognormal

datablock_ln<-list( n=length(y), 
                    p=ncol(X_obs),
                    N_obs=nrow(X_obs),
                    N_cens=nrow(X_cens),
                    t_obs=y_obs,
                    t_cens=y_cens,
                    X_cens=X_cens, 
                    X_obs=X_obs,
                    beta0=rep(0,ncol(X_obs)),
                    Sigma0=diag(100, ncol(X_obs)),
                    alpha0=2,  
                    lambda0=2) 

model_lognormal<-stan(file="lognormal.stan", 
                      data=datablock_ln, 
                      chains=2, 
                      iter=5000, 
                      warmup=1000,
                      seed=123)
rstan::traceplot(model_lognormal, pars = c("beta", "tau"))
params_STAN_log<- As.mcmc.list(model_lognormal, pars = c("beta", "tau"))
summary(params_STAN_log)

effectiveSize(params_STAN_log)
ci(params_STAN_log, method = "ETI")
ci(params_STAN_log, method = "HDI")

#model evaluation
loglikcens<- extract(model_lognormal, pars = c("log_lik_cens"))[[1]]
loglikoss<-extract(model_lognormal, pars = c("log_lik_obs"))[[1]]
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
loo(out_pred_dens_STAN) 
waic(out_pred_dens_STAN)

loglikcens<- exp(extract(model_lognormal, pars = c("log_lik_cens"))[[1]])
loglikoss<-exp(extract(model_lognormal, pars = c("log_lik_obs"))[[1]])
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
my_lpml(out_pred_dens_STAN) 
my_waic(out_pred_dens_STAN) 


#Plots
trace<-rstan::traceplot(model_lognormal, pars = c("beta", "tau"))
trace + facet_wrap(~ parameter, labeller = labeller(
  parameter = c(beta = "Coefficiente β", tau = "Precisione τ")
))


posterior <- as.array(model_lognormal)

original_param_names <- dimnames(posterior)$parameters

new_param_names <- original_param_names
beta_labels <- c("age", "size", "nodes", "pgr", "meno", "grade", "hormon")
for (i in seq_along(beta_labels)) {
  new_param_names <- gsub(paste0("beta\\[", i, "\\]"), beta_labels[i], new_param_names)
}

dimnames(posterior)$parameters <- new_param_names

mcmc_trace(posterior, pars = c(beta_labels, "tau"))

posterior <- as.matrix(params_STAN_log)

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars=c("beta[1]","beta[2]", "beta[3]", "beta[4]"),
           prob = 0.8) + plot_title

mcmc_areas(posterior,
           pars=c("beta[5]","beta[6]", "beta[7]"),
           prob = 0.8) + plot_title


#Variable selection with spike and slab method on Lognormal Model
c <- 100
k <- 0.05
tau2 <- k / sqrt(2 * log(c) * c^2 / (c^2 - 1))


datablock_ln_vs<-list( n=length(y), 
                       p=ncol(X_obs),
                       N_obs=nrow(X_obs),
                       N_cens=nrow(X_cens),
                       t_obs=y_obs,
                       t_cens=y_cens,
                       X_cens=X_cens, 
                       X_obs=X_obs,
                       beta0=rep(0,ncol(X_obs)),
                       alpha0=2,  #2
                       lambda0=2, 
                       alpha1 = 1, 
                       alpha2 = 1,
                       tau2 = tau2^2, 
                       c = c^2)

model_lognormal_vs<-stan(file="log_normal_varsel.stan", 
                         data=datablock_ln_vs, 
                         chains=1, 
                         iter=5000, 
                         warmup=1000,
                         seed=123)
rstan::traceplot(model_lognormal, pars = c("beta"))
params <- As.mcmc.list(extract(model_lognormal, pars = c("beta"))[[1]])
summary(params)
CI(params)
geweke.diag(params)

# -----------------------------
# -----------------------------

gamma_matrix <- ifelse(abs(params) > k, 1, 0)

# HPD 

unique_model <- unique(gamma_matrix, MARGIN = 1)
dim(unique_model)
freq <- apply(unique_model, 1, function(b)
  sum(apply(gamma_matrix, MARGIN = 1, function(a) mean(a == b) == 1)))
HPD_model <- unique_model[which.max(freq),]
HPD_model
# 1 1 1 1 1 1 1 1
# MPM
MPM_model <- as.numeric(colMeans(gamma_matrix) > 0.5)
#1 1 1 1 1 1 1 1

# HS
HS_model <- as.numeric(colMeans(gamma_matrix) == 1)
HS_model
#1 0 1 1 0 1 1 1


#Survival plots for lognormal 
beta<- as.matrix(model_lognormal, pars = c("beta"))
tau<-as.matrix(model_lognormal, pars = c("tau"))
grid <- seq(0, 2700, by=20)
f_t<-function(t, X, beta, tau){
  f<-c()
  for ( i in 1:nrow(beta) ){
    mu<-X%*%beta[i,]
    sigma<-1/sqrt(tau[i])
    f[i]<-dlnorm( t, mu, sigma)
  }
  out<-c(mean(f),quantile(f, 0.025), quantile(f, 0.975))
  return(out)
}
gamma[i]
s_t<- function(t, X, beta, tau){
  s<-c()
  for ( i in 1:nrow(beta) ){
    mu<-X%*%beta[i,]
    sigma<-1/sqrt(tau[i])
    s[i]<-1-pnorm((log(t)-mu)/sigma)
  }
  out<-c(mean(s),quantile(s, 0.025), quantile(s, 0.975))
  return(out)
}
cov1<-c(median(dati_s$age),median(dati_s$size), median(dati_s$nodes),median(dati_s$pgr), 0,2,0)
F_mean <- numeric(length(grid))
F_upper <- numeric(length(grid))
F_lower <- numeric(length(grid))

for (i in 1:length(grid)) {
  F_mean[i] <- f_t(grid[i],cov1, beta, tau)[1]
  F_lower[i] <- f_t(grid[i],cov1, beta, tau)[2]
  F_upper[i] <- f_t(grid[i],cov1, beta, tau)[3]
}

cov1<-c(median(dati_s$age),median(dati_s$size), median(dati_s$nodes),median(dati_s$pgr), 0,2,0)
S_mean <- numeric(length(grid))
S_upper <- numeric(length(grid))
S_lower <- numeric(length(grid))

for (i in 1:length(grid)) {
  S_mean[i] <- s_t(grid[i],cov1, beta, tau)[1]
  S_lower[i] <- s_t(grid[i],cov1, beta, tau)[2]
  S_upper[i] <- s_t(grid[i],cov1, beta, tau)[3]
}


data_plot1 <- data.frame(Time = grid, Mean = S_mean, Upper = S_upper, Lower = S_lower)
data_plotf <- data.frame(Time = grid, Mean = F_mean/S_mean, Upper = F_upper/S_upper, Lower = F_lower/S_lower)


#Curva 2 per mediane e 0 2 1 
cov2<-c(median(dati_s$age),median(dati_s$size), median(dati_s$nodes),median(dati_s$pgr), 0,2,1)
S_mean <- numeric(length(grid))
S_upper <- numeric(length(grid))
S_lower <- numeric(length(grid))

for (i in 1:length(grid)) {
  S_mean[i] <- s_t(grid[i],cov2, beta, tau)[1]
  S_lower[i] <- s_t(grid[i],cov2, beta, tau)[2]
  S_upper[i] <- s_t(grid[i],cov2, beta, tau)[3]
}

data_plot2 <- data.frame(Time = grid, Mean = S_mean, Upper = S_upper, Lower = S_lower)

library(ggplot2)

ggplot() +
  # Prima curva (data_plot)
  geom_line(data = data_plotf, aes(x = Time, y = Mean, color = "Hormon 0"), size = 1) +
  geom_ribbon(data = data_plotf, aes(x = Time, ymin = Lower, ymax = Upper, fill = "Hormon 0"), alpha = 0.2) +
  
  # Seconda curva (data_plot1)
  geom_line(data = data_plot2, aes(x = Time, y = Mean, color = "Hormon 1"), size = 1) +
  geom_ribbon(data = data_plot2, aes(x = Time, ymin = Lower, ymax = Upper, fill = "Hormon 1"), alpha = 0.2) +
  
  # Temi e assi
  theme_bw() +
  ylab("Funzione di sopravvivenza") +
  xlab("Tempo") +
  ylim(0, 1) +
  ggtitle("Curve di sopravvivenza stimate (modello Lognormal)") +
  
  # Personalizzazione colori legenda
  scale_color_manual(name = "Trattamento", values = c("Hormon 0" = "magenta", "Hormon 1" = "orange")) +
  scale_fill_manual(name = "Trattamento", values = c("Hormon 0" = "magenta", "Hormon 1" = "orange"))

#############
##Bayes Factor between Lognormal and Weibull models
logsumexp <- function(x){
  M <- max(x)
  lse <- M + log(sum(exp(x - M)))
  return(lse)
}

#estraggo log lik del weibull
out_pred_dens_w <- exp(extract(model_weib, pars = c("log_lik_cens", "log_lik_obs"))[[1]])
#estraggo log lik del log-normal
out_pred_dens_log <- exp(extract(model_lognormal, pars = c("log_lik_cens", "log_lik_obs"))[[1]])

mean(apply(out_pred_dens_w, 1, prod))
mean(apply(out_pred_dens_log, 1, prod))

M0<-log(nrow(out_pred_dens_w)) - logsumexp(-apply(out_pred_dens_w, 1, sum)) #mod weibull M0
M1<-log(nrow(out_pred_dens_log)) - logsumexp(-apply(out_pred_dens_log, 1, sum)) #mod log normal M1


bayes_factor<-M0/M1


############
#Lognormal model with random effect by "grade" 
y<- dati$rfstime
y_cens<-dati$rfstime[dati$status==0]
y_obs<-dati$rfstime[dati$status==1]
X<- dati[, -9]

X_s<- data.frame(scale(X[,c(1,3,5,6) ]), meno=X$meno, grade=X$grade, hormon=X$hormon, status=X$status)
X_obs_c<- X_s[ X_s$status==1, ] 
grade_obs<-X_obs_c$grade

summary(X_obs_c)
X_obs_c<-X_obs_c[, -c(6,8)] #without grade and status

X_cens_c<- X_s[X_s$status==0, ]
summary(X_cens_c)
grade_cens<-X_cens_c$grade
X_cens_c<-X_cens_c[, -c(6,8)] #without grade and status

ngr1<-3 #number of levels of random effect

datablock_ln<-list( n=length(y), 
                    p=ncol(X_obs_c),
                    N_obs=nrow(X_obs_c),
                    N_cens=nrow(X_cens_c),
                    t_obs=y_obs,
                    t_cens=y_cens,
                    X_cens=X_cens_c, 
                    X_obs=X_obs_c,
                    J=ngr1, 
                    group_obs=grade_obs, 
                    group_cens=grade_cens, 
                    beta0=rep(0,ncol(X_obs_c)),
                    Sigma0=diag(10^4, ncol(X_obs_c)),
                    alpha0=2,
                    lambda0=2, 
                    tau0=10^4)

model_logRE1<-stan(file="lognormal_RE_ns.stan", 
                   data=datablock_ln, 
                   chains=2, 
                   iter=5000, 
                   warmup=1000,
                   seed=123)
rstan::traceplot(model_logRE1, pars = c("beta", "u", "tau"))
params_STAN_logre<- As.mcmc.list(model_logRE1, pars = c("beta", "u", "tau"))
summary(params_STAN_logre)

ci(params_STAN_logre, method="HDI")

effectiveSize(params_STAN_logre)
geweke.diag(params_STAN_logre)

#model evaluation
loglikcens<- extract(model_logRE1, pars = c("log_lik_cens"))[[1]]
loglikoss<-extract(model_logRE1,  pars = c("log_lik_obs"))[[1]]
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
loo(out_pred_dens_STAN) 
waic(out_pred_dens_STAN)

loglikcens<- exp(extract(model_logRE1, pars = c("log_lik_cens"))[[1]])
loglikoss<-exp(extract(model_logRE1, pars = c("log_lik_obs"))[[1]])
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
my_lpml(out_pred_dens_STAN)
my_waic(out_pred_dens_STAN) 


#Traceplots
posterior <- as.array(model_logRE1)
original_param_names <- dimnames(posterior)$parameters
new_param_names <- original_param_names
beta_labels <- c("age", "size", "nodes", "pgr", "meno", "hormon")
for (i in 1:9) {
  new_param_names[i] <- gsub(paste0("beta\\[", i, "\\]"), beta_labels[i], new_param_names), 
  paste0("u\\[", i, "\\]"), u_labels[i], new_param_names)
}
u_labels <- c("grade1", "grade2", "grade3")
for (i in seq_along(u_labels)) {
  new_param_names_u<- gsub(paste0("u\\[", i, "\\]"), u_labels[i], new_param_names_u)
}
dimnames(posterior)$parameters <-  new_param_names_u

mcmc_trace(posterior, pars = c(u_labels))

posterior <- as.matrix(params_STAN_logre)[,1:10]

plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars=c("beta[1]", "beta[2]", "beta[3]", "beta[4]"),
           prob = 0.8) + plot_title

mcmc_areas(posterior,
           pars=c("beta[5]","beta[6]", "beta[7]","tau"),
           prob = 0.8) + plot_title



############
# Weibull model with random effect by "grade"  

datablock<-list( p=ncol(X_obs_c),
                 N_obs=nrow(X_obs_c),
                 N_cens=nrow(X_cens_c),
                 t_obs=y_obs,
                 t_cens=y_cens,
                 X_obs=X_obs_c,
                 X_cens=X_cens_c, 
                 J=ngr1, 
                 group_obs=grade_obs, 
                 group_cens=grade_cens, 
                 beta0=rep(0,ncol(X_obs_c)),
                 Sigma0=diag(10^4, ncol(X_obs_c)),
                 a0=2,
                 b0=2, 
                 tau0=10^4)

model_weibRE<-stan(file="weibull_RE.stan", 
                   data=datablock, 
                   chains=2, 
                   iter=5000, 
                   warmup=1000,
                   seed=123)
library(coda)
rstan::traceplot(model_weibRE, pars = c("beta", "alpha", "u"))
params_STAN_weibre <- As.mcmc.list(model_weibRE, pars = c("beta", "alpha", "u"))
summary(params_STAN_weibre)
ci(params_STAN_weibre, method="HDI")
#Model evaluation
loglikcens<- extract(model_weibRE, pars = c("log_lik_cens"))[[1]]
loglikoss<-extract(model_weibRE,  pars = c("log_lik_obs"))[[1]]
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
loo(out_pred_dens_STAN) 

waic(out_pred_dens_STAN)

loglikcens<- exp(extract(model_weibRE, pars = c("log_lik_cens"))[[1]])
loglikoss<-exp(extract(model_weibRE, pars = c("log_lik_obs"))[[1]])
out_pred_dens_STAN <-cbind(loglikcens, loglikoss)
my_lpml(out_pred_dens_STAN) 
my_waic(out_pred_dens_STAN) 

##########
#Piecewise exponential model

res <- BayesSurv(df = dati, #our data frame
                 time = "rfstime", #name of column with survival/censoring times
                 event = "status", #name of column with status indicator
                 prior = "Dependent",
                 #use dependent Gamma prior
)

summary(res)
names(res)

gg <- PlotBayesSurv(bayes.surv.object = res,
                    object = "survival", color = "magenta")

km <- survfit( Surv(rfstime, status) ~ 1, data = dati ) #Kaplan-Meier

df.km <- data.frame(t = km$time, km = km$surv)

gg <- gg + geom_line(data = df.km, aes(x = t, y = km), colour = "black", size = 1, lty = 6)
gg <- gg + labs(title = "With Kaplan-Meier + CI's")
gg

