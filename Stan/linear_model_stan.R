
# simple simulation + modeling, linear regression in Stan
# demonstrates: simple model fitting and troubleshooting
library(mvtnorm)
library(future)
library(future.apply)
library(rstan)

######### Data generating mechanism
# this is copied from CodeExamples/R/simulation_linear_regression.R
dgm_linmod <- function(
    # this function outputs a data frame with variables (y,x,z) that represent
  #  a continuous outcome, exposure, and potential confounder
  n=100, # sample size (default = 100)
  beta=c(0,1,1), # true beta coefficients (default = c(0,1,1)) - this will determine
  xzcorr=.9, # true correlation between x and z
  sigma=1.0 # standard deviation of normally distributed error term (default = 1)
){
  xz = rmvnorm(n, c(0,0), matrix(c(1, xzcorr, xzcorr, 1),nrow=2))
  x=xz[,1]
  z=xz[,2]
  y = beta[1] + xz %*% beta[2:3] + rnorm(n, 0, sigma)
  # now we can introduce things like missing data or measurement error
  # e.g.
  x_measured = x + rnorm(n, 0, .2) # classical exposure measurement error
  dat = data.frame(x=x,z=z, y=y, x_measured=x_measured)
  attr(dat, "truth") = beta # attributes can be used to keep track of how the data were generated
  return(dat)
}

# create some data
data = dgm_linmod(n=100, c(5,1,1))

######### Data analysis methods

# stnandard linear model
mod_adj = lm(y ~ x + z, data = data)
  

############################################################
# Bayesian linear model using Stan
############################################################
  
# rstan functions use a list as input for data and constants
standata = as.list(data) # data
standata$N = nrow(data) # constants
standata$X = cbind(rep(1,nrow(data)), data$x, data$z) # alternative way to express independent variables as a matrix
standata$J = ncol(standata$X)
    
##############
#  stan code for linear model
##############

stanmod <- "
data {
  int<lower=0> N;          // sample size
  int<lower=0> J;          // design matrix columns
  real y[N];               // estimated treatment effects
  matrix [N,J]X;               // estimated treatment effects
}
parameters {
  real<lower=0> sigma;
  vector[J] beta;
}
model {
  sigma ~ cauchy(0, 10); // half cauchy prior on scale
  beta[1] ~ normal(0,1000); // vague prior on intercept
  for(j in 2:J){
    beta[j] ~ normal(0,1); // shrinkage prior
  }
  for(i in 1:N){
    real mu =  X[i,] * beta;
    target += normal_lpdf(y[i] | mu, sigma);  // alternative to <<y[i] ~ normal(mu, sigma);>>
  }
}
generated quantities{

}
"

# compile model
compiledfit <- stan(model_code=stanmod, data = standata, chains=1, iter=10)

# run more samples, check for convergence
postsamples = stan(fit=compiledfit, iter=1000, data = standata)

print(postsamples)   # Bayesian fit
summary(mod_adj)  # MLE


##############
# alternative stan code for linear model
# includes extra trouble-shooting steps
##############

stanmod2 <- "
data {
  int<lower=0> N;          // sample size
  int<lower=0> J;          // design matrix columns
  real y[N];               // estimated treatment effects
  matrix [N,J]X;               // estimated treatment effects
}
parameters {
  real<lower=0> sigma;
  vector[J] beta;
}
model {
  sigma ~ cauchy(0, 10); // half cauchy prior on scale
  beta[1] ~ normal(0,1000); // vague prior on intercept
  beta[2:J] ~ normal(0,1); // shrinkage prior
  vector[N] mu = X * beta;
  print(\"\\n beta[x] =\",beta[1]);
  print(\"log density before y =\", target());
  y ~ normal(mu, sigma); // 'vectorized' form
  print(\"log density after y =\", target());
}
generated quantities{
  real ypred = [1, 0.1, 0.9] * beta; // predicting y from a test exposure, x = 0.1, z = 0.9
}
"

# compile model, get debug messages
compiledfit2 <- stan(model_code=stanmod2, data = standata, chains=1, iter=10)

# run more samples, check for convergence
postsamples2 = stan(fit=compiledfit2, iter=1000, data = standata)

print(postsamples2)   # Bayesian fit
summary(mod_adj)  # MLE
