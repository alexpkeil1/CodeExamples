######################################################################################################################
# Author: Alex Keil and Juwel Rana
# Program: 
# Language: R
# Date: 
# Project: Coding examples
# Tasks: 
# Data in: 
# Description: 
# Keywords: 
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################
library(mvtnorm)
library(clusterGeneration) # just for random covariance matrix creation


# step 0: helper functions

# inverse_logit: take log odds and convert to probability
inverse_logit <- function(mu){
  1/(1+exp(-mu))  
}


# step 1: generate data
data_generation_simple <- function(
    n # sample size
    ){
  # todo: delete this later
  b0 = -1.2
  beta = c(1,-1,0)
  ### end: will delete this later
  x <- mvtnorm::rmvnorm(n, mean=c(0,0,0), sigma = diag(rep(1, 3))) # todo: remove hard-coding, make sigma include off center
  # cor(x)
  
  mu <- b0 + x %*% beta  # equivalent to log odds given by: b0 + beta[1]*x[,1] + beta[2]*x[,2] + beta[2]*x[,2]
  y <- rbinom(n, 1, inverse_logit(mu))
  data.frame(y=y, x=x)
}
# data_generation_simple(10)


data_generation_complex <- function(
    n, # sample size
    maxT # maximum time points (integer)
){
  # create containers for all time-varying quantities
  keep_index    <- rep(1, n*maxT)
  id_container    <- rep(1:n, each=maxT)
  time_container    <- rep(1:maxT, times=n)
  income_container    <- numeric(n*maxT)
  pm_container        <- matrix(data=NA, nrow=n*maxT, ncol=7)
  mortality_container <- numeric(n*maxT) 
  
  # generate baseline confounders
  SES_u <- rnorm(n)
  # expand baseline confounders into containers
  SES_container = rep(SES_u, each=maxT)
  
  container_row = 1 # keep track of the row within a container
  for(i in 1:n){
    # loop over individuals
    dead = 0
    for(j in 1:maxT){
      # loop over time
      if(dead==1){
        keep_index[container_row] = 0
      }
      ### start: generate income
      mu_income = SES_container[container_row] * 200
      if(j==1){
        mu_income = mu_income + 100000
      } else{
        # income stays fairly constant over time
        mu_income = j*1000 + mu_income + income_container[container_row-1]
      }
      incomet = rnorm(1, mean=mu_income, sd=100)
      income_container[container_row] = incomet
      ### end: generate income
      ### start: generate pm componenents
      # generate log-pm component values as multivariate normal
      mu_pm <- rep(0,7) - 100000/1000 + (.01/100)*income_container[container_row]         # todo: use more reasonable value for how income impacts individual componentns
      Sigma_pm <- clusterGeneration::genPositiveDefMat(7)$Sigma   # todo: use empirically generated values
      pmt = rmvnorm(1, mu_pm, Sigma_pm)
      pm_container[container_row,] = pmt
      ### end: generate pm components
      
      ### start: generate mortality
      mu_death = 3*j + 6 + pmt %*% c(.1, .1, 0, 0 , 0, 0, 0) +  0.0001 * incomet # log odds of death at time j
      dead = rbinom(1, 1, inverse_logit(mu_death))
      mortality_container[container_row] = dead
      ### end: generate mortality
      #nothing should happen below this for a given row
      container_row = container_row + 1    
    }
  }
  
data.frame(id=id_container, time=time_container, y=mortality_container, pm=pm_container, income=income_container, ses=SES_container)[which(keep_index==1),]
}

#data_generation_complex(40, maxT=5)


# step 2: analyze a single set of data

# data_analysis_simple: logistic regression model for a small mixture
data_analysis_simple <- function(
    data
){
  glmfit = glm(y~x.1 + x.2 + x.3, data=data, family=binomial(link = "logit"))
  coef_estimates <- coef(glmfit)
  truth <- c(-1.2, 1, -1, 0)
  bias <- coef_estimates - truth
  variance <- diag(vcov(glmfit))
  c(bias = bias, variance=variance)
}
# dat <- data_generation_simple(1000000)
# data_analysis_simple()


data_analysis_complex <- function(
    data
){
  # fill in a dummy method that should be replaced (e.g. your extension to qgcomp)
  glmfit = glm(y ~ pm.1 + pm.2 + pm.3 + pm.4 + time + income, data=data, family=binomial(link = "logit"))
  coef_estimates <- coef(glmfit)
  #truth <- c(-1.2, 1, -1, 0) ? we will have to generate the truth
  #bias <- coef_estimates - truth
  variance <- diag(vcov(glmfit))
  c(coefs = coef_estimates, variance=variance)
  
}
#data = data_generation_complex(40, maxT=5)

# step 3: wrapper function to generate data and analyze it
wrapper_fn <- function(
    generation_function=data_generation, 
    analysis_function=data_analysis, 
    ...
    ){
  dat = generation_function(...)
  analysis_function(dat)
}
# wrapper_fn(data_generation_simple, data_analysis_simple, n=100)
# wrapper_fn(data_generation_complex, data_analysis_complex, n=100, maxT=5)

# step 4: iterate
simulation_iterate <- function(
    iterations=5,                       # number of simulation iterations (ideally 1000+)
    wrapper_function=wrapper_fn,
    ...                                 # arguments passed on to wrapper_function
    ){
  # todo: add in parallel processing via the future.apply::future_lapply function
  res_list = lapply(1:iterations, function(x) wrapper_function(...))
  res = do.call(rbind, res_list)
  res
}



# step 5: summarize simulations across iterations
summmarize_simulations <- function(
    sim_results                        # matrix of results (rows are iterations)
    ){
  meanres = apply(sim_results, 2, mean)
  stdres = apply(sim_results, 2, sd)
  cbind(
    mean = meanres,                   # Mean across all iterations
    #std = stdres,                     # Monte Carlo standard deviation
    var = stdres^2                    # Monte Carlo variance
  )
}
# simres <- simulation_iterate(iterations=1000, wrapper_function=wrapper_fn, generation_function=data_generation_simple, analysis_function=data_analysis_simple, n=1000)
# summmarize_simulations(simres)
# plot(density(simres[,"bias.(Intercept)"]))

# simres <- simulation_iterate(iterations=10, wrapper_function=wrapper_fn, generation_function=data_generation_complex, analysis_function=data_analysis_complex, n=1000, maxT=5)
# summmarize_simulations(simres)
# plot(density(simres[,"bias.(Intercept)"]))
