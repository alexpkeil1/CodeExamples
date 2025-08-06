# simulation based power for logistic regression with a continuous, lognormally distributed exposure
######################################################################################################################
# Author: Alex Keil
# Program: case_control_lognormal_exposure_power.R
# Language: R
# Date: Aug 2025
# Project: Coding examples
# Tasks: Calculate power for a crude logistic regression analysis of a case-control study
# Description: This calculates simulation based power for a case-control study with a 
#               single, log-normally distributed exposure. Since closed-form estimators of
#               power generally require strong assumptions, this is an improved way to 
#               calculate power. It could be extended to include covariates/confounders
#               or other design elements. The downside of simulation based power is the 
#               computational time required
# Keywords: power, simulation
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

# two packages to install via: install.packages(c("future", "future.apply")) to help with parallel processing
library(future)
library(future.apply)


dgm <- function(trueor = 1.1,         # true odds ratio for a 1 unit increase in exposure
                ncases=5000,          # expected number of cases in the study
                ncontrols = 500000,   # expected number of non-cases in the study
                logmean=-2,           # exposure distribution parameter: mean of log(x)
                logsd=0.75,           # exposure distribution  parameter: sd of log(x)
                xmax = 10             # exposure distribution  parameter: maximum exposure (log-normally distribted exposures can sometimes yield really large, unrealistic values, and this caps those)
                ){
  n = ncases + ncontrols
  popodds = (ncases/ncontrols) # log odds in the study population
  xmean = exp(logmean + logsd^2/2)
  #expected_cases = plogis(mean_logodds)*n
  base_odds = popodds/((trueor)^xmean)
  x = rlnorm(n, logmean, logsd)
  while(any(x > xmax)){
    x[which(x>xmax)] =  rlnorm(sum(x>xmax), logmean, logsd)
  }
  y = rbinom(n, 1, base_odds*trueor^x)
  data.frame(x,y)
}


analysis <- function(dat){
  ft <- glm(y~x, data = dat, family=binomial())
  var = vcov(ft)[2,2]
  z =coef(ft)[2]/sqrt(var)
  abs(z)
}


power_logistic <- function(alpha = 0.05,     # type 1 error rate
                           ...,              # this will pass arguments to the dgm function
                           iterations= 1000, # number of simulation iterations (higher is more accurate, but takes longer - 1000 is usually reasonable)
                           verbose = TRUE    # print out power and estimate of computational time
                           ){
  r1 = system.time({
    dat = dgm(...)
    zs1 = analysis(dat)
  })
  qx = quantile(dat$x)
  if(verbose){
    #cat("Quantiles of exposure in example dataset \n")
    #print(qx)
    exptime = as.numeric(r1[3])*iterations/60
    cat("\nPower calculation expected to take up to", round(exptime, 2), "minutes on a single core\n")
  }
  zstat = numeric(iterations)
  
  zs2toend = future_lapply(2:iterations, function(x){
    dat = dgm(...)
    absz = analysis(dat)
  }, future.seed=TRUE)
  zs = c(zs1, do.call(c,zs2toend))
  zcrit = qnorm(1-0.05/2)
  res = list(power=mean(zs > zcrit))
  if(verbose){
    #cat("\nPower at alpha =", alpha,"(%):", round(res*100))
  }
  attr(res, "zvalues") <- zs
  attr(res, "qx") <- qx
  attr(res, "alpha") <- alpha
  class(res) <- "simpow"
  res
}

# this function is just to make the printing of the power look nicer
print.simpow <- function(x,...){
  cat("\nPower:", round(x[["power"]]*100), "%")
  cat("\nalpha =", attr(x, "alpha") )
  cat("\nNumber of simulations: ", length(attr(x, "zvalues")))
  cat("\nExample exposure quantiles (single simulation):\n")
  print(attr(x, "qx"))
}

# example of calculating power for a case-control study with:
# 5000 cases
# 500000 controls
# exposure with a log normal distribution with parameters -2, 0.75 (mean, sd of the log-transformed values)
# maximum exposure of 10
# 1000 iterations (more will give better accuracy, but take longer)
# set parallel processing, see help for "plan" function for details - this version uses half of the available workers
future::plan(strategy=multisession, workers = availableCores(constraints = "connections-16")/2) 
message("Number of parallel processes: ", nbrOfWorkers())
pow_1 = power_logistic(trueor = 1.1, ncases=5000, ncontrols = 500000, logmean=-2, logsd=0.75, xmax = 10, iterations = 1000, verbose=TRUE)
pow_1
