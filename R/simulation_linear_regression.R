
# simple simulation framework, linear regression
# based on the following ideas in simulations:
#  - first generate data where the truth is known
#  - analyze those data with a method of interest (e.g. estimate a parameter), 
#     compare to the truth
#  - do this repeatedly - bias equals the average difference between the estimate
#     and the truth
#    
library(mvtnorm)
library(future)
library(future.apply)

######### Data generating mechanism
#  This generates data under a true linear model for Y, given some exposures
#  useful to use R functions so that this is easy to do repeatedly. Note that
#  this chooses a very specific distribution of the exposure (x) and 
#  covariate (z), and adaptation
#  to other settings requires generating X differently (e.g. binary variables)
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

# test it
data = dgm_linmod(n=15)

######### Data analysis methods
# it's helpful to put all analysis methods into one function that just needs the data
# Here, I've also included "bias" of each approach, but you can also just output
# the raw estimates and look at bias later
analyze <- function(data){
  # this function can be used to do the analysis/analyses you would like to compare
  mod_adj = lm(y ~ x + z, data = data)
  mod_crude = lm(y ~ x, data = data)
  # coefficients
  coef_adj = as.numeric(coef(mod_adj))
  coef_crude = as.numeric(coef(mod_crude))
  # variance/standard errors
  var_adj = as.numeric(diag(vcov(mod_adj)))
  var_crude = as.numeric(diag(vcov(mod_crude)))
  # now output the things you want to look at (here, just )
  truth = attr(data, "truth") # get the truth back
  trueIntercept = truth[1]
  trueX = truth[2]
  # output a named vector (note that we don't use the return() function - an
  #  r function will automatically return the last object)
  c( 
    # bias of intercept
    # single iteration is not "bias", but averaging this over many iterations gives you estimated bias
    biasIntercept_adj = coef_adj[1] - trueIntercept, 
    biasIntercept_crude = coef_crude[1] - trueIntercept,
    # bias of effect estimate of interest
    biasX_adj = coef_adj[2] - trueX,
    biasX_crude = coef_crude[2] - trueX,
    # 95% confidence interval coverage (should equal 95% for unbiased estimator)
    coverage_adj = as.numeric(
        (coef_adj[2] + qnorm(.025)*sqrt(var_adj[2])) < trueX & 
        (coef_adj[2] + qnorm(.975)*sqrt(var_adj[2])) > trueX ),
    coverage_crude = as.numeric(
        (coef_crude[2] + qnorm(.025)*sqrt(var_crude[2])) < trueX & 
        (coef_crude[2] + qnorm(.975)*sqrt(var_crude[2])) > trueX ),
    # power (or type 1 error under the null)
    power_adj = as.numeric(
        (abs(coef_adj[2])/sqrt(var_adj[2]) > qnorm(0.975))),
    power_crude = as.numeric(
        (abs(coef_crude[2])/sqrt(var_crude[2]) > qnorm(0.975)))
  )
}


######### test both functions
dat = dgm_linmod(n=100, beta = c(1,1,1), xzcorr = .8, sigma=2.0)
summary(dat)
cor(dat$x, dat$z) # should be around the true xzcorr value 
analyze(dat) # with only one iterations, these will look biased except with very large n


######## now analyze these repeatedly
# first we need a "wrapper" function that does the data generation
#  and analysis all at once (for purposes explained below, this
#  function needs to have at least one parameter that corresponds
#  to an iteration number) 
sim_wrap <- function(iter, beta = c(1,1,1)){
  dat = dgm_linmod(n=100, beta = beta, xzcorr = .8, sigma=2.0)
  results = analyze(dat)
  results = c(iter=iter, results)
  results
}
# one strategy is to create a loop where we do this over and over and save the 
#  results at each iteration - this has historically been very slow in R, but
# may be faster than alternative approaches, depending on the problem:

do_simulations_loop <- function(numiters){
  result1 = sim_wrap(1) # run it once to figure out storage size
  allresults = matrix(nrow=numiters, ncol = length(result1))
  colnames(allresults) = names(result1)
  allresults[1,] = result1
  for(i in 2:numiters){
    allresults[i,] = sim_wrap(i)
  }
  allresults
}



# alternatively, we can use the "apply" family of functions

# e.g. sapply will output a list where each column corresponds to a result from 
# a single simulation - here we repeat them 100 times (good to do more than that!)
# test to see what this looks like

sim_results1 = sapply(1:1000, sim_wrap, beta=c(0,0,1)) # note that this gives each iteration result in a column
sim_results2 = do_simulations_loop(100)# note that this gives each iteration result in a row

# the first few rows of the simulated results: every row is
head(sim_results2)


##### checking run-time to see if you can pick most efficient code
  # can compare times of the apply and looping approaches
  # look here for meaning of output of "system.time"
  system.time(sapply(1:500, sim_wrap))
  system.time(do_simulations_loop(500))
  
  # you can sometimes buy some extra time by running simulations in parallel
  # here's the easy way
  future::plan(multicore) # setup has some overhead
  future::nbrOfWorkers() # check number of "worker" processes that you can use in parallel
  system.time(sim_results3 <- future_sapply(1:500, sim_wrap, future.seed=TRUE))
    


# now we get the simulation means using the apply function (it "applies" the same 
  #  function to every row in sim_results: here that is the mean)
  # e.g. bias = mean(estimated - true)
  # we could also get standard deviations/variance if interested in something like mean squared error
  print(apply(sim_results1, 1, mean), 2) # rowwise
  print(apply(sim_results2, 2, mean), 2) # columnwise 

# Is any approach approximately unbiased?
# does the confidence interval coverage look appropriate?
# which approach has more statistical power?
# do the results change if you run multiple times? (seed value)
# mean squared error = bias^2 + variance. Can you modify the "analyze" function to get that for both models?
# what happens as you increase the number of iterations for unbiased estimator (bias -> 0)
  
