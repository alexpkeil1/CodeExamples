# bayesian logistic regression via polya gamma latent variables
# this is the Bayesian logit model to use!
#library(devtools)
#install_github("jwindle/BayesLogit", subdir="Code/BLPackage/BayesLogit")
#library(BayesLogit) # did not get completely installed, but only needed the polya gamma sampler (had to kludge together install from github download)
library(mvtnorm)
library(arm)
library(devtools)
library(coda)
#install_github('jgscott/helloPG')
#library(helloPG)
#library(pgdraw)
#poly gamma distribution
#rpg.devroye(1, 1, 0)
n = 80
p = 20
X <- rmvnorm(n, rep(0, p), diag(rep(1, p)))
beta = sample(c(1,rep(0, 14)), p, replace=TRUE)
y <- rbinom(n, 1, 1/(1+exp(0.5 - X %*% beta)))

#checking poly gamma rngs
pt = matrix(NA, ncol=3, nrow = 10000)
#for(i in 1:dim(pt)[1]){
  pt[,1] = BayesLogit::rpg(dim(pt)[1], 1, .3)
  pt[,2] = helloPG::helloPG(dim(pt)[1], .3)$draws[1]
  pt[,3] = pgdraw::pgdraw(rep(1, dim(pt)[1]), rep(0.3, dim(pt)[1])) # only one installable from CRAN
  
#n = 1e7
#system.time(BayesLogit::rpg(n, 1, .3))
#system.time(pgdraw::pgdraw(rep(1, n), rep(0.3, n))) # >100% speedup
  
  
#}
plot(density(pt[,1]))
lines(density(pt[,2]), col='red') # this one is nonsense
lines(density(pt[,3]), col='red')

gibbs.logit <- function(iter = 10000, p=p, prior.sig=0.8, burnin=100, package=c("pgdraw", "BayesLogit")){
  package = package[1]
  cat(paste("Using", package, "to draw from Polya-Gamma distribution\n"))
  B = matrix(NA, nrow = iter, ncol = p+1)
  # auxilliary
  om = 0
  # beta
  beta = runif(p+1)*4-2 + 20
  B[1,] = beta
  # hyper priors
  Bmu = rep(0, p+1)
  Bsig = diag(c(1e4^2, rep(prior.sig^2, p)))
  n = nrow(X)
  Xi = cbind(rep(1,n), X)
  tXi = t(Xi)
  invBsig = gnm::MPinv(Bsig)
  for(i in 2:iter){
    mu = Xi %*% beta
    if(package=="BayesLogit") {
      om = BayesLogit::rpg(n, 1, mu)
    }else om = pgdraw::pgdraw(rep(1, n), mu)
    V = gnm::MPinv(tXi %*% diag(om) %*% Xi + invBsig)
    kappa = (y-1/2)
    m = V %*% (tXi %*% kappa + invBsig %*% Bmu)
    beta  = t(mvtnorm::rmvnorm(1,m,V))
    B[i,] = beta
  }

  data = B[-c(1:burnin),]
  attr(data, "mcpar") <- c(start=burnin+1, end=iter, thin=1)
  attr(data, "class") <- "mcmc"
  data
}


#bayes
#system.time(gl <- gibbs.logit(iter = 2000, p=p, prior.sig=1, burnin=100, package="pgdraw"))
#system.time(gl2 <- gibbs.logit(iter = 2000, p=p, prior.sig=1, burnin=100, package="BayesLogit"))

bayeslogit <- gibbs.logit(iter = 20000, p=p, prior.sig=.75, burnin=1000, package="pgdraw")
summary(bayeslogit)
#plot(bayeslogit)
autocorr(bayeslogit)

# armbayes
bayes.t = bayesglm(y ~ X, family=binomial(), prior.df=Inf, prior.scale = 1, prior.scale.for.intercept = 1e4, scaled=FALSE)
# mle
mle = glm(y ~X, family=binomial())
# pg bayes
bayes.pg = apply(bayeslogit, 2, mean)

cbind(truth=c(-0.5, beta), bglm=bayes.t$coefficients, bayes.pg, mle=mle$coefficients)

