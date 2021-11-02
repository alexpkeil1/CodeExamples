# interventional analogue of natural effects for a single time-fixed exposure
# and a single time-fixed mediator
library(ggplot2)


#############################################
# R functions for simulating data for testing
#############################################

expit <- function(logodds) 1/(1+exp(logodds))

dgm <- function(n){
  z <- rnorm(n)
  # x takes on values in {0,1,2}
  x <- rbinom(n, 1, expit(-1+z)) + rbinom(n, 1, expit(-1+z))
  # stochastic potential mediators
  m0 <- rnorm(n, -1 + 0 + 0.5*z, 1)
  m1 <- rnorm(n, -1 + 1 + 0.5*z, 1)
  m2 <- rnorm(n, -1 - 2 + 0.5*z, 1)
  m <- ifelse(x==1, m1, ifelse(x==2, m2, m0))
  y00 <- rnorm(n,0,2)
  permute_idx = sample(n)
  yx1m1t <- y00 + 1 + m1
  yx2m2t <- y00 + 2 + 2*m2
  yx0m0t <- y00 + 0 + m0
  # true randomized analogue effects
  yx1m1 <- y00 + 1 + m1[permute_idx]
  yx2m2 <- y00 + 2 + 2*m2[permute_idx]
  yx0m0 <- y00 + 0 + m0[permute_idx]
  #
  yx1m0 <- y00 + 1 + m0[permute_idx]
  yx2m0 <- y00 + 2 + 2*m0[permute_idx]
  yx0m2 <- y00 + 0 + m2[permute_idx]
  yx0m1 <- y00 + 0 + m1[permute_idx]
  y <- ifelse(x==1, yx1m1t, ifelse(x==2, yx2m2t, yx0m0t))
  data.frame(
    z,x,m,y,m0,m1, yx1m1, yx2m2, yx2m0, yx1m0,yx0m2, yx0m1,yx0m0
  )
}


#############################################
# R functions for mediational g-formula with
# 1 time-fixed exposure, 1 time-fixed mediator
#############################################


fit_mediational_gformula <- function(
  dat,   # input data
  M=1000 # Monte Carlo sample size
  ){
  ########################
  # fit models
  ########################
  
  m_m <- glm(m ~ factor(x) + z, data=dat)   # mediator (continuous)
  m_y <- glm(y ~ factor(x)*m + z, data=dat) # outcome (continuous)
  # each of these could also be modeled assuming linear conditional effect of X
  # or could remove/add interactions
  
  ########################
  # MC sampling
  ########################
  MCids = sort(sample(1:nrow(dat), M, replace = TRUE)) # sorting not necessary
  MCdat <- dat[MCids,]
  
  ########################
  # simulate mediators under interventions on X
  ########################
  # Fix X = 1
  # predict full distribution of M | X=1, Z
  pm_a1 <- predict(m_m, newdata = data.frame(         x=1, MCdat[, c("z"), drop=FALSE]), se.fit=TRUE)
  pm_a1 <- pm_a1$fit + rnorm(M, 0, pm_a1$residual.scale)
  
  # Fix X = 2
  # predict full distribution of M | X=2, Z
  pm_a2 <- predict(m_m, newdata = data.frame(         x=2, MCdat[, c("z"), drop=FALSE]), se.fit=TRUE)
  pm_a2 <- pm_a2$fit + rnorm(M, 0, pm_a2$residual.scale)
  
  # Fix X = 0
  # predict full distribution of M | X=0, Z
  pm_a0 <- predict(m_m, newdata = data.frame(         x=0, MCdat[, c("z"), drop=FALSE]), se.fit=TRUE)
  pm_a0 <- pm_a0$fit + rnorm(M, 0, pm_a0$residual.scale)
  
  
  ########################
  # simulate outcomes under interventions on X, M
  ########################
  # two options discussed in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5285457/
  #Q(x*,x) or Q(x*,x*), x = 1,2 x^* = 0
  #option 1:  predict mean of Y | G(M(X=2)|Z), X=0, Z
  #mperm = 1:M
  #option 2:  predict mean of Y | G(M(X=0)), X=0, Z (permute M from marginal simulated distribution)
  mperm = sample(1:M)
  #
  py_a0g0 <- predict(m_y, newdata = data.frame(m=pm_a0[mperm], x=0, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  #py_a0g1 <- predict(m_y, newdata = data.frame(m=pm_a1[mperm], x=0, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  #py_a0g2 <- predict(m_y, newdata = data.frame(m=pm_a2[mperm], x=0, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  
  #Q(x,x) or Q(x,x*), x = 1,2 x^* = 0
  # i.e. Q(0,m(0)) or Q(0,m(1)) or Q(0,m(2))
  #option 1:  predict mean of Y | G(M(X=0)|Z), X=0, Z
  #mperm = 1:M
  #option 2:  predict mean of Y | G(M(X=0)), X=0, Z (permute M from marginal simulated distribution)
  mperm = sample(1:M)
  #
  py_a1g1 <- predict(m_y, newdata = data.frame(m=pm_a1[mperm], x=1, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  py_a1g0 <- predict(m_y, newdata = data.frame(m=pm_a0[mperm], x=1, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  py_a2g2 <- predict(m_y, newdata = data.frame(m=pm_a2[mperm], x=2, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  py_a2g0 <- predict(m_y, newdata = data.frame(m=pm_a0[mperm], x=2, MCdat[, c("z"), drop=FALSE]), se.fit=FALSE)
  
  # randomized analogue effect estimates (output these)
  data.frame(
    rTE1  = mean(py_a1g1) - mean(py_a0g0),
    rTE2  = mean(py_a2g2) - mean(py_a0g0),
    rNDE1 = mean(py_a1g0) - mean(py_a0g0),
    rNDE2 = mean(py_a2g0) - mean(py_a0g0),
    rNIE1 = mean(py_a1g1) - mean(py_a1g0),
    rNIE2 = mean(py_a2g2) - mean(py_a2g0)
  )
}

bootstrap_mediational_gformula <- function(
  dat,     # input data
  iter=200, # number of bootstrap iterations
  ... # additional arguments to fit_mediational_gformula
  ){
  nr = nrow(dat)
  doIter <- function(i, dat){
    idxi = sample(1:nr, nr, replace=TRUE)
    dati = dat[idxi,,drop=FALSE]
    fit_mediational_gformula(dati, ...)
  }
  res = lapply(X=1:iter, FUN = doIter, dat=dat)
  do.call(rbind, res)
}

printEsts <- function(pointests, bootests, alpha=0.05){
  se <- apply(bootests, 2, function(x) sd(x))
  lci <- apply(bootests, 2, function(x) quantile(x, p=alpha/2))
  uci <- apply(bootests, 2, function(x) quantile(x, p=1-alpha/2))
  z = as.numeric(pointests)/se
  x = cbind(
    estimate=as.numeric(pointests), 
    std.error=se, 
    lower.ci=lci, 
    upper.ci=uci,
    Z =z,
    `Pr(>Z)` = (1-pnorm(abs(z), lower.tail = TRUE))*2
  )
  printCoefmat(x, digits = 3, signif.stars = FALSE, P.values = TRUE, has.Pvalue = TRUE)
  invisible(x)
}


#############################################
# data simulation for testing
#############################################
# generate data
dat = dgm(1000)


#############################################
# analysis of simulated data
#############################################
# get point estimate
pointests <- fit_mediational_gformula(dat, M=5000)

# get bootstrap estimates
bootests <- bootstrap_mediational_gformula(dat, iter = 1000, M=5000)

# printing out
printEsts(pointests, bootests, alpha=0.05)

# comparing to large simulation truth
trdat = dgm(100000)
# truths
truths = c(
  true_rTE1  = mean(trdat$yx1m1) - mean(trdat$yx0m0), #rTE  (x=1 vs 0)
  true_rTE2  = mean(trdat$yx2m2) - mean(trdat$yx0m0), #rTE  (x=2 vs 0)
  true_rNDE1 = mean(trdat$yx1m0) - mean(trdat$yx0m0), #rNDE (x=1 vs 0)
  true_rNDE2 = mean(trdat$yx2m0) - mean(trdat$yx0m0), #rNDE (x=2 vs 0)
  true_rNIE1 = mean(trdat$yx1m1) - mean(trdat$yx1m0), #rNIE (x=1 vs 0)
  true_rNIE2 = mean(trdat$yx2m2) - mean(trdat$yx2m0) #rNIE (x=2 vs 0)
)
print(truths)

# visualizing bootstrap distribution
ggplot() + 
  geom_density(aes(x=rTE1, color="Total Effect (X=1 v X=0)"), data=bootests)+
  geom_density(aes(x=rNIE1, color="Indirect Effect (X=1 v X=0)"),data=bootests)+
  geom_density(aes(x=rNDE1, color="Direct Effect (X=1 v X=0)"), data=bootests)+
  geom_vline(aes(xintercept=rTE1, color="Total Effect (X=1 v X=0)"), data=pointests, linetype=3)+
  geom_vline(aes(xintercept=rNIE1, color="Indirect Effect (X=1 v X=0)"), data=pointests, linetype=3)+
  geom_vline(aes(xintercept=rNDE1, color="Direct Effect (X=1 v X=0)"), data=pointests, linetype=3)+
  theme_classic() + 
  scale_color_discrete(name="") + 
  scale_x_continuous(name="Estimate") + 
  scale_y_continuous(name="Bootstrap density") 
  
ggplot() + 
  geom_density(aes(x=rTE2, color="Total Effect (X=2 v X=0)"), data=bootests)+
  geom_density(aes(x=rNIE2, color="Indirect Effect (X=2 v X=0)"),data=bootests)+
  geom_density(aes(x=rNDE2, color="Direct Effect (X=2 v X=0)"), data=bootests)+
  geom_vline(aes(xintercept=rTE2, color="Total Effect (X=2 v X=0)"), data=pointests, linetype=3)+
  geom_vline(aes(xintercept=rNIE2, color="Indirect Effect (X=2 v X=0)"), data=pointests, linetype=3)+
  geom_vline(aes(xintercept=rNDE2, color="Direct Effect (X=2 v X=0)"), data=pointests, linetype=3)+
  theme_classic() + 
  scale_color_discrete(name="") + 
  scale_x_continuous(name="Estimate") + 
  scale_y_continuous(name="Bootstrap density") 



