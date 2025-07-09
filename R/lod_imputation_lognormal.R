# Imputation of log-normally distributed variables


# multiple imputation based on Lubin et al (2004)
# as well as a generalization to more standard multiple imputation schemes
# this can also be implemented within the `qgcomp' package via the
# mice.impute.tobit function

library(survival)


# function to do a single imputation of new values below the LOD
# note: y is assumed lognormally distributed
draw_tobit <- function(y, lod, x=NULL){
  # draw from log-normal distribution to impute y values below an LOD
  ry = !is.na(y)
  wy = !ry
  nmiss = sum(wy)
  if (length(lod)==1){
    LOD = rep(lod, length(y))
  } else  LOD = lod

  if(!is.null(x)){
    x = cbind(rep(1,length(y)), x)
  } else
    x = cbind(rep(1,length(y)))
  ylod = y
  ylod[wy] = LOD[wy]
  # fit distribution
  fit <- survreg(Surv(time = ylod, event = ry, type = "left") ~ -1+x, dist = "lognormal", control = survreg.control())
  # draw parameters from posterior
  draw = qgcomp:::.rmvnorm(1, c(fit$coefficients, `Log(Scale)` = log(fit$scale)), fit$var)
  # refit to get linear predictor under parameter draw
  fit2 <- survreg(Surv(time = ylod, event = ry, type = "left") ~ -1+x, dist = "lognormal", init = draw, control = survreg.control(maxiter = 0))
  fit2$linear.predictors[wy] = ifelse(is.na(fit2$linear.predictors[wy]), -20, fit2$linear.predictors[wy])
  # draw
  fub = plnorm(LOD[wy], fit2$linear.predictors[wy], fit2$scale)
  u <- runif(nmiss) * fub
  returny <- qlnorm(u, fit2$linear.predictors[wy], fit2$scale)
  impy = y
  impy[wy] = returny
  impy
}


# Example 
# first, generate data subject to LOD for the outcome
dgm <- function(n){
  x <- rnorm(n,0,1)
  y <- rnorm(n,x,1)
  lod <- 0
  yobs = y
  lodidx = which(yobs<lod)
  yobs[lodidx] <- NA
  ylod <- yobs
  ylod[lodidx] <- lod
  data.frame(x,y,yobs,ylod, lod)
}
dat = dgm(500)

M = 100



# Lubin's approach (no covariates used in imputation): will be biased if using imputed values in a regression
margimps = lapply(1:M, function(i) 
  # create M imputed datasets
  log(draw_tobit(exp(dat$yobs), exp(dat$lod), x=NULL))
  )
margfits <- lapply(margimps, function(vals) 
  # apply linear regression to every imputed sample
  lm(yimp ~ x, data = data.frame(yimp=vals, x=dat$x))
  )
# extract coefficients and variances
coef_margimps = sapply(margfits, function(x) coef(x)[2])
var_margimps = sapply(margfits, function(x) vcov(x)[2,2])
# Rubin's rules for estimate and standard error
marg_est = mean(coef_margimps)
marg_se = sqrt(var(coef_margimps)*(M+1)/(M-1) + mean(var_margimps))





# Conditional approach (covariates used in imputation): helps avoid null bias
condimps = lapply(1:M, function(i) 
  # create M imputed datasets
  log(draw_tobit(exp(dat$yobs), exp(dat$lod), x=dat$x))
  )
condfits <- lapply(condimps, function(vals) 
  # apply linear regression to every imputed sample
  lm(yimp ~ x, data = data.frame(yimp=vals, x=dat$x))
  )
# extract coefficients and variances
coef_condimps = sapply(condfits, function(x) coef(x)[2])
var_condimps = sapply(condfits, function(x) vcov(x)[2,2])
# Rubin's rules for estimate and standard error
cond_est = mean(coef_condimps)
cond_se = sqrt(var(coef_condimps)*(M+1)/(M-1) + mean(var_condimps))



# bias for the mean
mean(sapply(margimps, mean))
mean(sapply(condimps, mean))

# bias for regression coefficients
bias_cond = mean(coef_condimps) - 1
bias_marg = mean(coef_margimps) - 1

bias_marg
bias_cond
# standard error (Rubin's rules estimate)
marg_se
cond_se

# estimate with true values
trfit = lm(y ~ x, data = dat)
true_est = coef(trfit)[2]
true_se =  sqrt(vcov(trfit)[2,2])


# estimate and 95% confidence intervals
marg_est + qnorm(c(0.5, 0.025, 0.975))*marg_se
cond_est + qnorm(c(0.5, 0.025, 0.975))*cond_se
true_est + qnorm(c(0.5, 0.025, 0.975))*true_se




plot(density(coef_margimps), xlim = c(0.5, 1.5), main=NA)
lines(density(coef_condimps), col="red")
abline(v=1, lty=3)
legend("bottomleft", c("Marginal imp.", "Conditional imp."), col=c("red", "black"), lty=1)

