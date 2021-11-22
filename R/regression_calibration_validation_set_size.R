# regression calibration simple example (loss of variance with a small validation set)

Nstudy = 100       # number in analytic data set
Nvalid = 20       # number in validation data set
me_variance = 0.5 # measurement error variance 

# don't change anything below this
xtrue <- rnorm(Nstudy)
ytrue <- rnorm(Nstudy) + xtrue # true coefficient = 1
xmeasured <- xtrue + rnorm(Nstudy, 0, sqrt(me_variance))
validset <- 1:Nvalid
print(paste("Attenuation coefficient = ", var(xtrue)/var(xmeasured)))
xtrue[-validset] <- NA
# 

reg_calibration <- function(xtrue, xmeasured, ytrue, nboots){
  # estimate + bootstrap variance 
  alpha = lm(xtrue ~ xmeasured)$coefficients[2]
  beta = lm(ytrue~xmeasured)$coefficients[2]
  n = length(xtrue)
  nvalid = sum(!is.na(xtrue))
  rcals = numeric(nboots)
  for(i in 1:nboots){
    idx = sample(1:n, size=n,replace=TRUE)
    idxv = sample(1:nvalid, size=nvalid,replace=TRUE)
    alphab = lm(xtrue[idxv] ~ xmeasured[idxv])$coefficients[2]
    betab = lm(ytrue[idx]~xmeasured[idx])$coefficients[2]
    rcals[i] = betab/alphab
  }
  z = (beta/alpha)/sd(rcals)
  printCoefmat(cbind(estimate = beta/alpha,
             std.error = sd(rcals),
             `t value` = z,
             `Pr(>|t|) ` = 2*(1-pnorm(abs(z)))
               ), signif.stars = FALSE, P.values = TRUE)
}

doall <- function(xtrue, xmeasured, ytrue){
  cat("\nEstimate in validation data set\n")
  print(summary(lm(ytrue~xtrue))$coefficients[2,,drop=FALSE])
  cat("\nEstimate with error prone measure\n")
  print(summary(lm(ytrue~xmeasured))$coefficients[2,,drop=FALSE])
  cat("\nRegression calibration estimate\n")
  reg_calibration(xtrue, xmeasured, ytrue, nboots=200)
  
}

doall(xtrue, xmeasured, ytrue)
