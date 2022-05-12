# cumulative risk with competing events in R
library(survival)
library(boot)
# some example data
dat = structure(list(
  # confounder
  z = c(0.6404, 3.4271, 1.4669, 2.3247, 1.4541, 
        1.121, 4.3327, 2.2207, 2.4344, 2.4264, 3.1343, 2.2454, 4.5132, 
        0.1525, 2.9201, 3.8828, 0.0527, 0.5397, 2.401, 0.4918, 3.8984, 
        0.0763, 2.0559, 0.1441, 0.4752, 4.3646, 3.4368, 2.319, 2.7647, 
        2.9145, 2.4415, 1.4893, 0.887, 1.7802, 0.484, 4.1863, 3.9136, 
        1.3015, 0.0822, 1.9567, 1.9624, 4.5027, 2.5258, 3.4387, 0.9406, 
        3.7496, 3.1775, 3.0166, 1.598, 1.0061, 0.3052, 0.3272, 0.0219, 
        1.587, 1.6079, 2.5279, 3.438, 0.7045, 3.7353, 0.8537, 1.1939, 
        4.146, 4.4421, 0.4723, 0.5816, 0.6541, 4.5024, 2.1109, 4.009, 
        4.2143, 1.4815, 0.4459, 0.2995, 1.4923, 3.8026, 0.1972, 3.9699, 
        1.1596, 0.0317, 3.7168, 1.5888, 3.9596, 4.0694, 0.2931, 0.9599, 
        4.026, 1.1861, 1.7855, 1.6138, 0.7544, 4.2854, 0.0189, 2.6171, 
        2.8315, 4.8379, 1.9569, 4.9944, 3.5785, 0.3768, 1.1604), 
  # main binary exposure
  x = c(0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 
        0L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L, 
        1L, 0L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 
        0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 
        1L, 1L, 1L, 1L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 
        0L, 1L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L, 
        0L, 0L, 1L, 0L, 0L), 
  # event/censoring time
  t = c(1, 1, 0.408, 1, 1, 1, 1, 
        1, 0.3931, 1, 1, 1, 1, 1, 1, 1, 0.9797, 1, 1, 1, 0.3958, 1, 1, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.6086, 1, 0.6343, 
        1, 1, 1, 0.4194, 1, 1, 1, 1, 1, 0.9595, 1, 1, 0.2245, 1, 1, 1, 
        1, 0.7205, 1, 1, 1, 1, 1, 0.2587, 1, 1, 1, 1, 1, 0.0338, 0.2411, 
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.4594, 1, 1, 1, 1, 1, 1, 1, 
        0.0486, 1, 1, 1, 1, 1, 1, 0.0021, 1), 
  # any event (1) vs censored (0)
  d = c(0, 0, 1, 0, 
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 
        0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0), 
  # censored (0) vs event type 1 (1) vs event type 2 (2)
  # a "multi state" outcome
  d_multi = c(0, 0, 1, 0, 
        0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 1, 0, 
        0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 
        0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
        0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0)), 
  class = "data.frame", row.names = c(NA, -100L
  ))

# recode outcome type as factor for use in coxph function
dat$event <- factor(dat$d_multi, 0:2, labels=c("censor", "d_cause1", "d_other"))
dat$agein <- rep(0, nrow(dat)) # late entry
dat$id <- 1:nrow(dat)


# check event times (may need to implement some fixes in real data)
(sc <- survcheck(Surv(agein, t, event) ~ 1, data = dat, id=id, timefix=TRUE))
#sc$gap
#dat[sort(c(sc$gap$row-1, sc$gap$row)),]

# multistate cox model (equivalent to fitting separate Cox models for event of interest + competing causes)
mscox = coxph(Surv(agein, t, event) ~ x+z, data = dat, id=id)
# set up data ahead of time to look at contrast of exposed vs. unexposed
dat1 = dat0 = dat
# x is your exposure
dat1$x = 1
dat0$x = 0
msmod = survfit(mscox, newdata = dat) # get survival curves, population
msmod1 = survfit(mscox, newdata = dat1) # get survival curves, if exposed
msmod0 = survfit(mscox, newdata = dat0) # get survival curves, if unexposed
tidx = which(msmod$time > .7)[1] # find first time after the time point of interest
#predrisk_cens = cbind(msmod$pstate[tidx,,1], msmod1$pstate[tidx,,1], msmod0$pstate[tidx,,1]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
predrisk_cause1 = cbind(avgrisk = msmod$pstate[tidx,,2], risk_exposed = msmod1$pstate[tidx,,2], risk_unexposed = msmod0$pstate[tidx,,2]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
#predrisk_other = cbind(msmod$pstate[tidx,,3], msmod1$pstate[tidx,,3], msmod0$pstate[tidx,,3]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
# average risk, risk in exposed, risk in unexposed
apply(predrisk_cause1, 2, mean)



# now do bootstrapping (also gives original point estimate)
bootrisk = function(dat, idx){
  keepids = unique(dat$id)[idx]
  newids = 1:length(idx)
  # resample, allowing for multiple obs per individual
  relist = lapply(1:length(keepids), function(x) {
    tdat= dat[dat$id==keepids[x],,drop=FALSE]
    tdat$id = x
    tdat
    }
    )
  bootdat = do.call(rbind, relist)
  mscox = coxph(Surv(agein, t, event) ~ x+z, data = bootdat, id=id)
  dat1 = dat0 = bootdat
  # x is your exposure
  dat1$x = 1
  dat0$x = 0
  msmod = survfit(mscox, newdata = dat) # get survival curves, population
  msmod1 = survfit(mscox, newdata = dat1) # get survival curves, if exposed
  msmod0 = survfit(mscox, newdata = dat0) # get survival curves, if unexposed
  tidx = which(msmod$time > .7)[1] # find first time after the time point of interest
  #predrisk_cens = cbind(msmod$pstate[tidx,,1], msmod1$pstate[tidx,,1], msmod0$pstate[tidx,,1]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
  predrisk_cause1 = cbind(avgrisk = msmod$pstate[tidx,,2], risk_exposed = msmod1$pstate[tidx,,2], risk_unexposed = msmod0$pstate[tidx,,2]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
  #predrisk_other = cbind(msmod$pstate[tidx,,3], msmod1$pstate[tidx,,3], msmod0$pstate[tidx,,3]) # individual predicted risk for cause of interest (times X sample size X event types + 1)
  # average risk, risk in exposed, risk in unexposed
  apply(predrisk_cause1, 2, mean)
}

#
bootstraps = boot(dat, bootrisk, R=1000)

# original estimates
bootstraps$t0
# risk difference 
(rd <- bootstraps$t0["risk_exposed"] - bootstraps$t0["risk_unexposed"])
# number needed to harm/treat
(nnh = 1/rd)


# bootstrap risk difference
(rdb <- bootstraps$t[,2] - bootstraps$t[,3])
quantile(rdb, c(0.025, 0.975))     # percentile method (generally preferred)
rd + qnorm(c(.025, .975))*sd(rdb) # asymptotic normality based method
# bootstrap number needed to harm/treat
quantile(1/rdb, c(0.025, 0.975))     # percentile method (generally preferred)
1/(rd + qnorm(c(.025, .975))*sd(rdb)) # asymptotic normality based method (probably not good to use with the nnh)

hist(rdb)
hist(1/rdb) # not even close to normal

