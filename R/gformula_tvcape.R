######################################################################################################################
# Author: Alex Keil
# Program: gformula_tvcape.R
# Language: R
# Date: Wednesday, January 22, 2014 at 6:14:27 PM
# Project: Coding examples
# Tasks: Simulations using the G-formula over a 2 time point study with 1 baseline and 1 (possibly unmeasured) time-
#        varying confounder affected by prior exposure (TVCAPE)
# Data in: NA (simulated in programf)
# Description: Several R functions to generate (binary) data over two time points and analyze according to the 
#              parametric g-formula algorithm (though with adequate sample size can be done non-parametrically)
#             Note that this program does not implement variance estimation using bootstrapping
# Keywords: g-formula, unmeasured confounding, simulation, sensitivity analysis, toy analysis
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################

# data generation ----
#  # DAG:
#  {v, w1} -> {x1}
#  {v, w1, x1} -> {y1}
#  {v, w2} -> {x2}
#  {v, w1, x1, x2} -> {y2}
#
makedata <- function(N=100, tvcape=TRUE){
  ntimes = 2 # hard-coded number of time points
  id = rep(1:N, each=ntimes)
  time = rep(1:ntimes, N)
  
  v = rbinom(N, 1, 0.5)
  w1 = rbinom(N, 1, 0.5)
  x1 = rbinom(N, 1, 1/(1+exp(-(-2 + w1 + v)))) 
  y1 = rbinom(N, 1, 1/(1+exp(-(-1 + x1 + w1 + v))))
  
  if(tvcape) {
    w2 = rbinom(N, 1, 1/(1+exp(-(-1 + x1 + v))))
  } else{
    w2 = rbinom(N, 1, 0.5)
  }
  x2 = rbinom(N, 1, 1/(1+exp(-(-2 + w2 + v)))) 
  y2 = rbinom(N, 1, 1/(1+exp(-(-2 + 3*x1 + x2 + w2 + v))))
  cx = x1 + x2
  cu = w1 + w2
  d = apply(cbind(y1,y2), 1, max)
  # "flat" data
  dat = data.frame(v, w1, w2, x1, y1, x2, y2, cx, cu, d)
  # person period data (this is an ugly trick for 2 time points)
  ppdat = data.frame(id, time, v=c(rbind(v, v)), x=c(rbind(x1, x2)), xlag=c(rbind(x1*NA, x1)), w=c(rbind(w1, w2)), y=c(rbind(y1, y2)))
  DROP = c(rbind(ifelse(y1,FALSE,FALSE), ifelse(y1,TRUE,FALSE)))
  list(widedat = dat, longdat = ppdat[!DROP,,drop=FALSE])
}



# modeling ----
gmod <- function(D=dat, tvcape=TRUE){
  #D = makedata(N=100)$longdat
  # play around with these to see impacts of model misspecification
  # time 1
  mx1 = glm(x ~ w + v, data=D[D$time==1,,drop=FALSE], family=binomial())
  my1 = glm(y ~ x + w + v, data=D[D$time==1,,drop=FALSE], family=binomial())
  # time 2
  if(tvcape) {
    mw2 = glm(w ~ xlag + v, data=D[D$time==2,,drop=FALSE], family=binomial())
  } else{
    mw2 = glm(w ~ 1, data=D[D$time==2,,drop=FALSE], family=binomial())
  }
  mx2 = glm(x ~ w + v, data=D[D$time==2,,drop=FALSE], family=binomial())
  my2 = glm(y ~ x + xlag + w + v, data=D[D$time==2,,drop=FALSE], family=binomial())
  res = list(mx1=mx1, my1=my1, mw2=mw2, mx2=mx2, my2=my2)
  return(res)
}

# simulation/prediction ----
predstep <- function(D=dat, MCSize=1000, int=1, moddat=gmoddat){
  # simulate in target population
  maxtime = 2 # hardcode this into this example
  # first create emtpy data frame for target population
  allbase = which(D$time==1)
  keepobs = sample(allbase, size = MCSize, replace=TRUE)
  keepobs2 = c(rbind(keepobs, keepobs)) # another dirty trick for 2 time points
  outdf = data.frame(
    id = rep(1:MCSize, each=maxtime), 
    time = rep(1:maxtime, MCSize), 
    x=ifelse(is.null(int), NA, int), 
    xlag=NA, 
    w=D[keepobs2,"w"],
    v=D[keepobs2,"v"],
    py=NA, # p(Y=1), for survival outcome only need predicted probability (need to fully simulate outcome otherwise)
    y=NA, # p(Y=1), for survival outcome only need predicted probability (need to fully simulate outcome otherwise)
    ylag=0 # p(Y=1), for survival outcome only need predicted probability (need to fully simulate outcome otherwise)
  )
  outdf$w = ifelse(outdf$time==2, NA, outdf$w)
  # loop over observations in simulated dataset (not always efficient in R, especially when repeatedly taking rows of the data.frame)
  for(i in 1:nrow(outdf)){
    # TIME 1 ----
    if(outdf$time[i] == 1){
      # SIMULATE EXPOSURE ----
      #intervene, time 1 (simulate from natural course if int=NULL)
      if(is.null(int)) {
        px1 = predict(moddat$mx1, type="response", newdata=outdf[i,,drop=FALSE])  
        outdf$x[i] = rbinom(1, 1, px1)
      } else{
        outdf$x[i] = int
      }
      # SIMULATE OUTCOME ----
      outdf$py[i] = predict(moddat$my1, type="response", newdata=outdf[i,,drop=FALSE])
      outdf$y[i] = rbinom(1, 1, outdf$py[i])
    }
    # TIME 2 ----
    if(outdf$time[i] == 2){
      outdf$xlag[i] = outdf$x[i-1]
      outdf$ylag[i] = outdf$y[i-1]
      # SIMULATE TIME VARYING CONFOUNDER ----
      pw2 = predict(moddat$mw2, type="response", newdata=outdf[i,,drop=FALSE])  
      outdf$w[i] = rbinom(1, 1, pw2)
      # SIMULATE EXPOSURE ----
      #intervene, time 1 (simulate from natural course if int=NULL)
      if(is.null(int)) {
        px2 = predict(moddat$mx2, type="response", newdata=outdf[i,,drop=FALSE])  
        outdf$x[i] = rbinom(1, 1, px2)
      } else{
        outdf$x[i] = int
      }
      # SIMULATE OUTCOME ----
      outdf$py[i] = predict(moddat$my2, type="response", newdata=outdf[i,,drop=FALSE])
      outdf$y[i] = rbinom(1, 1, outdf$py[i])
    }
    
  }
  outdf
}

predstep_en <- function(...){
  alw <- predstep(int=1,...)
  nev <- predstep(int=0,...)
  nev$id = nev$id+max(alw$id)
  dat = rbind(alw,nev)
}


summary_gdat <- function(data, usepy=FALSE){
  message("average hazards at t1, t2 and cumulative incidence at t2")
  # hazards
  if(usepy){
    m1 = mean(data[data$time==1,"py"])
    m2 = mean(data[data$time==2,"py"]) 
    everY = data[data$time==1,"py"] + (1-data[data$time==1,"py"])*(data[data$time==2,"py"])
  }
  if(!usepy){
    m1 = mean(data[data$time==1,"y"])
    m2 = mean(data[data$time==2,"y"]) 
    everY = tapply(data$y, data$id, max)
  }
  # cumulative incidence (kaplan meier estimator)
  r = mean(everY)
  list(h1=m1, h2=m2, ci=r)
}

# worked example with single simulated data set ----
require("survival")
popsize=800
dat = makedata(N=popsize)$longdat
head(dat)
# g-formula model fitting ----
mods = gmod(D=dat, tvcape=TRUE)
# g-formula simulations ----
# natural course (can take a while) MCSize should ideally be large enough to reduce simualation error to negligible
simdat = predstep(D=dat, MCSize = 10000, int=NULL, moddat=mods)
#interventions: always or never exposed (can take a while)
simdat2 = predstep_en(D=dat, MCSize = 10000, moddat=mods)
# interventions compared to natcourse
simdat$intervention = 0
simdat2$intervention = 1
simdat3 = rbind(simdat, simdat2[simdat2$x==0,])

# compare hazard ratio estimates for different interventions ----
# crude
coxph(Surv(time,y)~x, data=dat, method="breslow")                                      #observed, crude
coxph(Surv(time,y)~x, data=simdat[simdat$ylag==0,], method="breslow")                  #natural course, crude
# inappropriate adjustment for TVCAPE
coxph(Surv(time,y)~x+w+v, data=dat, method="breslow")                                  # observed, adjusted
coxph(Surv(time,y)~x+w+v, data=simdat[simdat$ylag==0,], method="breslow")              #natural course, adjusted
# marginal estimates from g-computation
coxph(Surv(time,y)~x, data=simdat2[simdat2$ylag==0,], method="breslow")                #intervention (causal estimate: always vs. never exposed)
coxph(Surv(time,y)~intervention, data=simdat3[simdat3$ylag==0,], method="breslow")     #intervention (causal estimate: natural course vs. never)

# get discrete time hazards
summary_gdat(dat)              # observed
summary_gdat(simdat)           # natural course
summary_gdat(simdat, usepy=TRUE)           # better estimate using predicted probabilities of binary outcome
summary_gdat(simdat2)          # 2 cohorts combined: 1 exposed, 1 unexposed
summary_gdat(simdat2, usepy=TRUE)           # better estimate using predicted probabilities of binary outcome
# differences between simdat and simdat2 have a lot to do with exposure prevalence
mean(dat$x)
mean(simdat$x)
mean(simdat2$x)
