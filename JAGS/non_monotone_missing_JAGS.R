library("tibble")
library("dplyr")
library("R2jags")
library("coda")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

##########################
# Define some parameters #
##########################

# Expected marginal probabilities of each missingness pattern
ePR2 = 0.35 # 0.20
ePR3 = 0.20 # 0.15
ePR4 = 0.15 # 0.10
# Set parameters of missingness models
gamma_21 <- log(1.5)
gamma_22 <- log(.5)
gamma_31 <- log(1.2)
gamma_41 <- log(2)
gamma_42 <- log(0.5)

#################
# Generate data #
#################

set.seed(5)
n <- 200

# Full data
full <- tibble( 
  id = c(1:n),
  z = rbinom(n, size=1, prob=0.5) ,
  x = rbinom(n, size=1, prob=0.2 + 0.2*z),
  py = 1/(1+exp(-1*(log(0.1/0.9) + log(2)*z))),
  y0 = rbinom(n, size=1, prob=py),
  y1 = rbinom(n, size=1, prob=py*2),
  y = x*y1 + (1-x)*y0)

#lapply(full, mean)

# Calculate intercepts for missingness models
# When example is static, change these to large pop means
gamma_20 <- with(full, -log(1/ePR2-1) - gamma_21*mean(x) - gamma_22*mean(z))
gamma_30 <- with(full, -log(1/ePR3-1) - gamma_31*mean(x))
gamma_40 <- with(full, -log(1/ePR4-1) - gamma_41*mean(x) - gamma_42*mean(y))

# Data with missingness
withmiss <- full %>%
  mutate(pR2 = plogis(gamma_20 + gamma_21*x + gamma_22*z),
         pR3 = plogis(gamma_30 + gamma_31*x ),
         pR4 = plogis(gamma_40 + gamma_41*x + gamma_42*y),
         pR1 = 1 - pR2 - pR3 - pR4) %>%
  rowwise() %>%
  mutate(R = sample(c(1,2,3,4),size=n(),replace=T,prob=c(pR1,pR2,pR3,pR4)),
         z = ifelse(R %in% c(3,4),NA,z),
         y = ifelse(R %in% c(2,3),NA,y)) %>%
  select(id, z, x, y, R) %>%
  arrange(R)

prop.table(table(withmiss$R))

################################################################################
# constraint
################################################################################
SIGMASTAR = 0.01


################################################################################
#
# jags/bugs implementation (no need for rejection sampling step)
#
################################################################################
# data prep
dat = as.list(withmiss[,c("R")])
dat$L = as.matrix(withmiss[,c("z", "x", "y")])
dat$N = length(dat$R)
dat$L[is.na(dat$L)] = -9999 # stan doesn't accept NA, so set all NA to Inf
dat$sigmastar = SIGMASTAR # constraint
# some helpers for the JAGS quirks
dat$z = rep(1, dat$N)
dat$Nc = sum(dat$R==1)
dat$onesc = rep(1, dat$Nc)


jmod <- function() {
  for(i in 1:N){
    z[i] ~ dbern(pmiss[i,R[i]]) # z = 1 for all obs
    #
    logit(pmiss[i, 2]) <- g0[1] + gamma[1]*L[i,1] + gamma[2]*L[i,2]                   
    logit(pmiss[i, 3]) <- g0[2] +                   gamma[3]*L[i,2]                   
    logit(pmiss[i, 4]) <- g0[3] +                   gamma[4]*L[i,2] + gamma[5]*L[i,3] 
    #
    pmiss[i,1] <- 1 - pmiss[i,2] - pmiss[i,3] - pmiss[i,4]
  }
  # constraint (first Nc of N observations must be complete cases)
  for (i in 1:Nc){
    onesc[i] ~ dbern(C[i])
    C[i] <- step(pmiss[i,1]-sigmastar)
  }
  
  for(k in 1:3){
    g0[k] ~ dnorm(0,0.01)
  }
  for(k in 1:5){
    gamma[k] ~ dnorm(0,0.01)
  }
}


# specify parameters to track
parms <- c('gamma', 'g0')
# initial values
jagschains=4
initlist = lapply(1:stanchains,
                  #crucial to make this consistent with a valid first estimate of pi
                  function(x) list(
                    g0=runif(3, min=-7,max=-3), gamma=runif(5, min=-.05,max=.05)
                  )
) 



set.seed(1) #MCMC generally is not deterministic, so use a seed for reproducibility (but not testing)
bayes.m1 <- R2jags::jags(data=dat, inits=initlist,
                         parameters.to.save=parms,
                         n.chains=jagschains, n.iter=20000, n.burnin=0,n.thin=1,
                         model.file=jmod)
bayes.m1
#ft2 # versus Stan
jagsfit <- as.mcmc(bayes.m1)
#traceplot(jagsfit)

# convert to data frame
res = data.frame(as.matrix(jagsfit), check.names = FALSE)
burnin = bayes.m1$BUGSoutput$n.burnin
iters = bayes.m1$BUGSoutput$n.iter
nchains = bayes.m1$BUGSoutput$n.chains
res$iter = rep((burnin+1):iters, times=nchains)
res$chain = rep(1L:nchains, each=iters-burnin)

ggplot()+ theme_classic() + scale_color_discrete(name="") +
  geom_line(aes(y=`gamma[1]`, x=iter, color=factor(chain)), data = res, size=0.1)
