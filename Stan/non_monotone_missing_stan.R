library("rstan")
library("tibble")
library("dplyr")
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
# stan implementation (utilizes post fitting rejection step)
# This is a really iffy implementation
#
################################################################################

LbyR = lapply(unique(withmiss$R),function(x) as.data.frame(withmiss[withmiss$R==x,c("id", "z", "x", "y")]))
whichobsbyR = lapply(LbyR, function(x) which(!is.na(x[1,]))) # just to see which variables in L are present for each unique value of R
whichobsbyR

# data prep
dat = as.list(withmiss[,c("R")])
dat$L = as.matrix(withmiss[,c("z", "x", "y")])
dat$N = length(dat$R)
dat$NR = length(unique(dat$R)) # number of patterns
dat$L[is.na(dat$L)] = -9999 # stan doesn't accept NA, so set all NA to Inf
totparms = sapply(whichobsbyR, length) # number of gamma parameters per model
dat$g = sum(totparms) # total number of gamma parameters
dat$ncols = ncol(dat$L) # total number of variables in L
dat$sigmastar = SIGMASTAR # constraint

stanmod <- "
data{
  real sigmastar;
  int<lower=0> N;
  int<lower=0> NR;
  int<lower=0> g;
  int<lower=0> ncols;
  matrix[N,ncols] L;
  int R[N];
}
transformed data{
  vector[NR-1] ones;
  for(i in 1:(NR-1))
    ones[i] = 1.0;
}
parameters{
 real g0[NR-1];
 real gamma[g - NR - ncols];
}
transformed parameters{
}
model{
  g0 ~ normal(0,10);
  gamma ~ normal(0,10);
  { // local
    matrix[N,NR] pmiss;
    pmiss[:, 2] = inv_logit( g0[1] + gamma[1]*L[:,1] + gamma[2]*L[:,2]                   );
    pmiss[:, 3] = inv_logit( g0[2] +                   gamma[3]*L[:,2]                   );
    pmiss[:, 4] = inv_logit( g0[3] +                   gamma[4]*L[:,2] + gamma[5]*L[:,3] );
    pmiss[:, 1] = 1.0 - (pmiss[:, 2:4] * ones); // matrix trick to get row sums complete cases
    for(i in 1:N){
      target += log(pmiss[i, R[i]]);
    }
  } // end local
}
generated quantities{
  // use this for rejection sampling (stan operates on log probabilities and cannot handle a discrete 1/0 increment to the posterior)
  int<lower=0,upper=1> constraint = 0;
  { // local
  real pmissi[NR];
  for(i in 1:N){
   // this could be more efficient by looping only over those with R[i]==1
     pmissi[2] = inv_logit(g0[1] + gamma[1]*L[i,1] + gamma[2]*L[i,2]                   );
     pmissi[3] = inv_logit(g0[2] +                   gamma[3]*L[i,2]                   );
     pmissi[4] = inv_logit(g0[3] +                   gamma[4]*L[i,2] + gamma[5]*L[i,3] );
    if ( R[i] == 1 && sum(pmissi[2:4]) >= 1.0-sigmastar){
       constraint = 1;
       break;
     }
  }
  } // end local
}
"


# initial values (one set for each chain)
#crucial to make this consistent with a valid first estimate of pi
stanchains=4
initlist = lapply(1:stanchains,
  function(x) list(g0=runif(3, min=-7,max=-3), gamma=runif(5, min=-.05,max=.05))
) 

#compile model (ingore all warnings)
testfit = stan(model_code = stanmod, data=dat, chains=1, iter=1, init = list(initlist[[1]]))
# fit model
ft2 = stan(fit=testfit, data=dat, chains=stanchains, iter=20000, warmup=1000, init=initlist)
ft2
stan_trace(ft2, pars='g0')
stan_trace(ft2, pars='gamma')
stan_dens(ft2, pars='g0')
stan_dens(ft2, pars='gamma')

# rejection sampling
allsamples = as.data.frame(ft2)
allsamples$divergent = as.numeric(rstan::get_divergent_iterations(ft2))
#with(allsamples, prop.table(table(divergent, constraint), margin = 2))

constrained_samples <- allsamples[allsamples$constraint==0,]
dim(allsamples)
dim(constrained_samples)

ggplot()+ theme_classic() + scale_color_discrete(name="")+
  geom_density(aes(x=`lp__`, color="Unconstrained"), data = allsamples)+
  geom_density(aes(x=`lp__`, color="Constrained"), data = constrained_samples)

ggplot()+ theme_classic() + scale_color_discrete(name="")+
  geom_density(aes(x=`g0[1]`, color="Unconstrained"), data = allsamples)+
  geom_density(aes(x=`g0[1]`, color="Constrained"), data = constrained_samples)

ggplot()+ theme_classic() + scale_color_discrete(name="")+
  geom_density(aes(x=`gamma[4]`, color="Unconstrained"), data = allsamples)+
  geom_density(aes(x=`gamma[4]`, color="Constrained"), data = constrained_samples)
