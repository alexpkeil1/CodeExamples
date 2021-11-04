# transporting causal effects from randomized trials and non-randomized studies
# with identified effects

#https://ftp.cs.ucla.edu/pub/stat_ser/r400-reprint.pdf


######################################## 
# example 1: direct transport unconditionally
########################################

sem1 <- function(n){
  # study sample
  # error terms
  ux <- rnorm(n)
  uz1 <- rnorm(n)
  uu <- rnorm(n)
  uy <- rnorm(n)
  # structural model
  fz1 <- function(uz1) uz1
  fu <- function(uu) uu
  #fx <- function(z1,u,ux) z1 + u + ux
  fx <- function(z1,u,ux) ux    # randomized
  fy <- function(x,z1,u,uy) x + z1 + u + x*z1 + uy
  # realizations
  z1 <- fz1(uz1)
  u <- fu(uu)
  x <- fx(z1,u,ux)
  y <- fy(x,z1,u,uy)
  y0 <- fy(0,z1,u,uy)
  y1 <- fy(1,z1,u,uy)
  data.frame(z1,u,x,y,y0,y1)
}

sem2 <- function(n){
  # target population
  # error terms
  ux <- rnorm(n)
  uz1 <- rnorm(n)
  uu <- rnorm(n)
  uy <- rnorm(n)
  # structural model
  fz1 <- function(uz1) uz1
  fu <- function(uu) uu
  fx <- function(z1,u,ux) z1 + u + ux # non-randomized
  #fx <- function(z1,u,ux) ux    # randomized
  fy <- function(x,z1,u,uy) x + z1 + u + x*z1 + uy
  # realizations
  z1 <- fz1(uz1)
  u <- fu(uu)
  x <- fx(z1,u,ux)
  y <- fy(x,z1,u,uy)
  y0 <- fy(0,z1,u,uy)
  y1 <- fy(1,z1,u,uy)
  data.frame(z1,u,x,y,y0,y1)
}

# notation:
# sample domain is Pi, target domain is Pi*
# <M,M*> structural causal models for sample and target (scm1, scm2)
# G=G* is the causal diagram shared by both SCMs with nodes V_1, V_2, ..., V_i
# If G not shared then start with G*
# D (selection diagram) starts with G* and then adds S_i -> V_i for each 
#  V_i in D where f_i /neq f*_i or P(U_i) \neq P*(U_i)

# transporting 
sample <- sem1(1000)
target <- sem2(1000)

# transporting of causal effect R(Pi*) = P*(y|do(x),z) is transportable from 
#  sample (Pi) to target (Pi*) in D  if R(Pi*) is uniquely computable from P,P*,I
# in any model that induces D

# fx varies across domains, but causal DAGs (G) / selection diagrams (D) would be
# G       z1,u -> y; x -> y;              # identified (randomized)
# G* x <- z1,u -> y; x -> y;              # non identified (u unmeasured)
# D  x <- z1,u -> y; x -> y; s -> x

#P(y | do(x), z) is transportable from G to G* if z,x d-separates y from s in D
# trivially, this is true here because P(y | do(x)) = p(y | do(x), s)
# so the marginal and z-specific effect of x on y is identical in the sample and the target
(samplemod <- lm(y~x, data = sample)) # directly applicable (and sample is a trial)
# using simple model, transporting P(y|do(x)=1)-P(y|do(x)=0)
preds1triv <- predict(samplemod, newdata = data.frame(x=1, target[,c("z1"), drop=FALSE])) 
preds0triv <- predict(samplemod, newdata = data.frame(x=0, target[,c("z1"), drop=FALSE])) 
mean(preds1triv) - mean(preds0triv) # trivially transporting single causal effect (valid)

(targetmod  <- lm(y~x, data = target)) # confounded (if u is measured it's trivial transport)
# using flexible model, transporting P(y|do(x)=1)-P(y|do(x)=0)
flexmod <- lm(y~x*z1, data = sample)
preds1 <- predict(flexmod, newdata = data.frame(x=1, target[,c("z1"), drop=FALSE])) 
preds0 <- predict(flexmod, newdata = data.frame(x=0, target[,c("z1"), drop=FALSE])) 
mean(preds1) - mean(preds0)  # directly transporting single causal effect (valid)
mean(target$y1 - target$y0) # estimate from potential outcomes under data generation


######################################## 
# example 1: direct transport (on z1)
########################################

sem2b <- function(n){
  # target population
  # error terms
  ux <- rnorm(n)
  uz1 <- rnorm(n)
  uu <- rnorm(n)
  uy <- rnorm(n)
  # structural model
  fz1 <- function(uz1) exp(uz1)
  fu <- function(uu) uu
  fx <- function(z1,u,ux) z1 + u + ux # non-randomized
  #fx <- function(z1,u,ux) ux    # randomized
  fy <- function(x,z1,u,uy) x + z1 + u + x*z1 + uy
  # realizations
  z1 <- fz1(uz1)
  u <- fu(uu)
  x <- fx(z1,u,ux)
  y <- fy(x,z1,u,uy)
  y0 <- fy(0,z1,u,uy)
  y1 <- fy(1,z1,u,uy)
  data.frame(z1,u,x,y,y0,y1)
}

# fx,fz vary across domains, but causal DAGs (G) / selection diagrams (D) would be
# G       z1,u -> y; x -> y;              # identified (randomized)
# G* x <- z1,u -> y; x -> y;              # non identified
# D  x <- z1,u -> y; x -> y; s -> x,z1    # added s -> z1 because marginal z1 distribution varies over domains


target_b <- sem2b(1000)

# attempting transport unconditional on z
(samplemod <- lm(y~x, data = sample)) # directly applicable (and sample is a trial)
# using simple model, transporting P(y|do(x)=1)-P(y|do(x)=0)
preds1triv <- predict(samplemod, newdata = data.frame(x=1, target_b[,c("z1"), drop=FALSE])) 
preds0triv <- predict(samplemod, newdata = data.frame(x=0, target_b[,c("z1"), drop=FALSE])) 
mean(preds1triv) - mean(preds0triv) # trivially transporting single causal effect (not valid)
mean(target_b$y1 - target_b$y0) # estimate from potential outcomes under data generation

# using direct transport conditional on z
(targetmod  <- lm(y~x, data = target_b)) # confounded
# using flexible model, transporting P(y|do(x)=1,z1)-P(y|do(x)=0,z1)
flexmod <- lm(y~x*z1, data = sample)
preds1 <- predict(flexmod, newdata = data.frame(x=1, target_b[,c("z1"), drop=FALSE])) 
preds0 <- predict(flexmod, newdata = data.frame(x=0, target_b[,c("z1"), drop=FALSE])) 
mean(preds1) - mean(preds0)  # directly transporting single causal effect (valid)
mean(target_b$y1 - target_b$y0) # estimate from potential outcomes under data generation

######################################## 
# example 3: direct transport (on z1) from non-randomized trial
########################################
sem1b <- function(n){
  # study sample
  # error terms
  ux <- rnorm(n)
  uz1 <- rnorm(n)
  uu <- rnorm(n)
  uy <- rnorm(n)
  # structural model
  fz1 <- function(uz1) uz1
  fu <- function(uu) uu
  fx <- function(z1,u,ux) 4*z1 + ux # non-randomized, but identifiable 
  #fx <- function(z1,u,ux) ux    # randomized
  fy <- function(x,z1,u,uy) x + z1 + u + x*z1 + uy
  # realizations
  z1 <- fz1(uz1)
  u <- fu(uu)
  x <- fx(z1,u,ux)
  y <- fy(x,z1,u,uy)
  y0 <- fy(0,z1,u,uy)
  y1 <- fy(1,z1,u,uy)
  data.frame(z1,u,x,y,y0,y1)
}

# fx varies across domains, but causal DAGs (G) / selection diagrams (D) would be
# G  x <- z1   -> y; x -> y;              # identified (z1 is sufficient)
# G* x <- z1,u -> y; x -> y;              # non identified
# D  x <- z1,u -> y; x -> y; s -> x,z1


sample_b <- sem1(1000)
target_b <- sem2b(1000)

# attempting transport unconditional on z
(samplemod <- lm(y~x*z1, data = sample_b)) 
# using simple model, transporting P(y|do(x)=1,z1)-P(y|do(x)=0,z1)
preds1triv <- predict(samplemod, newdata = data.frame(x=1, target_b[,c("z1"), drop=FALSE])) 
preds0triv <- predict(samplemod, newdata = data.frame(x=0, target_b[,c("z1"), drop=FALSE])) 
mean(preds1triv) - mean(preds0triv) # trivially transporting single causal effect (not valid)
mean(target_b$y1 - target_b$y0) # estimate from potential outcomes under data generation


######################################## 
# example 4: trivial transport (identifiability in target - we don't even need study sample)
########################################

target_b$z2 = target_b$u

# G* x <- z1,z2 -> y; x -> y;              # identified (z2, z1 sufficient)
# D  x <- z1,z2 -> y; x -> y; s -> x

# trivial transport conditional on z1, u (i.e. we meet identifiability in the target)
(targetmod  <- lm(y~x, data = target_b)) # confounded
(targetmod  <- lm(y~x+z1, data = target_b)) # still confounded
# using flexible model, transporting P(y|do(x)=1)-P(y|do(x)=0)
flexmod <- lm(y~x*z1*z2, data = target_b)
preds1 <- predict(flexmod, newdata = data.frame(x=1, target_b[,c("z1", "z2"), drop=FALSE])) 
preds0 <- predict(flexmod, newdata = data.frame(x=0, target_b[,c("z1", "z2"), drop=FALSE])) 
mean(preds1) - mean(preds0)  # directly transporting single causal effect (valid)
mean(target_b$y1 - target_b$y0) # estimate from potential outcomes under data generation
