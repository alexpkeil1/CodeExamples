# direct effects using logistic regression under assumption of no interaction on 
# additive scale

dgm <- function(N, setX=NULL, setXkeepz=NULL, setZ=NULL, seed=NULL){
  if(is.null(seed)) seed = sample(1:10e6, 1)
  set.seed(seed)
  x <- rbinom(N, 1, 0.5)
  if(!is.null(setX)) x <- rep(setX, times=N)
  z <- rbinom(N, 1, 1/(1+exp(-0.5-x)))
  if(!is.null(setZ)) z <- rep(setZ, times=N)
  if(!is.null(setXkeepz)) x <- rep(setXkeepz, times=N)
  # the toggle - is there XY modification by Z on the additive effect scale?
  #y <- rbinom(N, 1, 0.4 + .2*x + .2*z)
  y <- rbinom(N, 1, 0.4 + .2*x + .2*z - .2*z*x)
  data.frame(x,y,z)
}

analyze <- function(SEED){
  N=500
  SEED=as.numeric(SEED)
  dat <- dgm(N, seed=SEED)
  # setting Z to the observed value
  dat1 <- dgm(N, seed=SEED, setXkeepz=1)
  dat0 <- dgm(N, seed=SEED, setXkeepz=0)
  truecde <- mean(dat1$y)-mean(dat0$y)
  mod <- glm(y ~ x + z + x*z, data=dat, family=binomial())
  r1 <- predict(mod, newdata=dat1, type='response')
  r0 <- predict(mod, newdata=dat0, type='response')
  estcde = mean(r1-r0)

  # setting Z to a constant value
  dat1 <- dgm(N, seed=SEED, setX=1, setZ=1)
  dat0 <- dgm(N, seed=SEED, setX=0, setZ=1)
  truecde2 <- mean(dat1$y)-mean(dat0$y)
  mod <- glm(y ~ x + z + x*z, data=dat, family=binomial())
  r1 <- predict(mod, newdata=dat1, type='response')
  r0 <- predict(mod, newdata=dat0, type='response')
  estcde2 = mean(r1-r0)

  # setting only X
  dat1 <- dgm(N, seed=SEED, setX=1)
  dat0 <- dgm(N, seed=SEED, setX=0)
  trueme <- mean(dat1$y)-mean(dat0$y)
  mod <- glm(y ~ x, data=dat, family=binomial())
  r1 <- predict(mod, newdata=dat1, type='response')
  r0 <- predict(mod, newdata=dat0, type='response')
  estme = mean(r1-r0)
  c(truecde,truecde2,trueme,estcde,estcde2,estme)
}

system.time(results <- sapply(1:1000, analyze))

pr <- apply(results, 1, mean)
names(pr) <- c('truecde','truecde2', 'trueme','estcde','estcde2', 'estme') 
print(pr, 4)
