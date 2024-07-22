
# simple simulation + multivariate normal model in Stan
# also see: seemingly unrelated regression (here, an example using estimating equations)
# demonstrates: simple model fitting and troubleshooting
library(mvtnorm)
library(future)
library(future.apply)
library(rstan)
library(mvtnorm)
library(rootSolve)
library(numDeriv)

### simulate data from a multivariate model
dgm <- function(n,b1=c(1,1,0), b2=c(0,0,1)){
  X1 <- rmvnorm(n, rep(0,3), diag(c(1,2,3)))
  Y1 <- X1 %*% b1 + rnorm(n)
  X2 <- X1
  Y2 <- X1 %*% b2 + rnorm(n)
  data.frame(x1=X1,y1=Y1, x2=X2, y2=Y2)
}




linreg_ee <- function(Y,X,p,theta){
  if(length(Y)==0) return(rep(0,length(theta)))
  mu = X %*% theta
  do.call(c,lapply(1:p, function(x){
    t((Y - mu)) %*% X[,x]
  } ))
}

normreg_ee <- function(Y,X,ppluss,zeta){
  p = ppluss - 1
  theta = zeta[1:p]
  sigma2 = zeta[p+1]
  if(length(Y)==0) return(rep(0,length(theta)+1))
  mu = X %*% theta
  fun = do.call(c,lapply(1:p, function(x){
    t((Y - mu)) %*% X[,x]
  } ))
  fun = c(fun,
          log(mean((Y - X %*% theta)^2)) - sigma2
  )
}

sur_esteq <- function(theta, Ylist, Xlist){
  # linear model gradient-based estimating equations
  dim1 = dim(Xlist[[1]])
  dim2 = dim(Xlist[[2]])
  fun1 = linreg_ee(Ylist[[1]], Xlist[[1]], dim1[2], theta[1:dim1[2]])
  fun2 = linreg_ee(Ylist[[2]], Xlist[[2]], dim2[2], theta[(dim1[2]+1):(dim1[2]+dim2[2])])
  c(
    fun1,fun2
  )
}

sur_esteq_df <- function(theta, df,y1col,y2col,x1cols,x2cols){
  sur_esteq(theta,list(df[,y1col,drop=TRUE],df[,y2col,drop=TRUE]),list(as.matrix(df[,x1cols,drop=FALSE]),as.matrix(df[,x2cols,drop=FALSE])))
}

datl = dgm(100)
y1col = which(names(datl)=="y1")
y2col = which(names(datl)=="y2")
x1col = grep("x1", names(datl))
x2col = grep("x2", names(datl))

theta = rep(0, length(x1col)+length(x2col))

eqfit <- rootSolve::multiroot(sur_esteq_df, start=theta, df=datl, y1col=y1col,y2col=y2col, x1col=x1col,x2col=x2col)
#
A = numDeriv::jacobian(func=sur_esteq_df, x=eqfit$root, df=datl, y1col=y1col,y2col=y2col, x1col=x1col,x2col=x2col)


uid =   1:nrow(datl)
psii = lapply(uid, function(x){
  selidx=x
  sur_esteq_df(eqfit$root, datl[x,,drop=FALSE],y1col,y2col,x1col,x2col) 
} )
Bi = lapply(psii, function(x) x%*%t(x))
n = length(uid)
B = Bi[[1]]
for(i in 2:length(Bi)){
  B = B + Bi[[i]]
}

ibread = solve(A)
(fullcovmat = ibread %*% B %*% t(ibread))


############################################################
# Bayesian multivariate model using Stan
############################################################
  
# rstan functions use a list as input for data and constants
standata = list() # data
standata$N = nrow(datl) # constants
standata$X = datl[,x1col] # just ignore that there are two predictor matrixes (which are identical here)
standata$Y = datl[,c(y1col,y2col)] # multivariate outcome
standata$J = ncol(standata$X)
standata$K = ncol(standata$Y)

##############
#  stan code for multivariate linear model
##############

stanmod <- "
data {
  int<lower=0> N;          // sample size
  int<lower=0> K;          // number of outcomes
  int<lower=0> J;          // design matrix columns
  matrix[N,K] Y;               // estimated treatment effects
  matrix [N,J] X;               // estimated treatment effects
}
parameters {
  cov_matrix[K] Sigma; // ensures positive semi-definite
  matrix[K,J] beta;
}
model {
  //Sigma ~ inv_wishart(2.0, [[2, 0], [0, 2]]); // use vague prior by commenting this line
  for(k in 1:K){
    beta[k,1] ~ normal(0,1000); // vague prior on intercepts
    for(j in 2:J){
      beta[k,j] ~ normal(0,1); // shrinkage prior
    }
  }
  matrix [N,K] mu;
  for(i in 1:N){
    for(k in 1:K){
      mu[i,k] =  X[i,] * beta[k,]';                 // beta is structured this way just so that estimates are printed in an order that matches estimating equation approach
    }
    target += multi_normal_lpdf(Y[i,] | mu[i,], Sigma); 
  }
}
generated quantities{
}
"

# compile model
compiledfit <- stan(model_code=stanmod, data = standata, chains=1, iter=10)

# run more samples, check for convergence
postsamples = stan(fit=compiledfit, iter=10000, data = standata)

# all parameter estimates
print(postsamples)   # Bayesian fit


# model coefficient parameter estimates
print(postsamples, pars="beta")   # Bayesian fit
cbind(est=eqfit$root, stderr = sqrt(diag(fullcovmat))) # estimating equations approach


# covariance matrix for coefficients
allsamps = as.matrix(postsamples)
betapost = allsamps[,c("beta[1,1]","beta[1,2]","beta[1,3]","beta[2,1]","beta[2,2]","beta[2,3]")]

round(100*cov(betapost), 3)
round(100*fullcovmat,3) # estimating equations approach
