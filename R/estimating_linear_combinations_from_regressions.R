# getting estimates for linear combinations of variables in a regression
# e.g. contrasting doubly exposed with unexposed

est_lincom <- function (terms, coefs, grad = NULL) {
  weightvec <- rep(0, length(coefs))
  if (is.null(grad)){
    weightvec[which(names(coefs) %in% terms)] <- 1
  } else if (!is.null(grad) && length(grad)==1){
    weightvec[which(names(coefs) %in% terms)] <- grad
  } else if (!is.null(grad) && length(grad)==length(weightvec)){
    weightvec <- grad
  }
  eff <- weightvec %*% coefs
  eff[1]
}

se_lincom <- function (terms, covmat, grad = NULL) {
  if (!is.matrix(covmat)) {
    nm <- names(covmat)
    covmat = matrix(covmat)
    colnames(covmat) <- nm
  }
  weightvec <- rep(0, dim(covmat)[1])
  if (is.null(grad)){
    weightvec[which(colnames(as.matrix(covmat)) %in% terms)] <- 1
  } else if (!is.null(grad) && length(grad)==1){
    weightvec[which(colnames(as.matrix(covmat)) %in% terms)] <- grad
  } else if (!is.null(grad) && length(grad)==length(weightvec)){
    weightvec <- grad
  }
  var <- weightvec %*% covmat %*% weightvec
  sqrt(var)[1, , drop = TRUE]
}

# simulate data for testing
n=100
z <- rnorm(n)
x <- rnorm(n) + z
y <- rbinom(n, 1, exp(-2 + x + z*x))

fit <- glm(y~x+z*x, family=binomial(link="log"), data = data.frame(y,x,z), start=c(-1,0,0,0))
coef = coef(fit)
vcov = vcov(fit)
summary(fit) # for comparison

coef["x"]
se_lincom(c("x"), covmat=vcov)         # standard error of beta_x, the hard way


# estimate E(y | x = 1, z=1)
# beta_0 + beta_x + beta_z + beta_x*z
est_lincom(c("(Intercept)", "x", "z", "x:z"), coefs =coef)   # linear combination of beta_0 + beta_x + beta_z + beta_x*z
se_lincom(c("(Intercept)", "x", "z", "x:z"), covmat=vcov)    # standard error of beta_0 + beta_x + beta_z + beta_x*z
# equivalent (ignores terms argument)
est_lincom(c("(Intercept)", "x", "z", "x:z"), coefs =coef, grad = c(1,1,1,1))
se_lincom(c("(Intercept)", "x", "z", "x:z"), covmat=vcov, grad = c(1,1,1,1))
# also equivalent
est_lincom(coefs = coef, grad = c(1,1,1,1))
se_lincom(covmat = vcov, grad = c(1,1,1,1)) 

# contrast E(y | x = 1, z=1) with E(y | x=0, z=3)
# = beta_0 + beta_x + beta_z + beta_x*z - (beta_0 + 3*beta_z)
# = -beta_x - 2*beta_z 
contr = c(0,-1,-2,0)
est_lincom(coefs =coef, grad = contr)
se_lincom(covmat=vcov,  grad = contr)  

