# two issues in modified poisson modeling to estimate a risk ratio

# 1: estimates from modified poisson may not maximize the binomial likelihood 
# 2: estimates from modified poisson may imply a binomial likelihood of zero due to in-sample predictions > 1.0


loglik_bin <- function(beta, X, Y){
  D =  cbind(rep(1, length(X)), X)
  lmu = D%*%beta # log probability y = 1
  pred = signif(exp(lmu), digits=8) # account for machine error
  if(any(((pred > 1) + (pred < 0)))){
    message("Predicted probabilities outside valid range")
  } else   if(any(((pred==1) + (pred == 0))))
    message("Estimates on edge of parameter space")
  lli = suppressWarnings(ifelse(Y==1, lmu , log(1-exp(lmu))))
  ll = sum(lli)
  ll
}




# example 1: MLE on edge of parameter space, estimates are similar for log binomial and modified poisson
y = c(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0)
x = c(0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, 1.0)
mean(y)

(m1 <- glm(y~x, family=binomial(link = "log"), start=c(-2,0)))  # log binomial
(m2 <- glm(y~x, family=poisson(link = "log")))    # modified poisson (point estimates don't rely on robust variance)

loglik_bin(coef(m1), x, y) # MLE is at the boundary of the parameter space
loglik_bin(coef(m2), x, y) # note that this estimate isn't the MLE for a binomial outcome!



summary(predict(m1, type = "response"))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1874  0.1874  0.1874  0.2978  0.2848  1.0000
# results in predicted probabilities for the observed data <= 1

# modified poisson (valid probabilities)
summary(predict(m2, type = "response"))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1172  0.1172  0.1172  0.3000  0.1978  0.9505
# results in predicted probabilities for the observed data < 1

# example 2: parameter boundary issues result in problems for log binomial model
# from https://academic-oup-com.libproxy.lib.unc.edu/aje/article/159/2/213/166397?login=true
datstr <- textConnection("
x   y
1 	0 
2 	0 
3 	0 
4 	0 
5 	1 
6 	0 
7 	1 
8 	1 
9 	1 
10 	1 
")
dat = read.table(datstr, header=TRUE)
mean(dat$y)

# true MLE: intercept  = -2.09, slope = .2094
            
(m1b <- glm(y~x, data = dat, family=binomial(link = "log"), start=c(-1,0)))  # log binomial
warnings() # should see a note that algorithm stopped at boundary value
(m2b <- glm(y~x, data = dat, family=poisson(link = "log"), start=c(0,0)))    # modified poisson (point estimates don't rely on robust variance)
sascoef = c(-2.8634, 0.2863) # with same starting values in PROC GENMOD, SAS gives up and yields these estimates with non-convergence message
sascoef = c(-3.6373, 0.3637) # with same starting values in PROC GENMOD, SAS gives up and yields these estimates with non-convergence message

loglik_bin(coef(m1b), dat$x, dat$y)          # likelihood at estimated MLE - converged at parameter space boundary
loglik_bin(coef(m2b), dat$x, dat$y)          # higher apparent likelihood for Poisson estimate, but log-likelihood is -Inf because this implies probabilities > 1.0 in the data
loglik_bin(sascoef, dat$x, dat$y)            # not the MLE, implying non-convergence


summary(predict(m1b, type = "response"))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.1520  0.2445  0.3920  0.4642  0.6269  1.0000 
# results in predicted probabilities for the observed data <= 1

# modified poisson
summary(predict(m2b, type = "response"))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0774  0.1625  0.3387  0.5000  0.7012  1.4437 
# results in predicted probabilities for the observed data > 1 (meaning the likelihood should be 0)
