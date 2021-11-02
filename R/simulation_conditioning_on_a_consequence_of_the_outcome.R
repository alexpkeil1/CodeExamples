#shows how bias occurs from conditioning on a collider that is a consequence of the outcome, even if you
# block the path

dgm <- function(N, beta=0){
 x <- rnorm(N, 0, 1)
 m <- rnorm(N, x, 1)
 y <- rnorm(N, beta*x, 1)
 c <- rnorm(N, m + y, 1)
 data.frame(x,m,c,y)  
}

analyze <- function(dat){
 c <- lm(y ~ x, data=dat)$coefficients[2]
 ma <- lm(y ~ x + m, data=dat)$coefficients[2]
 ca <- lm(y ~ x + c, data=dat)$coefficients[2]
 mca <- lm(y ~ x + c + m, data=dat)$coefficients[2]
 c(c, ma, ca, mca)
}

popsize=200
iterate <- function(niter, ...){
 res <- matrix(nrow=niter, ncol=4)
 for(i in 1:niter){
  dat <- dgm(popsize, ...)
  res[i,] <- analyze(dat)
 }
 colnames(res) <- c("crude", "m-adjusted", "collider adjusted", "m and collider adjusted")
 res
}




#bias occurs in the M and collider adjusted model only when beta is not equal to 0
B=3

res <- iterate(1000, beta=B)
options(scipen=10)

apply(res,2, function(x) c(mean=mean(x), 
                           bias=mean(x)-B, 
                           var=var(x), 
                           T_nobias=(T <- (mean(x)-B)/(sqrt(var(x)/popsize))), 
                           p_nobias = (1-pnorm(abs(T)))*2))


