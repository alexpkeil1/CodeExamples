dgm <- function(N, setX=NULL, seed=NULL){
 if(!is.null(seed)) set.seed(seed)
 u1 <- rbinom(N, 1, 0.5) 
 u2 <- rbinom(N, 1, 0.5) 
 z  <- rbinom(N, 1, plogis(-(-2 + 2*u1 + 2*u2))) 
 x  <- rbinom(N, 1, plogis(-(-1 + 2*u1)))
 if(!is.null(setX)) x=setX
 y  <- rbinom(N, 1, plogis(-(-1.5 + 2*u2 + x))) 
 data.frame(x,y,z)
}

# pretending it's a confounder
gformula_zconfounder<- function(dat){
 py11 = with(dat, mean(y[x==1 & z==1])) 
 py01 = with(dat, mean(y[x==0 & z==1])) 
 py10 = with(dat, mean(y[x==1 & z==0])) 
 py00 = with(dat, mean(y[x==0 & z==0])) 
 pz = with(dat, mean(z))
 
 y1 = py11*pz + py10*(1-pz)
 y0 = py01*pz + py00*(1-pz)
 y1-y0
 
}

# pretending it's a mediator
gformula_zmediator <- function(dat){
 py11 = with(dat, mean(y[x==1 & z==1])) 
 py01 = with(dat, mean(y[x==0 & z==1])) 
 py10 = with(dat, mean(y[x==1 & z==0])) 
 py00 = with(dat, mean(y[x==0 & z==0])) 
 pz1 = with(dat, mean(z[x==1]))
 pz0 = with(dat, mean(z[x==0]))
 
 y1 = py11*pz1 + py10*(1-pz1)
 y0 = py01*pz0 + py00*(1-pz0)
 y1-y0
}

gformula_xrandomized <- function(dat){
 y1 = with(dat, mean(y[x==1])) 
 y0 = with(dat, mean(y[x==0])) 
 y1-y0
}

seed = 12321
#truth ATE,  -0.1952
round(mean(dgm(100000000, setX=1, seed)$y) - mean(dgm(100000000, setX=0, seed)$y), 4)

dat = dgm(100000000, seed=seed)
#if Z treated like a confounder, -0.2256 (biased for ATE)
round(gformula_zconfounder(dat), 4)

#if Z treated like a mediator, -0.1951 (unbiased for ATE)
round(gformula_zmediator(dat), 4)

#if Z treated like a collider (x randomized) -0.1951 (unbiased for ATE)
round(gformula_xrandomized(dat), 4)
