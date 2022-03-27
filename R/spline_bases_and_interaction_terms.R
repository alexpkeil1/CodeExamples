# restricted cubic spline
rcs <- function(var, k=NULL, nk=NULL){
  #restricted cubic spline
  if(!is.null(nk)){
    #knot placement taken from sas %DASPLINE macro
    if( nk==3 ) p = c( 5, 50, 95)
    if( nk==4 ) p = c( 5, 35, 65, 95)
    if( nk==5 ) p = c( 5, 27.5, 50, 72.5, 95)
    if( nk==6 ) p = c( 5, 23, 41, 59, 77, 95)
    if( nk==7 ) p = c( 2.5, 18.3333, 34.1667, 50, 65.8333, 81.6667, 97.5)
    if( nk==8 ) p = c( 1, 15, 29, 43, 57, 71, 85, 99)
    if( nk==9 ) p = c( 1, 13.25, 25.5, 37.75, 50, 62.25, 74.5, 86.75, 99)
    if( nk==10 ) p = c( 1, 11.88889, 22.77778, 33.66667, 44.55556, 55.44444, 66.33333, 77.22222, 88.11111, 99)
    if( nk==11 ) p = c( 1.0, 10.8, 20.6, 30.4, 40.2, 50.0, 59.8, 69.6, 79.4, 89.2, 99.0)
    k = quantile(x=var, p=p/100)
    cat("\nknots:", k, "\n")
  }
  numKnots <- length(k)
  numVar <- numKnots-2
  if(numVar<1) stop("number of knots must be at least 3")
  mat <- matrix(ncol=numVar, nrow=length(var))
  for(i in 1:numVar){
    # following "norm=2" from frank harrel's daspline
    mat[,i] <- ((var>k[i]) * (var-k[i])/((k[numKnots] - k[1])**(2/3)))**3 +
      ((k[numKnots-1]-k[i])* ((var>k[numKnots]) * (var-k[numKnots])/((k[numKnots] - k[1])**(2/3)))**3 -
         (k[numKnots]-k[i]) * ((var>k[numKnots-1]) * (var-k[numKnots-1])/((k[numKnots] - k[1])**(2/3)))**3)/(k[numKnots]-k[numKnots-1])
  }
  colnames(mat) <- paste0(deparse(substitute(var)), "_sp", 1:numVar)
  mat
}

n = 1000
# country: France or US
fr = rbinom(n, 1, 0.5)
us = 1-fr
# continuous exposure with 3 sd difference in mean b/w countries
x = rnorm(n,3*fr)
# country specific knots/spline basis calculation
xsp1 = rcs(x[fr==1], nk=4)
xsp2 = rcs(x[us==1], nk=4)
# combined knots
xspa = rcs(x, nk=4)

# combine spline basis terms into single matrix
xsp = matrix(NA, nrow=n, ncol=ncol(xsp1))
xsp[which(fr==1),] = xsp1
xsp[which(us==1),] = xsp2

# generate outcomes with true coefficients all equal to 1.0 or -1.0, country specific variables underly true model
mu = us*x + us*xsp%*%c(1,1)  +
     fr*x + (fr*xsp)%*%c(-1,-1) +
     fr
y = rnorm(n, mu, 1)


# generative model
# this returns unbiased estimates - using only country*spline interactions rather than spline + country*spline
f1 <- lm(y~ fr +
          fr:x+fr:xsp+
          us:x+us:xsp
          )

# this is an equivalent model fit to generative model, but parameterizes the model differently (spline + country*spline)
f2 <- lm(y~fr +
             us:x+us:xsp +
             x+xsp
          )

# this is a different model, since spline variables are created using combined knots
# at face value, it looks like the generative model
f3 <- lm(y~fr +
           fr:x+fr:xspa+
           us:x+us:xspa
)


summary(f1)
summary(f2)
summary(f3)
