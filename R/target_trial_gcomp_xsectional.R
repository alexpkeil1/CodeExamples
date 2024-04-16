# an example of using g-computation with a single index exposure
# contrasting two hypothetical scenarios

# step 0) 
# bring in some data
data(metals, package="qgcomp")

mixture = c("cadmium", "lead", "chromium", "mercury", "arsenic")

summary(metals[,mixture]) # these are already centered
apply(metals[,mixture], 2, sd) # these are already scaled

# step 1) 
# Perform a data reduction to get a single index exposure.
# Many methods possible: here I use PCA and take the first principal component

pcafit = prcomp(metals[,mixture])
metals$idxexposure = predict(pcafit)[,1]

# step 2) 
# now fit our fully adusted model (only one confounder here)

gcompfit = lm(y ~ idxexposure + mage35, data = metals)


# step 3) 
# Make predictions for the target population. There are several options here.

# 3a) the "natural course"
py_nc = predict(gcompfit)

# 3b) a hypothetical intervention to decrease exposures by 60%
#  i.e. X_int = X * (1-pct_decrease)
# Example: if an individual has a level of 0.9 ug/L for cadmium,
#  we would decrease their exposure to 0.9 *(1-0.6) = 0.36 ug/L.
#  Other exposures could be reduced by the same, or different amounts.
# Tricky part: the data are already standardized
#   For a z-score transform like we've done, we have
#   X_z = (X-mn)/std                  # standardized exposure
#   X_int = X * (1-pct_decrease)      # post-intervention exposure
#   X_intz = (X_int-mn)/std           # post-intervention exposure after standardizing
# Then:
#   X_intz = (X_int-mn)/std
#   X_intz = (X * (1-pct_decrease)-mn)/std
#   X_intz = (X-X*pct_decrease -mn)/std
#   X_intz = (X -mn)/std -(X*pct_decrease)/std
#   X_intz = X_z - X*pct_decrease/std
#   X_intz = X_z - (X_z*std+mn)*pct_decrease/std         
#   X_intz = X_z - (X_z+mn/std)*pct_decrease
#   X_intz = X_z - X_z*pct_decrease - mn/std*pct_decrease
#   X_intz = X_z*(1 - pct_decrease) - mn/std*pct_decrease    # gives you the appropriate post-intervention value if you know the original mean & std
#   

#   For a log transform done, we have
#   X_l = log(X)                  # standardized exposure
#   X_int = X * (1-pct_decrease)      # post-intervention exposure
#   X_intl = log(X_int)
# Then:
#   X_intl = log(X_int)
#   X_intl = log(X * (1-pct_decrease))
#   X_intl = log(X) + log(1-pct_decrease)
#   X_intl = X_l + log(1-pct_decrease)


# quick example with standardized exposures
x = rexp(100)
mn = mean(x)
std = sd(x)
x_z = (x-mn)/std
pct_decrease = .8

x_int = x * (1-pct_decrease)
(x_intz = (x_int - mn)/std)
x_z*(1 - pct_decrease) - mn/std*pct_decrease






# back to original example, 60% reduction
pct_decrease = 0.6
# 3b.1) create a copy of the original, standardized data under the intervention

# original means:
metalmeans = c(cadmium=0.001, lead=0.03, chromium=0.1, mercury=0.01, arsenic=0.1)
# original standard deviations
metalstds = c(cadmium=0.1, lead=0.3, chromium=0.1, mercury=0.2, arsenic=0.4)

metals60 = metals # make a copy of the data
metals60[,mixture] = metals[,mixture] * (1-pct_decrease) - metalmeans/metalstds * pct_decrease

# 3b.2) predict the index variable given the post-intervention exposures

metals60$idxexposure = predict(pcafit, newdata = metals60[,mixture])[,1]

# 3b.3) make predictions under the interventions data
py_int60 = predict(gcompfit, newdata = metals60)

mean(py_int60- py_nc) # mean difference: effect of joint reduction in metals that go into the index is to increase the average outcome


# 3b.3alt) alternative way to do g-computation here is to fit a saturated model to the combined data
metals$intervention = 0
metals60$intervention = 1
metals$py = py_nc
metals60$py = py_int60
alldata = rbind(metals, metals60)
lm(py~intervention, data = alldata)






# 4) Bootstrapping for variance
# Note, if we just try to get CI by doing this: confint(lm(py~intervention, data = alldata))
# We get invalid CI because individuals are in the data multiple times.
# Bootstrapping is just doing the whole analysis on multiple draws from the original dataset

# 4a) In R, it helps to create a function that just does the entire analysis - here we'll just modify the code from above

gcomp_pca <- function(data, pct_decrease=0.6){
  mixture = c("cadmium", "lead", "chromium", "mercury", "arsenic")
  # original means:
  metalmeans = c(cadmium=0.001, lead=0.03, chromium=0.1, mercury=0.01, arsenic=0.1)
  # original standard deviations
  metalstds = c(cadmium=0.1, lead=0.3, chromium=0.1, mercury=0.2, arsenic=0.4)
  # step 1)  
  pcafit_boot = prcomp(data[,mixture])
  data$idxexposure = predict(pcafit)[,1]
  # step 2) 
  gcompfit_boot = lm(y ~ idxexposure + mage35, data = metals)
  # step 3) 
  # 3a) the "natural course"
  py_nc = predict(gcompfit_boot)
  # 3b) percent reduction
  # 3b.1) create a copy
  data_reduced = data
  data_reduced[,mixture] = data_reduced[,mixture] * (1-pct_decrease) - metalmeans/metalstds * pct_decrease
  
  # 3b.2) predict the index variable given the post-intervention exposures
  data_reduced$idxexposure = predict(pcafit_boot, newdata = data_reduced[,mixture])[,1]
  
  # 3b.3) make predictions under the interventions data
  py_int = predict(gcompfit_boot, newdata = data_reduced)
  
  # output the effect estimate
  mean(py_int- py_nc) 
}


# confirm this function gives the same estimate as above, when applied to the data
estimate = gcomp_pca(metals)


# bootstrapping: ideally, do 1000+ bootstrap samples
# a loop is not always ideal in R (for speed) but it is helpful here
nbootsamples = 1000
bootres = numeric(nbootsamples)  # pre-allocate a vector for results
for(i in 1:nbootsamples){
  bootdata = metals[sample(1:nrow(metals), replace=TRUE),]
  bootres[i] = gcomp_pca(bootdata)
}


# the standard deviation of the bootstrap samples is the bootstrap standard error
boot_se = sd(bootres)

# 95% confidence intervals
ci = estimate + c(-1.96, 1.96) * boot_se

# final estimates
cat("Mean difference\n", estimate, "\n")
cat(paste("95% CI\n"),ci, "\n")


