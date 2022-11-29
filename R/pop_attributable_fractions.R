######################################################################################################################
# Author: Alex Keil
# Program: pop_attributable_fractions.R
# System notes: R 4.0.2, MacOS
# Date: Nov 29, 2022
# Tasks: valid calculation of population attributable fractions in uncensored data, under confounding
# Description: program simulates data for a small cross-sectional/single time point analysis
# Using a potential outcome defintion of population attributable fractions (PAF), it calculates a PAF 
# as:
#   PAF = (r-r0)/r
# where:
#   r0 = average risk of disease, had everyone been unexposed
# r  = average risk of disease in the population
# 
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
# Next steps: this yields a point estimate only. confidence intervals can be obtained via bootstrapping
######################################################################################################################/
# simulate some data for an example (use your real data);
library(boot)

dgm <- function(n){
  z <- rbinom(n, 1, 0.5)
  x <- rbinom(n, 1, 0.2+0.2*z)
  y <- rbinom(n, 1, 0.2 + 0.2*x + 0.2*z)
  data.frame(exposure=x, confounder=z, outcome=y)
}

original_data <- exposed_data <-  unexposed_data <-  dgm(1000)
exposed_data$exposure=1
unexposed_data$exposure=0

# model on simulated data (use your real model);
mfit <- glm(outcome ~ exposure + confounder, data=original_data, family=binomial(link="log")) #POISSON COULD BE USED, but check that predictions are < 1.0 - can also use LOGIT link
pred_population <- predict(mfit, type="response")
pred_exposed <- predict(mfit, exposed_data, type="response")
pred_unexposed <- predict(mfit, unexposed_data, type="response")

risk_population = mean(pred_population)
risk_exposed = mean(pred_exposed)
risk_unexposed = mean(pred_unexposed)

# marginal risk difference
(rd <- (risk_exposed - risk_unexposed))

# Population attributable fraction
(par <- (risk_population - risk_unexposed)/risk_population)

