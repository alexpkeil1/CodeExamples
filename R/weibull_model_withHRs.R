library(survival)


# simulating potential outcome times for exposed/unexposed where
# the exposed potential time (t1) is a fixed function of the unexposed potential time (t0)


####### exponential model #######
shape0 = 1
scale0 = 5

n = 10000
HR = 2 # hazard ratio for X on the outcome
x <- rbinom(n, 1, 0.5)
t0 <- rweibull(n, shape = shape0, scale = scale0)
t1 <- exp(log(t0) - log(HR))
t <- ifelse(x, t1, t0)

# should give HR = 2
log(HR)
coxph(Surv(t)~x)

# log aft parameter = log(1/HR) for exponential
survreg(Surv(t)~x, dist="exponential")$coefficients[2]
log(1/HR)



####### weibull model #######

shape0 = 3
scale0 = 5

n = 10000
HR = 2 # hazard ratio for X on the outcome
x <- rbinom(n, 1, 0.5)
t0 <- rweibull(n, shape = shape0, scale = scale0)
t1 <- exp(log(t0) - log(HR)/(shape0))
t <- ifelse(x, t1, t0)

log(HR)
coxph(Surv(t)~x)

# log aft parameter = HR^(-1/shape0)
survreg(Surv(t)~x, dist="weibull")$coefficients[2]
log(HR^(-1/shape0))
