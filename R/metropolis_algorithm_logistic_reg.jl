######################################################################################################################
# Author: Alex Keil
# Program: metropolis_algorithm_logistic_reg.jl
# Language: Julia v1.7
# Date: Wednesday, February 16, 2022 at 11:03:12 AM
# Project: Coding examples
# Tasks: Simple example of using Metropolis algorithm for Bayesian fitting of logistic 
#  regression model
# Released under the GNU General Public License: http://www.gnu.org/copyleft/gpl.html
######################################################################################################################


using Random, DataFrames, Distributions, GLM, StatsBase, LinearAlgebra



function expit(mu::T) where T
 p::T = one(T) /(one(T) + exp(- mu))
 p
end

function dgm(N;linear=true,seed=1231)
  rng = Random.MersenneTwister(seed)
  z = rand(rng, Float64,N)
  x = rand.(rng, Bernoulli.(expit.(-0.5 .+ z)))
  lpy = -1.0 .+ x + z
  y = rand.(rng, Bernoulli.(expit.(lpy)))
  z,x,y
end


function dgm2(N;linear=true,seed=1231)
  rng = Random.MersenneTwister(seed)
  z = rand(rng, Float64,N)
  x = rand.(rng, LogNormal.(-0.5 .+ z, 0.3)) # distribution prone to outliers
  lpy = -1.0 .+ x + z
  y = rand.(rng, Bernoulli.(expit.(lpy)))
  z,x,y
end


function mydiag(X) where T
  if size(X,2)==1
    return X
  else 
    return diag(X)
  end
end

function mhlogit(y,X;iter=100, burnin=0, thin=1, printrat=true,bs=0.3,binits=randn(3)*2, adapt=true, keepwarmup=false, seed=sample([Int(i) for i in 1:1e8]))
  rng = Random.MersenneTwister(seed)
  n = size(y,1)
  p = size(X,2) 
  # inits
  _beta = binits
  # acceptors
  saB = zero(Float64)
  #priors
  
  # proposal sds
  covB = Symmetric(Array(I,p,p)*bs)
  #storage
  lp_store = Array{Float64, 1}(undef, iter)             # proportional to log posterior
  beta_store = Array{Float64, 2}(undef, iter, p)
  beta_store[1,:]    = _beta
  #allocation
  mu = expit.(X * _beta)
  muc = expit.(X * _beta)
  propB = zero(_beta)  
  _betac = zero(_beta)  
  adaptscale = sqrt(5.76/p)
  for j in 2:iter
    ############################
    # adaptation phase
    ############################
    if adapt && j < burnin && (j % 5)==0 && j > 100
       covB .= Symmetric(adaptscale .* cov(beta_store[1:(j-1),:]))
    end
    ############################
    # proposal values
    ############################
    if size(covB,2)>1 && isposdef(covB)
      propB[:] .= rand(rng,MvNormal(zeros(p), covB))
    else
      propB[:] .= rand.(rng,Normal.(zeros(p), sqrt.(mydiag(covB))))
    end
    ############################
    # block update model parameters
    ############################
    _betac[:] .= _beta + propB
    muc[:] .= expit.(X * _betac)
    mu[:] .= expit.(X * _beta )
    # diffuse normal prior
    lpc  = sum(logpdf.(Bernoulli.(muc), y)) + sum(logpdf.(Normal(0,10000.), _betac))
    lpo  = sum(logpdf.(Bernoulli.(mu), y))  + sum(logpdf.(Normal(0,10000.), _beta))
    # uniform prior
    #lpc  = sum(logpdf.(Bernoulli.(muc), y)) 
    #lpo  = sum(logpdf.(Bernoulli.(mu), y)) 
    ## acceptance probability ##
    aB = min(1.0, exp(lpc-lpo))
    saB += aB 
    # accept/reject
    if aB > rand(rng)
      lp_store[j] = lpc
      _beta[:] .= _betac
    else
      lp_store[j] = lpo
    end
    beta_store[j,:]    = _beta
  end # MCMC iteration loop
  if printrat
    println("Beta acceptance ratio ", saB/iter)
  end
  startkeep = keepwarmup ? 1 : burnin+1
  keepidx = range(startkeep, iter, step=thin)
  df = DataFrame(hcat([i for i in 1:iter], lp_store, beta_store), vcat(
       :iter, :logpost,[Symbol("b" * "_$i") for i in 0:(p-1)]))
  return (
      df[keepidx,:], 
      (:B => saB/iter)
    )
end


function summarycol(col)
   means = mean(col)
   medians = median(col)
   pl = quantile(col, 0.025)[1]
   pu = quantile(col,  0.975)[1]
   stds =  std(col)
   ac1, ac5 =  autocor(col, [1,5])
   lens = length(col)
   means, medians, pl, pu, stds, ac1, ac5, lens
end


function summarymcmc(results::DataFrame)
 ns = size(results)
 outres = Array{Float64,2}(undef, ns[2], 8)
 nm = names(results)
 for i in 1:ns[2]
   outres[i,:] .= summarycol(results[:,i])
 end
 DataFrame(hcat(nm, outres), 
           [:nm, :mean, :std, :median, :lower2_5, :upper97_5, :autocor_1, :autocor_5, :length])
end

# generate data
#z,x,y = dgm(150, seed=1232) # simple example with binary exposure
z,x,y = dgm2(150, seed=1232) # this highlights some difference in Bayes vs. ML
X = hcat(ones(length(x)), x, z)

# fit models
mx = glm(@formula(y~x+z), DataFrame(x=x,y=y,z=z), Binomial());
res, acc = mhlogit(y,X, seed=1232, iter=100000, burnin=10000, thin=10,bs=0.15,binits=zeros(size(X,2)));
reso = summarymcmc(res[:,2:end]);

### Maximum likelihood
println("### ML")
println(mx)
println("log-likelihood: $(loglikelihood(mx))")
# y ~ 1 + x + z
# 
# Coefficients:
# ──────────────────────────────────────────────────────────────────────────
#                  Coef.  Std. Error      z  Pr(>|z|)   Lower 95%  Upper 95%
# ──────────────────────────────────────────────────────────────────────────
# (Intercept)  -0.738985    0.533878  -1.38    0.1663  -1.78537     0.307397
# x             1.37353     0.716591   1.92    0.0553  -0.0309611   2.77803
# z            -0.128342    0.887125  -0.14    0.8850  -1.86707     1.61039
# ──────────────────────────────────────────────────────────────────────────
# log-likelihood: -94.22718845397189


### Bayes
println("### Bayes")
println(reso)
#  Row │ nm       mean       std        median    lower2_5  upper97_5  autocor_1  autocor_5    length 
#      │ Any      Any        Any        Any       Any       Any        Any        Any          Any    
# ─────┼──────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ logpost  -126.129   -125.8     -129.371  -124.729  1.23374    0.177401   0.0070398    9000.0
#    2 │ b_0      -0.797824  -0.783011  -1.88873  0.224845  0.537868   0.131137   0.00813639   9000.0
#    3 │ b_1      1.46708    1.42894    0.127659  2.93323   0.720429   0.140184   -0.00184611  9000.0
#    4 │ b_2      -0.171078  -0.17837   -1.95199  1.57014   0.892912   0.156661   -0.0126926   9000.0


### Bayes maximum a posteriori
println("### MCMC based MAP")
println(res[argmax(res.logpost),:])
#   Row │ iter     logpost   b_0        b_1      b_2       
#       │ Float64  Float64   Float64    Float64  Float64   
# ──────┼──────────────────────────────────────────────────
#  8312 │ 93111.0  -124.616  -0.723186  1.35561  -0.106778


# uncomment for plotting
#using StatsPlots
#plot(res.b_1)
#density(res.b_1)