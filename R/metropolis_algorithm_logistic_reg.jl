# simulating multiplicative Berkson measurement error
# updated 2/22 to fix bug in exposure updates

#!julia -p 8
using Random, StatsPlots, DataFrames, Distributions, GLM, StatsBase, LinearAlgebra
plotly()


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
  x = rand.(rng, LogNormal.(-0.5 .+ z, 0.3))
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

function mhlogit(y,X;iter=100, burnin=0, thin=1, printrat=true,bs=0.3,binits=randn(3)*2, adapt=true, seed=sample([Int(i) for i in 1:1e8]))
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
  lp_store = Array{Float64, 1}(undef, iter)
  beta_store = Array{Float64, 2}(undef, iter, p)
  beta_store[1,:]    = _beta
  #allocation
  mu = expit.(X * _beta)
  muc = expit.(X * _beta)
  propB = zero(_beta)  
  _betac = zero(_beta)  
  
  for j in 2:iter
    ############################
    # adaptation phase
    ############################
    if adapt && j < burnin && (j % 5)==0 && j > 100
       covB .= Symmetric(sqrt(5.76/(p)) .* cov(beta_store[1:(j-1),:]))
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
    aB = min(1.0, exp(lpc-lpo))
    saB += aB 
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
  df = DataFrame(hcat([i for i in 1:iter], lp_store, beta_store), vcat(
       :iter, :logpost,[Symbol("b" * "_$i") for i in 0:(p-1)]))
  return (
      df[range(burnin+1, iter, step=thin),:], 
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


#z,x,y = dgm(150, seed=sample([Int(i) for i in 1:1e8]))
z,x,y = dgm2(150, seed=sample([Int(i) for i in 1:1e8])) # this highlights some difference in Bayes vs. ML
X = hcat(ones(length(x)), x, z)
println(sum(y)/length(y))

mx = glm(@formula(y~x+z), DataFrame(x=x,y=y,z=z), Binomial())

res, acc = mhlogit(y,X, iter=100000, burnin=10000, thin=10,bs=0.15,binits=zeros(size(X,2)));
reso = summarymcmc(res[:,2:end]);

# Maximum likelihood
println("ML")
println(mx)
println("log-likelihood: $(loglikelihood(mx))")
# Bayes
println("Bayes")
println(reso)
# Bayes maximum a posteriori
println("MCMC based MAP")
println(res[argmax(res.logpost),:])


plot(res.b_1)
density(res.b_1)