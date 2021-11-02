# parameter expansion as in https://wwwf.imperial.ac.uk/~dvandyk/Research/08-jcgs-pxda.pdf
# this example is a simple linear model (non-hierarchical)
#  there is no apparent benefit here - need hierarchical groups
# with high p, this suffers issues with cholesky decompositions
using Distributions, Random, DataFrames, GLM, StatsBase, LinearAlgebra, wellwisejl, SpecialFunctions, CSV
using Plots
plotlyjs()


function sampleinvsigma2(_eta, _sigmaa0, _sigmab0, y, X, _beta, Nf)
# sample precision from a gamma prior
  a = _sigmaa0 + Nf/2.
  res = y .- X * _beta .- _eta 
  se = transpose(res) * res #  permutedims(y .- X * _beta) * (y .- X * _beta)
  b = _sigmab0 + se/2.
  _invsigma2 = rand(rng, Gamma(a, 1. / b))
end

function sampleeta(_beta, _invsigma2, _mu_eta0, _tau_eta0, y, X, Nf)
  yr =  y .- X * _beta 
  iS = Nf * _invsigma2 + _tau_eta0^(-2)
  AMS = sum(yr) * _invsigma2 + _mu_eta0 * _tau_eta0^(-2)
  rand(rng, NormalCanon(AMS, iS)) # parameters mu*Sigma^-1, Sigma^-1
end

function updatebeta!(_beta, _eta, _invsigma2, iLam, _muvec, y, xtx, Xt)
  yr = y .- _eta
  _A = Symmetric(xtx .* _invsigma2 + iLam)
  _A_ms = Xt * yr .* _invsigma2 + iLam * _muvec
  _beta .= rand(rng, MvNormalCanon(_A_ms,_A)) # parameters mu*Sigma^-1, Sigma^-1
end

function updatezeta!(_zeta, _eta, _invsigma2, iLamz, _muvecz, y, _alpha, Wb)
  yr = y .- _eta
  Xs = X * diagm(Wb * _alpha)
  Xst = transpose(Xs)
  xstx = Xst * Xs
  _A = Symmetric(xstx .* _invsigma2 .+ iLamz)
  _A_ms = Xst * yr .* _invsigma2 .+ iLamz * _muvecz
  _zeta .= rand(rng, MvNormalCanon(_A_ms,_A)) # parameters mu*Sigma^-1, Sigma^-1
end

function updatezeta!(_zeta, _invsigma2, iLamz, _muvecz, y, _alpha, Wb)
  yr = y
  Xs = hcat(ones(size(X, 1)), X) * diagm(Wb * _alpha)
  Xst = transpose(Xs)
  xstx = Xst * Xs
  _A = Symmetric(xstx .* _invsigma2 .+ iLamz)
  _A_ms = Xst * yr .* _invsigma2 .+ iLamz * _muvecz
  _zeta .= rand(rng, MvNormalCanon(_A_ms,_A)) # parameters mu*Sigma^-1, Sigma^-1
end


function updatealpha!(_zeta, _eta, _invsigma2, iLama, _muveca, y, _alpha, Wb)
  yr = y .- _eta
  Xs = X * diagm(_zeta) * Wb
  Xst = permutedims(Xs) # use this to keep as a vector
  xstx = Xst * Xs
  _A = Symmetric(xstx .* _invsigma2 .+ iLama)
  _A_ms = Xst * yr .* _invsigma2 .+ iLama * _muveca
  _alpha .= rand(rng, MvNormalCanon(_A_ms,_A)) # parameters mu*Sigma^-1, Sigma^-1
end

function updatealpha!(_zeta, _invsigma2, iLama, _muveca, y, _alpha, Wb)
  yr = y
  Xs = hcat(ones(size(X, 1)), X) * diagm(_zeta) * Wb
  Xst = permutedims(Xs) # use this to keep as a vector
  xstx = Xst * Xs
  _A = Symmetric(xstx .* _invsigma2 .+ iLama)
  _A_ms = Xst * yr .* _invsigma2 .+ iLama * _muveca
  _alpha .= rand(rng, MvNormalCanon(_A_ms,_A)) # parameters mu*Sigma^-1, Sigma^-1
end

function gibbs(y,X,iter)
  N,p = size(X)
  Nf = Float64(N)
  beta_store = zeros(iter, p+1)
  _beta = randn(p)
  _eta = randn()
  #
  Xt = transpose(X)
  xtx = Xt * X
  _muvec = zeros(p)
  iLamd = ones(p) ./ 50.0 
  iLam = Diagonal(iLamd)
  for m in 1:iter
    _invsigma2 = sampleinvsigma2(_eta, 0.0, 0.0, y, X, _beta, Nf)
    _eta = sampleeta(_beta, _invsigma2, 0.0, 100.0, y, X, Nf)
    updatebeta!(_beta, _eta, _invsigma2, iLam, _muvec, y, xtx, Xt)
    beta_store[m,:] = vcat(_eta, _beta)
  end
  DataFrame(beta_store)
end

function gibbs_expand(y,X,iter)
  N,p = size(X)
  pl = ones(Int64, p)
  Wb = I(p) # no batches
  Nf = Float64(N)
  beta_store = zeros(iter, p+1)
  _zeta = rand(p)
  _alpha = rand(size(pl, 1))
  _eta = randn()
  _muvecz = zeros(p)
  _muveca = zeros(size(pl, 1))
  #iLam = Diagonal(ones(p) ./ 1000.0) # original precision prior on beta
  iLamd = ones(p) ./ 50.0           # original precision  prior on beta
  iLama = Diagonal(ones(size(pl, 1)) ./ 100.0) # this doesn't seem to matter except for numerical stability
  iLamz = Diagonal(Array{Float64, 1}(undef, size(pl, 1)))
  _beta = (Wb * _alpha) .* _zeta
  for m in 1:iter
    _invsigma2 = sampleinvsigma2(_eta, 0.0, 0.0, y, X, _beta, Nf)
    _eta = sampleeta(_beta, _invsigma2, 0.0, 100.0, y, X, Nf)
    updatealpha!(_zeta, _eta, _invsigma2, iLama, _muveca, y, _alpha, Wb) 
    # group level variance for beta: 1/sigbeta2 = 1/alpha^2 * 1/sigz; 1/sigz = 1/sigmeta2 * alpha^2
    asq = vcat([fill(_alpha[ii] .^ (2.0), jj) for (ii, jj) in enumerate(pl) ]...)
    # group level variance for zeta in terms of original priors
    iLamz[diagind(iLamz)] = iLamd .* asq
    updatezeta!(_zeta, _eta, _invsigma2, iLamz, _muvecz, y, _alpha, Wb);
    _beta .= (Wb * _alpha) .* _zeta;
    beta_store[m,:] = vcat(_eta, _beta)
  end
  DataFrame(beta_store)
end


function gibbs_expand2(y,X,iter)
# blocking instituted with hardcoded example, includes intercept
  N,p = size(X)
  #pl = ones(Int64, p)
  #Wb =f I(p) # no batches
  pl = [6,5,p-11+1]
  Wb = zeros(p+1, 3)
  Wb[1:pl[1], 1] .= 1.0
  Wb[(pl[1]+1):sum(pl[1:2]), 2] .= 1.0
  Wb[(sum(pl[1:2])+1):sum(pl[1:3]), 3] .= 1.0
  Nf = Float64(N)
  beta_store = zeros(iter, p+1)
  _zeta = rand(p+1)
  _alpha = rand(size(pl, 1))
  _muvecz = zeros(p+1)
  _muveca = zeros(size(pl, 1))
  #iLam = Diagonal(ones(p) ./ 1000.0) # original precision prior on beta
  iLamd = ones(p+1) ./ 50.0           # original precision  prior on beta
  iLama = Diagonal(ones(size(pl, 1)) ./ 10.0) # this doesn't seem to matter except for numerical stability
  iLamz = Diagonal(Array{Float64, 1}(undef, p+1))
  _betaf = (Wb * _alpha) .* _zeta
  for m in 1:iter
    _eta = _betaf[1]
    _beta = _betaf[2:end]
    try
      _invsigma2 = sampleinvsigma2(_eta, 0.0, 0.0, y, X, _beta, Nf)
    catch
      print("s")
    end
    #_eta = sampleeta(_beta, _invsigma2, 0.0, 100.0, y, X, Nf)
    try
      updatealpha!(_zeta, _invsigma2, iLama, _muveca, y, _alpha, Wb) 
    catch
      continue
      print("a")
    end
    # group level variance for beta: 1/sigbeta2 = 1/alpha^2 * 1/sigz; 1/sigz = 1/sigmeta2 * alpha^2
    asq = vcat([fill(_alpha[ii] .^ (2.0), jj) for (ii, jj) in enumerate(pl) ]...)
    # group level variance for zeta in terms of original priors
    iLamz[diagind(iLamz)] = iLamd .* asq
    try
      updatezeta!(_zeta, _invsigma2, iLamz, _muvecz, y, _alpha, Wb);
     catch
      print("z")
    end
   _betaf .= (Wb * _alpha) .* _zeta
    beta_store[m,:] = _betaf
  end
  DataFrame(beta_store)
end

rng = MersenneTwister(sample([Int64(i) for i in 1:1e4]));
  N = 50
  p = 20
  cormat = fill(.8, p, p);
  cormat[diagind(cormat)] .= 1.0;
  X = permutedims(rand(rng, MvNormal(zeros(p), cormat), N));
  beta = sort(sample(rng, [0.1, 3.0, 0.0, 0.0], p, replace=true), rev=true);
  y = rand(rng, Normal(0, 1), N) .* 3 .+ X * beta;
  Xi = hcat(ones(N), X);

ols = fit(LinearModel,Xi,y)
#GLM.dispersion(ols)

res1 = gibbs(y,X,50000)
res2 = gibbs_expand(y,X,50000)
res3 = gibbs_expand2(y,X,50000)


rs1 = summarygibbs(res1[1000:end,:])
rs2 = summarygibbs(res2[1000:end,:])
rs3 = summarygibbs(res3[10:end,:])
ols

plot(hcat(res1.x2[1:100], res2.x2[1:100]))
plot(hcat(res1.x2[1000:5000], res2.x2[1000:5000]))


plot(hcat(rs1.autocor_1, rs2.autocor_1))
plot(hcat(rs1.autocor_5, rs2.autocor_5))
plot(hcat(rs1.ess, rs2.ess))

plot(hcat(rs1.mean, rs2.mean)); plot!( coef(ols)); plot!( vcat(0.0, beta), t=:scatter, label="truth")
plot(hcat(rs1.std, rs2.std)); plot!(stderror(ols))
plot(hcat(rs1.lower2_5, rs2.lower2_5))
plot(hcat(rs1.upper97_5, rs2.upper97_5))

