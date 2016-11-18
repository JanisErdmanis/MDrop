using JLD
using DataFrames

mu = 10

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

using DataFrames
dfg = readtable("$(dirname(@__FILE__))/guntars/elongfit.dat",separator=',',header=false)
dfg = dfg[2:end,:]

fig = plt.figure()
plt.xlabel(L"$ H_0^2 $")
plt.ylabel(L"$b/a$")

plt.plot(dfg[:x1].^2,dfg[:x2]./dfg[:x3],".")

mu = 11
alpha = 1/(mu-1)
function g(K)
    eps = sqrt(1-K^2)
    n = K^2 * ( -2*eps + log((1+eps)/(1-eps)))/2/eps^3

    8*pi*(n + alpha)^2 * eps^2 * K^(-4/3) * (1 + 2*K^2 + (1 - 4*K^2)/eps/K *asin(eps))/(-6 + (2 + K^2)/eps*log((1+eps)/(1-eps)))
end

K = 0.05:0.001:0.99
plt.plot(map(g,K)/0.24,1./collect(K),"-k")

#plt.xlim(0,200)

plt.savefig("figures/eksp-elonganation.pdf")
plt.close("all")
