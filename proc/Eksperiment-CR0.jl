using JLD

mu = 10

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

using DataFrames

#append!(dfg,dfg2)

plt.figure()
plt.xlabel(L"H_0^2,~\rm{Oe}^2")
plt.ylabel(L"c/R_0")
# 
# plt.plot(dfg[:x2],dfg[:x6])

dfg = readtable("data/experiment/3_2_rez_dat.dat",separator=',',header=false)

dfg1 = dfg[1:20,:]
dfg1[:c] = 1./dfg1[:x4]./dfg1[:x6]

dfg2 = dfg[21:41,:]
dfg2[:c] = 1./dfg2[:x4]./dfg2[:x6]*0.94

dfgg = readtable("data/experiment/4_2_rez_dat.dat",separator=',',header=false)

dfg3 = dfgg[1:20,:]
dfg3[:c] = 1./dfg3[:x4]./dfg3[:x6]


df = vcat(dfg1,dfg2,dfg3)
df[:error] = sqrt(df[:x5]./df[:x4].^2 + df[:x7]./df[:x6].^2).*df[:c]

#plt.plot(df[:x2],df[:c],".")

plt.errorbar(df[:x2],df[:c],yerr=df[:error]/20,fmt=".",markersize=3,label="experiment",capsize=1,lw=0.5)

H02 = 0:1:120
plt.plot(H02,1-H02/100*0.235,"-k",label="fit")

plt.xlim(0,110)
plt.ylim(0.75,1.05)

plt.legend()
plt.savefig("figures/experimentCR0.pdf")
plt.close("all")


