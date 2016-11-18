ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

using DataFrames

fig = plt.figure()

plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$b/a$")

# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=\infty")
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=4",linestyle=[1,(3,3)])
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=10",linestyle=[1,(1.5,1.5)])
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=100",linestyle=[1,(1.5,1.5,3,1.5)])

df = readtable("$(dirname(@__FILE__))/cebers/ab_5_short_en_1_150.dat",separator='\t',header=false)
plt.plot(df[:x1],df[:x2],label=L"\mu=6",linestyle=linestyle=[1,(3,3)])

df = readtable("$(dirname(@__FILE__))/cebers/ab_9_short_en_1_200.dat",separator='\t',header=false)
plt.plot(df[:x1],df[:x2],label=L"\mu=10",linestyle=linestyle=[1,(1.5,1.5)])

df = readtable("$(dirname(@__FILE__))/cebers/ab_16_short_en_1_200c.dat",separator='\t',header=false)
plt.plot(df[:x1],df[:x2],"-",label=L"\mu=17")

df = readtable("$(dirname(@__FILE__))/cebers/ab_16_short_en_1_200backc.dat",separator='\t',header=false)
plt.plot(df[:x1],df[:x2],"-")

plt.xlim(0,95)
plt.ylim(0,1.1)
plt.legend(loc=4)
plt.savefig("figures/cebersplot.pdf")
plt.close("all")
