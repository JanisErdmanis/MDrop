ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

### Need to look closer on case 58 with magnetic bond number 22.4
Bm = [50,72,35,23,20.3,29.26,50.4,68.6,89.6,    17.15, 22.4 , 12.6 ,  8.75, 22.2]
ab = [0.20,0.197,0.268,0.398,NaN,0.307,1,1,1,1,1,1,1, 0.452]
cb = [0.50,0.323,0.591,0.673,NaN,0.63,NaN,NaN,NaN,  0.552401, 0.423585, 0.672603, 0.774839, 0.645]

using DataFrames

df = readtable("raimonds/uabrotminkap9.dat",separator='\t',header=false)

using PyPlot; const plt = PyPlot

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
#plt.rc("legend",fontsize=8)

width = 3.487 
height = width / 1.618

fig = plt.figure(figsize=(width,height))

plt.plot(df[:x1],df[:x2])
plt.plot(Bm[1:6],ab[1:6],"g.")
plt.plot(Bm[10:14],ab[10:14],"g.")
plt.plot(Bm[7:9],ab[7:9],"r.")


plt.xlabel(L" H_0^2 R_0 / \gamma ")
plt.ylabel(L"b/a")
plt.ylim(0,1.1)



plt.grid("off")
plt.tight_layout()
fig[:savefig]("figures/raimondsba.pdf")
plt.close("all")


df = readtable("raimonds/ucbrotminkap9.dat",separator='\t',header=false)

fig = plt.figure(figsize=(width,height))

plt.plot(df[:x1],df[:x2])

plt.plot(Bm[1:6],cb[1:6],"g.")
plt.plot(Bm[10:14],cb[10:14],"g.")
#plt.plot(Bm[7:9],ab[7:9],"r.")


plt.xlabel(L" H_0^2 R_0 / \gamma ")
plt.ylabel(L"c/b")
plt.ylim(0,1.1)



plt.grid("off")
plt.tight_layout()
fig[:savefig]("figures/raimondsac.pdf")
plt.close("all")

