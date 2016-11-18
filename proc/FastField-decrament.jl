using JLD
using DataFrames

mu = 10

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

#fig, ax = plt.subplots()
fig = plt.figure()

plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L" \eta R_0/\gamma \tau")

### Now adding the exerimental data

#simulation = "sphere0.1:mu=4.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=6.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"
#simulation = "Psphere0.2:mu=10.0;Bm=158.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=50.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:84;omega=0.0"
simulation = "sphere0.1:mu=10.0;Bm=1:84;omega=0.0"

indir = "data/outdir/QstaticFastFieldEiler:$simulation"

Bm = Float64[]
tau = Float64[]
# E = []
# Etheor = []

for infile in readdir(indir)

    # points = load("$indir/$infile","points")
    # faces = load("$indir/$infile","faces")
    taui = load("$indir/$infile","taui")

    Bmi = parse(Float64,infile[1:end-4])
    push!(Bm,Bmi)

    # a,b,c = getabc(points)
    # factor = (volume(points,faces)/(4/3*pi*a*b*c))^(1/3)
    # println("factor = $factor")
    
    #Etheori = TheoreticalDropEnergy(factor*a,factor*b,factor*c,mu,Bmi)
    # Etheori = TheoreticalDropEnergy(a,b,c,mu,Bmi)
    # push!(Etheor,Etheori)
    
    #println("Bm=$Bmi a=$a b=$b c=$c Area=$(EllipsoidArea(a,b,c))")
    # I could add line which corresponds to ellipsoid for given parameters

    push!(tau,taui)
    #push!(E,Ei)    
end

#Bm, E, Etheor = Bm[sortperm(Bm)],E[sortperm(Bm)], Etheor[sortperm(Bm)]
#Bm, tau = Bm[sortperm(Bm)][1:3:end],tau[sortperm(Bm)][1:3:end]
Bm, tau = Bm[sortperm(Bm)][1:1:end],tau[sortperm(Bm)][1:1:end]
#plt.plot(Bm[1:end-11],1./tau[1:end-11],"*",markerfacecolor="w",label="simulation")
plt.plot(Bm[1:2:end],1./tau[1:2:end],"*",markerfacecolor="w",label="simulation")

Bm = collect(0:1:25)
plt.plot(Bm,0.49 - 0.025*Bm,"-",label="extrapolation")
#plt.plot(Bm,1./tau,".-k")

### Reference energy

# Here I would get one line from Cebers for comparssion
# using DataFrames
# df = readtable("$(dirname(@__FILE__))/guntars/20160613_4_2.dat",separator=',',header=false)
# plt.plot(df[:x1]*0.75,df[:x3],".r")

plt.legend()

plt.xlim(0,20)
#plt.xlim(0,30)
plt.ylim(0,0.5)

#plt.grid("off")
#plt.tight_layout()
fig[:savefig]("figures/decrament.pdf")
plt.close("all")
