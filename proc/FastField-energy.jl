using JLD
using DataFrames
using SurfaceGeometry

mu = 10

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

#fig, ax = plt.subplots()

fig = plt.figure()

plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$(E_M + E_S)/4\pi R_0^2 \gamma$")
#plt.legend()


### Now adding the exerimental data

function getabc(points)

    ar = 0
    al = 0
    br = 0
    bl = 0
    cr = 0
    cl = 0
    
    for xkey in 1:size(points,2)
        x = points[:,xkey]
        x[1]<al && (al=x[1])
        x[1]>ar && (ar=x[1])
        x[2]<bl && (bl=x[2])
        x[2]>br && (br=x[2])
        x[3]<cl && (cl=x[3])
        x[3]>cr && (cr=x[3])
    end

    a = (ar - al)/2
    b = (br - bl)/2
    c = (cr - cl)/2

    return a,b,c
end

function ellipsoid_demagnetization_coefficients(a,b,c)

    UP_BOUND = 1000

    Ru2(u) = (u+a^2)*(u+b^2)*(u+c^2)

    nx = 1/2 * a*b*c * quadgk(s -> 1/(s+a^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    ny = 1/2 * a*b*c * quadgk(s -> 1/(s+b^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    nz = 1/2 * a*b*c * quadgk(s -> 1/(s+c^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]

    return [nx, ny, nz]
end

function EllipsoidField(a,b,c,mu,H0)

    H0x, H0y, H0z = H0
    nx, ny, nz = ellipsoid_demagnetization_coefficients(a,b,c)

    Hix = H0x/(1 + (mu-1)*nx)
    Hiy = H0y/(1 + (mu-1)*ny)
    Hiz = H0z/(1 + (mu-1)*nz)

    return [Hix,Hiy,Hiz]
end

function TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Hx = EllipsoidField(a,b,c,mup,[H0,0,0])
    Hy = EllipsoidField(a,b,c,mup,[0,H0,0])

    Emt = 1/8/pi*(1 - mup) * (dot(Hx,H0*[1,0,0]) + dot(Hy,H0*[0,1,0]))/2 * 4/3*pi * a * b * c
    return Emt
end

import Elliptic
function EllipsoidArea(a,b,c)

    if a<b
        a,b = b,a
    end
    if b<c
        b,c = c,b
    end
    if a<b
        a,b = b,a
    end

    cosphi = c/a
    phi = acos(cosphi)
    sinphi = sqrt(1 - cosphi^2)
    k = sqrt(a^2*(b^2-c^2)/b^2/(a^2-c^2))

    S = 2*pi*c^2 + 2*pi*a*b/sinphi * (sinphi^2 * Elliptic.E(phi,k) + cosphi^2 * Elliptic.F(phi,k))

    return S
end

function TheoreticalDropEnergy(a,b,c,mup,Bm)
    ### I should scale a,b,c inside
    
    gammap = 1.
    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    Emt = TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Etotal = Emt + gammap*EllipsoidArea(a,b,c) # + surface area
end

#simulation = "sphere0.1:mu=4.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=6.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"
#simulation = "Psphere0.2:mu=10.0;Bm=158.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:36;omega=0.0"

#simulation = "sphere0.2:mu=10.0;Bm=50.0;omega=0.0"
simulation = "sphere0.2:mu=10.0;Bm=1:84;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:81;omega=0.0"
#indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEiler:$simulation"
indir = "data/outdir/QstaticFastFieldEiler:$simulation"

# simulation = "mu10Bm50:mu=10.0;Bm=50.0;omega=0.0"
# indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEilerDecreasing:$simulation"

Bm = []
E = []
Etheor = []

points = nothing
faces = nothing

for infile in readdir(indir)

    points = load("$indir/$infile","points")
    faces = load("$indir/$infile","faces")
    Ei = load("$indir/$infile","Ei")

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
    
    push!(E,Ei)    
end

#Bm, E, Etheor = Bm[sortperm(Bm)],E[sortperm(Bm)], Etheor[sortperm(Bm)]
Bm, E = Bm[sortperm(Bm)],E[sortperm(Bm)]

### Reference energy
E0 = 4*pi*(3*volume(points,faces)/4/pi)^(2/3)
#plt.plot(Bm,E/E0,"-r",label="simulation")

#####

# df1 = readtable("data/outdir/ellipsoid/ca_9_short.dat",separator='\t',header=false)
# df2 = readtable("data/outdir/ellipsoid/cb_9_short.dat",separator='\t',header=false)

# R0 = 1
# E0 = 4*pi*R0^2 ### gamma=1 (used in theoretical Drop energy calc)
# Bmtheor = []
# Etheor = []
# for i in 1:size(df1,1)
#     Bmi = df1[:x1][i]
#     ca = df1[:x2][i]
#     cb = df2[:x2][i]

#     c = R0*(ca*cb)^(1/3)
#     a = c/ca
#     b = c/cb

#     Etheori = TheoreticalDropEnergy(a,b,c,mu,Bmi)
#     push!(Etheor,Etheori)
#     push!(Bmtheor,Bmi)
# end

#plt.plot(Bmtheor,Etheor/E0,label="ellipsodial energy")

###

df3 = readtable("data/ellipsoid/energy_9_short_en_1_200.dat",separator='\t',header=false)

plt.plot(df3[:x1],df3[:x2]./2,label="minimal ellipsodial energy",lw=1.5)


function EofBm(Bm)
    i = 1
    for i in 1:size(df3,1)
        if df3[:x1][i]>Bm
            break
        end
    end
    #println("Bm is $(df3[:x1][i])")
    return df3[:x2][i]./2
end

                # Bmtheor = []
# Etheor = []

# df4 = readtable("$(dirname(@__FILE__))/raimonds/myatempt.dat",separator='\t',header=true)
# plt.plot(df4[:Bm],df4[:E],"-k")


# Here I would get one line from Cebers for comparssion
# using DataFrames
# df = readtable("$(dirname(@__FILE__))/guntars/20160613_4_2.dat",separator=',',header=false)
# plt.plot(df[:x1]*0.75,df[:x3],".r")

### Plotting simulation afterwards

plt.plot(Bm[1:3:end],E[1:3:end]/E0,"*",markerfacecolor="w",label="simulation")
#plt.errorbar(Bm[1:7:end],E[1:7:end]/E0,yerr=0.17,fmt="*",markerfacecolor="w",label="simulation")

Bm = [30 40 18 19 20 60 27 22.5 25][:] #47
Efin = [-0.2217 -0.8356 0.4060 0.3632 0.3168 -2.1905 -0.0551 0.1919 0.0581][:] #-1.4512
Ein = [-0.2172 -0.837 0.4059 0.3628 0.3165 -2.2269 -0.0567 0.1918 0.0580][:] #-1.4512

yerr = [abs(EofBm(Bm[i]) - Ein[i]) for i in 1:length(Bm)]


plt.errorbar(Bm,Efin,fmt="o",yerr=yerr,markerfacecolor=(1, 1, 1, 0.7),label="ellipsoid relaxation")

#plt.errorbar(Bm,Efin,fmt=".",markersize=3,label="experiment",capsize=1,lw=0.5)


#plt.plot([],[],"o",markerfacecolor=(1, 1, 1, 0.7),label="ellipsoid relaxation")

# Bm = [30 40 18 19 20 60 27 22.5 25]
# Ein = [-0.2172 -0.837 0.4059 0.3628 0.3165 0.7288 -2.2269 -0.0567 0.1918 0.0580]
# plt.plot(Bm,Efin,"o",markerfacecolor=(1, 1, 1, 0.7))
# plt.plot([],[],"o",markerfacecolor=(1, 1, 1, 0.7),label="numerical calculation of ellipsoid energy")

plt.legend()

plt.xlim(0,65)
plt.ylim(-2.5,1)

# plt.xlim(17,20)
# plt.ylim(0.3,0.5)


fig[:savefig]("figures/rotatingfieldenergy.pdf")
plt.close("all")
