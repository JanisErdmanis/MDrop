using JLD

mu = 10

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

# plt.rc("text",usetex=true)
# plt.rc("font", family="serif", serif="Times")
# plt.rc("xtick",labelsize=10)
# plt.rc("ytick",labelsize=10)
# plt.rc("axes",labelsize=10)
# plt.rc("legend",fontsize=8)

# width = 3.487 
# height = width / 1.618
fig = plt.figure(figsize=(4.1844,2.586*(1 + 0.970-0.835)))
ax1 = fig[:add_subplot](111)
ax2 = ax1[:twiny]()
#ax2[:set_xlim](0,50/0.75)
#ax2[:set_xlabel]
# ax1[:set_xlabel](L"$ H_0^2 R_0/\gamma$")
# ax1[:set_ylabel](L"$b/a$")
# plt.xlabel(L"$ H_0^2 R_0/\gamma$")
# plt.ylabel(L"$b/a$")

### Now adding the exerimental data

plt.sca(ax2)
plt.xlabel(L"H_0^2,~ \rm{Oe}^2")
#x2lim = 30
#x2lim = 260
#x2lim = 260

x1lim = 80
x2lim = x1lim/0.148
plt.xlim(0,x2lim)

plt.sca(ax1)


using DataFrames
dfg = readtable("data/experiment/4_2_rez_dat.dat",separator=',',header=false)
dfg2 = readtable("data/experiment/3_2_rez_dat.dat",separator=',',header=false)
append!(dfg,dfg2)

#plt.plot(df[:x1],df[:x3],".r")
#plt.plot(df[:x1]*0.75,df[:x3],".g")
#plt.plot(df[:x2],df[:x6]./df[:x4])
yerr = dfg[:x6]./dfg[:x4] .* sqrt((dfg[:x5]./dfg[:x4]).^2 + (dfg[:x7]./dfg[:x6]).^2)
scale = x1lim/x2lim
plt.errorbar(dfg[:x2]*scale,dfg[:x6]./dfg[:x4],xerr=dfg[:x3]*scale,yerr=yerr,fmt=".",markersize=3,label="experiment",capsize=1,lw=0.5)

### Theory and simulation

plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$b/a$")

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

#simulation = "sphere0.1:mu=4.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=6.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"
#simulation = "Psphere0.2:mu=10.0;Bm=158.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:36;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:81;omega=0.0"
simulation = "sphere0.2:mu=10.0;Bm=1:84;omega=0.0"

#simulation = "sphere0.2:mu=10.0;Bm=50.0;omega=0.0"
indir = "data/outdir/QstaticFastFieldEiler:$simulation"
#indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEiler:$simulation"

Bm = []
ab = []
cc = []

for infile in readdir(indir)

    points = load("$indir/$infile","points")
    faces = load("$indir/$infile","faces")
    Ei = load("$indir/$infile","Ei")

    Bmi = parse(Float64,infile[1:end-4])
    push!(Bm,Bmi)

    a,b,c = getabc(points)
    #println("Bm=$Bmi a=$a b=$b c=$c")
    push!(ab,b/a)
    push!(cc,c)
end

Bm, ab, cc = Bm[sortperm(Bm)],ab[sortperm(Bm)],cc[sortperm(Bm)]

### Going back
# simulation = "mu10Bm50:mu=10.0;Bm=50.0;omega=0.0"
# indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEilerDecreasing:$simulation"

# Bm = []
# ab = []
# cc = []

# for infile in readdir(indir)

#     points = load("$indir/$infile","points")
#     faces = load("$indir/$infile","faces")
#     Ei = load("$indir/$infile","Ei")

#     Bmi = parse(Float64,infile[1:end-4])
#     push!(Bm,Bmi)

#     a,b,c = getabc(points)
#     #println("Bm=$Bmi a=$a b=$b c=$c")
#     push!(ab,b/a)
#     push!(cc,c)
# end

# Bm, ab, cc = Bm[sortperm(Bm)],ab[sortperm(Bm)],cc[sortperm(Bm)]
# plt.plot(Bm,ab,"-k")


# df = readtable("$(dirname(@__FILE__))/raimonds/ab_9_short_en_1_200.dat",separator='\t',header=false)
# plt.plot(df[:x1],df[:x2],"-")

df = readtable("data/ellipsoid/ab_9_short_en_1_200.dat",separator='\t',header=false)
plt.plot(df[:x1],df[:x2],"-",label="ellipsodial approximation")

### PLot in front

plt.plot([Bm[1:5:15];Bm[16:30];Bm[30:5:end]],[ab[1:5:15];ab[16:30];ab[30:5:end]],"*",markerfacecolor="w",label="quasistatic simulation")
#plt.plot([Bm[1:4:15];Bm[16:26];Bm[27:1:end]],[ab[1:4:15];ab[16:26];ab[27:1:end]],"*",markerfacecolor="w",label="quasistatic simulation")

Bm = [30 40 18 19 20 25 60 27 22.5]
ba = [0.1950 0.1638 0.5676 0.4180 0.3545 0.2379 0.1578 0.2220 0.2772]
plt.plot(Bm,ba,"o",markerfacecolor=(1, 1, 1, 0.7))
plt.plot([],[],"o",markerfacecolor=(1, 1, 1, 0.7),label="ellipsoid relaxation")

plt.xlim(0,x1lim)
plt.ylim(0,1.1)

plt.legend()

plt.subplots_adjust(top=0.835)
plt.subplots_adjust(left=0.14,right=0.97)

fig[:savefig]("figures/rotatingfieldAB.pdf")
plt.close("all")
