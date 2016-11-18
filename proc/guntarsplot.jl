ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

using JLD
using DataFrames

fig = plt.figure()

plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"axis")

# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=\infty")
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=4",linestyle=[1,(3,3)])
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=10",linestyle=[1,(1.5,1.5)])
# plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=100",linestyle=[1,(1.5,1.5,3,1.5)])

# df = readtable("$(dirname(@__FILE__))/cebers/ab_5_short_en_1_150.dat",separator='\t',header=false)
# plt.plot(df[:x1],df[:x2],label=L"\mu=6",linestyle=linestyle=[1,(3,3)])

df1 = readtable("$(dirname(@__FILE__))/cebers/ca_9_short_en_1_200.dat",separator='\t',header=false)
df2 = readtable("$(dirname(@__FILE__))/cebers/cb_9_short_en_1_200.dat",separator='\t',header=false)

Bm = df1[:x1]
c = (df1[:x2].*df2[:x2]).^(1/3)
a = (c ./ df1[:x2])
b = (c ./ df2[:x2])

plt.plot(Bm,a,"-",label=L"a/R_0")
plt.plot(Bm,b,label=L"b/R_0",linestyle=[1,(1.5,1.5)])
plt.plot(Bm,c,label=L"c/R_0",linestyle=[1,(3,3)])


### Simulation

using SurfaceGeometry

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
simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"
#simulation = "Psphere0.2:mu=10.0;Bm=158.0;omega=0.0"

#simulation = "sphere0.2:mu=10.0;Bm=50.0;omega=0.0"
indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEiler:$simulation"

Bm = []
a = []
b = []
c = []

for infile in readdir(indir)

    points = load("$indir/$infile","points")
    faces = load("$indir/$infile","faces")
    Ei = load("$indir/$infile","Ei")

    Bmi = parse(Float64,infile[1:end-4])
    push!(Bm,Bmi)

    a_,b_,c_ = getabc(points)
    vol = volume(points,faces)
    R0 = (3*vol/4/pi)^(1/3)
    
    #println("Bm=$Bmi a=$a b=$b c=$c")
    push!(a,a_/R0)
    push!(b,b_/R0)
    push!(c,c_/R0)
    # push!(cc,c)
end

Bm, a, b, c = Bm[sortperm(Bm)],a[sortperm(Bm)],b[sortperm(Bm)],c[sortperm(Bm)]

Bm = [Bm[1:4:15];Bm[16:4:26];Bm[27:4:end]]
a = [a[1:4:15];a[16:4:26];a[27:4:end]]
b = [b[1:4:15];b[16:4:26];b[27:4:end]]
c = [c[1:4:15];c[16:4:26];c[27:4:end]]

plt.plot(Bm,a,"*",markerfacecolor="w",label="simulation")
plt.plot(Bm,b,"*",markerfacecolor="w")
plt.plot(Bm,c,"*",markerfacecolor="w")

plt.plot(Bm,a.*b.*c,"*",markerfacecolor="r",label=L"a*b*c/R_0^3")

plt.plot(Bm,1./a./b,"*",markerfacecolor="g",label=L"R_0^2/a b")

plt.plot(Bm,2./a./b/3,"*",markerfacecolor="b",label=L"2 R_0^2/ 3 a b")

# df = readtable("$(dirname(@__FILE__))/cebers/ab_16_short_en_1_200c.dat",separator='\t',header=false)
# plt.plot(df[:x1],df[:x2],"-",label=L"\mu=17")

# df = readtable("$(dirname(@__FILE__))/cebers/ab_16_short_en_1_200backc.dat",separator='\t',header=false)
# plt.plot(df[:x1],df[:x2],"-")

plt.xlim(0,50)
plt.ylim(0.0,1.1)
plt.legend(loc=3)
plt.savefig("figures/guntarsplot.pdf")
plt.close("all")
