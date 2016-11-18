#ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"
using JLD
using DataFrames

mu = nothing

function g(K)
    alpha = 1/(mu - 1)
    eps = sqrt(K^2-1)
    return  pi*( eps*(-1 + 2*alpha*eps^2) + K^2*atan(eps))^2 * (2*eps*K*(1 + 2*K^2) + (1 - 4*K^2)*log((K+eps)/(K-eps))) / (K^(7/3)*eps^4*(-3*eps + (2 + K^2)*atan(eps)) )
end

ENV["MPLCONFIGDIR"] = dirname(@__FILE__)
using PyPlot; const plt = PyPlot
info("Config file $(plt.matplotlib[:matplotlib_fname]())")

#using PyCall

#@pyimport matplotlib


# matplotlib[:rcParams]["figure.figsize"] = (5.5,5.5)
#matplotlib[:rcParams]["text.usetex"] = true
# matplotlib[:rcParams]["text.usetex"] = false
# matplotlib[:rcParams]["font.size"]  = 12

# matplotlib[:rcParams]["font.family"] = "sans-serif"
# matplotlib[:rcParams]["font.sans-serif"]  = "Helvetica"

# matplotlib[:rcParams]["font.family"] = "serif"
# matplotlib[:rcParams]["font.serif"]  = "Times"


#plt.rc("text",usetex=true)
# plt.rc("font", family="sans-", serif="Times")
#plt.rc("font", Dict("family"=>"serif","serif"=>"Times"))
# plt.rc("xtick",labelsize=9)
# plt.rc("ytick",labelsize=9)
# plt.rc("axes",labelsize=9)
# plt.rc("legend",fontsize=9)
#plt.rc("image",cmap="gray")
#plt.rc("style",use="ggplot")
# plt.rc("axes",prop_cycle="(cycler('color', 'kkkk') + cycler('linestyle',['-', '--', ':', '-.']))")
# plt.rc("axes",prop_cycle="(cycler('color', 'kkkk') + cycler('linestyle',['-.',':','--','-']))")

#plt.rc("axes",prop_cycle="(cycler('color','k'))")
#plt.rc("axes",prop_cycle="(cycler('color','bgcm'))")

# reds = plt.get_cmap("Reds")
# colstr = join([reds(x) for x in [0.3,0.5,0.7,0.9]],',')
# plt.rc("axes",prop_cycle="(cycler('color', [$colstr]))")



#plt.rc("axes",prop_cycle="(cycler('color', 'kkkk') + cycler('linestyle', ['-', (0,(3,1)),'--','-'] ))")

# using PyCall
# @pyimport pylab as pltfull


#plt.rc("axes",prop_cycle="(cycler('color', 'kkkk') + cycler('dashes',[(5,5),(5,5),(5,5),(5,5)]))")

#plt.rc("axes",prop_cycle="cycler('dashes',[(5,5),(5,5)])")

#['-', '--', ':', '-.']
#plt.rcParams['image.cmap'] = 'gray'  # change default colormap


# width = 3.487
# height = width / 1.618 
fig = plt.figure()

#fig, ax = plt.subplots()
plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$c/R_0$")

K = 1.00001:0.01:10000
mu = 10000
plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=\infty")

K = 1.00001:0.01:1000

mu = 4
plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=4",linestyle=[1,(3,3)])

mu = 10
plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=10",linestyle=[1,(1.5,1.5)])

mu = 100
plt.plot(map(g,K),collect(K).^(-2/3),lw=1,label=L"\mu=100",linestyle=[1,(1.5,1.5,3,1.5)])

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

# function deformability(p)
#     z = []
#     r = []

#     for i in 1:size(p,2)
#         pos = p[:,i]
#         push!(z,pos[1])
#         push!(r,sqrt(pos[2]^2 + pos[3]^2))
#     end

#     return (maximum(z) - maximum(r))/(maximum(r) + maximum(z))
# end

# function deformability(p,axis)
#     if axis==:x
#         p = [p[1,:]; p[2,:]; p[3,:]]
#     elseif axis==:y
#         p = [p[2,:]; p[1,:]; p[3,:]]
#     elseif axis==:z
#         p = [p[3,:]; p[1,:]; p[2,:]]
#     end
#     return deformability(p)
# end

#session = "rotatingfield"

# cases = [57,58,65,66,61,64,55,54]
# lastelement = [-1,-1,-1,-1,30,18,54,294]

#simulation = "sphere0.1:mu=4.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=6.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"

### A line from Cebers
df1 = readtable("$(dirname(@__FILE__))/raimonds/ca_9_short_en_1_200back.dat",separator='\t',header=false)
df2 = readtable("$(dirname(@__FILE__))/raimonds/cb_9_short_en_1_200back.dat",separator='\t',header=false)
#plt.plot(df1[:x1],(df1[:x2].*df2[:x2]).^(1/3),"--k",lw=2)

### From simulation

simulations = ["sphere0.1:mu=4.0;Bm=200.0;omega=0.0", "sphere0.2:mu=10.0;Bm=50.0;omega=0.0","sphere0.2:mu=100.0;Bm=50.0;omega=0.0"]

parameters = [500,21,9]
for (par,simulation) in zip(parameters,simulations)

    indir = "$datadir/ExtractedEqFigures/QstaticFastFieldEiler:$simulation"
    Bm = []
    cR0 = []

    case = 1
    for infile in readdir(indir)

        points = load("$indir/$infile","points")
        faces = load("$indir/$infile","faces")
        Ei = load("$indir/$infile","Ei")

        Bmi = parse(Float64,infile[1:end-4])
        push!(Bm,Bmi)
        
        # Di = deformability(points,:z)
        # push!(D,Di)

        a,b,c = getabc(points)
        push!(cR0,c) ### I should also normalise with volume when R0 is not 1
    end

    #ab = (D + 1)./(1 - D)
    #Dab = (1 + ab).*Derr ./(1 - D)
    #plt.plot(Bm,ab,".")
    #plt.errorbar(Bm[1:4],ab[1:4],yerr=Dab[1:4],fmt=".g")

    ### Pretty plots
    
    Bm, cR0 = Bm[sortperm(Bm)],cR0[sortperm(Bm)]
    #plt.plot(Bm,cR0,".-r",markersize=2)

    Bmf = Float64[]
    cR0f = Float64[]

    deltax = Inf
    deltay = Inf
    j = 1
    push!(Bmf,Bm[1])
    push!(cR0f,cR0[1])

    lg(x) = log(x)/log(10)
    
    for i in 1:length(Bm)
        deltax = lg(Bm[i]) - lg(Bm[j])
        
        if deltax>0.1
            push!(Bmf,Bm[i])
            push!(cR0f,cR0[i])
            j = i
        end
    end
    push!(Bmf,Bm[end])
    push!(cR0f,cR0[end])


    #println("size is $(size(Bmf)) and $(size(cR0f))")
    
    plt.plot(Bmf[Bmf.<par],cR0f[Bmf.<par],"*",markerfacecolor="w")
    plt.plot(Bmf[Bmf.>=par],cR0f[Bmf.>=par],"*",markerfacecolor="k")
    
    #plt.plot(Bmf,cR0f,"*",markersize=5,markerfacecolor=(0, 0, 0.0, 0.1),markeredgecolor=(0, 0, 0, 1.0))
end

plt.plot([],[],"*",markerfacecolor="w",label="simulation")

plt.legend()

#plt.xlim(0,200)
plt.xlim(1,200)
plt.ylim(0,1)
plt.xscale("log")

plt.grid("off")
#plt.tight_layout()
fig[:savefig]("figures/rotatingfieldCR0.pdf")
plt.close("all")
