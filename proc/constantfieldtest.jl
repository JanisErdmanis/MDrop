ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

alpha = 1/9

function g(K)
    eps = sqrt(1-K^2)
    n = K^2 * ( -2*eps + log((1+eps)/(1-eps)))/2/eps^3

    8*pi*(n + alpha)^2 * eps^2 * K^(-4/3) * (1 + 2*K^2 + (1 - 4*K^2)/eps/K *asin(eps))/(-6 + (2 + K^2)/eps*log((1+eps)/(1-eps)))
end


using PyPlot; const plt = PyPlot
# @pyimport matplotlib as mpl
# mpl.use("pdf")

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
#plt.rc("legend",fontsize=8)

width = 3.487 
height = width / 1.618
fig = plt.figure(figsize=(width,height))

K = 0.001:0.001:0.99
plt.plot(map(g,K),1./collect(K))
plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$c/a$")
#plt.legend()


### Now adding the exerimental data

function deformability(p)
    z = []
    r = []

    for i in 1:size(p,2)
        pos = p[:,i]
        push!(z,pos[1])
        push!(r,sqrt(pos[2]^2 + pos[3]^2))
    end

    return (maximum(z) - maximum(r))/(maximum(r) + maximum(z))
end

function deformability(p,axis)
    if axis==:x
        p = [p[1,:]; p[2,:]; p[3,:]]
    elseif axis==:y
        p = [p[2,:]; p[1,:]; p[3,:]]
    elseif axis==:z
        p = [p[3,:]; p[1,:]; p[2,:]]
    end
    return deformability(p)
end

session = "constantfield4"
calculate = false
viewmeshactive = true  ### avoiding of plotting

#cases = [1,2,3,4,5,6,7,9,10,11,12,13,14,15]
#cases = [1,2,3,15,4,5,9,10]
cases = [16,17,18,19,20,21]

Bm = []
D = []
Derr = []

case = 1
for case in cases
    include("diaryconst.jl")
    n = length(memory)
    push!(Bm,norm(H0)^2 * (a*b*c)^(1/3)/gammap)

    Di = deformability(memory[n][2],:x)
    h = tmax/N
    Derri = (Di - deformability(memory[n-1][2],:x))/h * tmax

    push!(Derr,Derri)
    push!(D,Di)
end

ab = (D + 1)./(1 - D)
Dab = (1 + ab).*Derr ./(1 - D)

#plt.plot(Bm,ab,".")
plt.errorbar(Bm,ab,yerr=Dab,fmt=".g")
plt.xlim(0,25)
plt.ylim(0,10)

plt.grid("off")
plt.tight_layout()
fig[:savefig]("figures/constantfield.pdf")
plt.close("all")
