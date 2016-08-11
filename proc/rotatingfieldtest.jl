ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

alpha = 1/9

function g(K)
    eps = sqrt(K^2-1)
    return  pi*( eps*(-1 + 2*alpha*eps^2) + K^2*atan(eps))^2 * (2*eps*K*(1 + 2*K^2) + (1 - 4*K^2)*log((K+eps)/(K-eps))) / (K^(7/3)*eps^4*(-3*eps + (2 + K^2)*atan(eps)) )
end

using PyPlot; const plt = PyPlot

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
plt.rc("legend",fontsize=8)

width = 3.487 
height = width / 1.618
fig = plt.figure(figsize=(width,height))

K = 1.001:0.1:30
plt.plot(map(g,K),1./collect(K))
plt.xlabel(L"$ H_0^2 R_0/\gamma$")
plt.ylabel(L"$c/a$")
plt.legend()


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

session = "rotatingfield"
calculate = false
viewmeshactive = true  ### avoiding of plotting

cases = [57,58,65,66,61,64,55,54]
lastelement = [-1,-1,-1,-1,30,18,54,294]

Bm = []
D = []
Derr = []

case = 1
for (case,n) in zip(cases,lastelement)
    include("diaryrot.jl")

    if n==-1
        n = length(memory)
    end

    push!(Bm,norm(H0)^2 * 1/gammap)
    
    Di = deformability(memory[n][2],:z)
    h = tmax/N
    Derri = (Di - deformability(memory[n-1][2],:z))/h * tmax
    
    push!(Derr,Derri)
    push!(D,Di)
end

ab = (D + 1)./(1 - D)
Dab = (1 + ab).*Derr ./(1 - D)
#plt.plot(Bm,ab,".")
plt.errorbar(Bm[1:4],ab[1:4],yerr=Dab[1:4],fmt=".g")
plt.plot(Bm[5:end],ab[5:end],".r")
plt.xlim(0,100)
plt.ylim(0,1)

plt.grid("off")
plt.tight_layout()
fig[:savefig]("figures/rotatingfield.pdf")
plt.close("all")
