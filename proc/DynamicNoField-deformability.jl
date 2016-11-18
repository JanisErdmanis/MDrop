ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using JLD
isdefined(:session) || (session="relaxation")


using SurfaceGeometry

using PyCall
import PyPlot; const plt = PyPlot
using LaTeXStrings


#plt.rc("font", family="serif", serif="Times")

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)

### Plotting
width = 3.487
height = width / 1.618
fig = plt.figure(figsize=(width,height))


# fig, ax = plt.subplots()
# fig[:subplots_adjust](left=.19, bottom=.19, right=0.97, top=.97)

#dire = Pkg.dir("Storage","simulations",session)
#PrintData(session="relaxation")

calculate = false
viewmeshactive = true
cases = [3,4]

case = nothing
for case in cases
    #st = load(dire*"/"*fi)
    include("diaryrelax.jl")
    
    #a,b,c = st["a"], st["b"], st["c"]
    #etap = st["etap"]
    #gammap = st["gammap"]

    tau = etap*(a*b*c)^(1/3)/gammap
    
    #memory = st["memory"]
    
    t = Float64[memory[i][1] for i in 1:length(memory)]
    D = Float64[]

    #println("The length is $(length(t))")
    
    for i in 1:length(t)
        p = memory[i][2]
        z = []
        r = []

        for i in 1:size(p,2)
            pos = p[:,i]
            push!(z,pos[1])
            push!(r,sqrt(pos[2]^2 + pos[3]^2))
        end

        D_ = (maximum(z) - maximum(r))/(maximum(r) + maximum(z))
        push!(D,D_)
    end

    plt.plot(t/tau,log(abs(D)))
end


#\mathrm{\alpha}
plt.ylabel(L"\ln D")
plt.xlabel(L"t \cdot \eta R_0 / \gamma")
#plt.ylabel("Logarithm of deformability")
#plt.xlabel("Dimensionless time")
#plt.xlim(0,maximum(t))

#fig[:set_size_inches](width,height)
plt.grid("off")
plt.tight_layout()

fig[:savefig]("figures/$session.pdf")



# fig = figure()

# plot(t,log(eccentricity),"k.-",lw=3)
# xlabel(L"t")
# ylabel(L"\log(e)")
# xlim(0,20)
