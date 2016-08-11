ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

### The target here is to make axial projection of the mesh so it would be easier to compare with other works

### Extracting mesh elements
using JLD
@load "compute/frames.jld"

function axialProjection(points)
    z = []
    r = []

    for i in 1:size(p,2)
        pos = p[:,i]
        push!(z,pos[1])
        push!(r,sqrt(pos[2]^2 + pos[3]^2))
    end

    s = sortperm(z)
    rs = s[end:-1:1]

    return [z[s];z[rs]], [r[s];-r[rs]]
end


#include("../utils.jl")
#t = data["t"]
#elpar = [ellipsoid_parameters(frames[i][1]) for i in 1:size(t,1)]


### Pcking out now a single frame

#n = 20
import PyPlot; const plt = PyPlot
plt.figure()

for n in [1,20,40,90]

    t = memory[n][1]
    p = memory[n][2]

    z = []
    r = []

    for i in 1:size(p,2)
        pos = p[:,i]
        push!(z,pos[1])
        push!(r,sqrt(pos[2]^2 + pos[3]^2))
    end

    println("max z is $(maximum(z)) and max r is $(maximum(r))")
    ### And now simple plotting

    s = sortperm(z)
    # plt.plot(z[s],r[s],"b-",lw=2)
    # plt.plot(z[s],-r[s],"b-",lw=2)
    plt.plot([z[s];z[s[end:-1:1]]],[r[s];-r[s[end:-1:1]]],lw=2,label="$(round(t,2))")

end
    
plt.axis("equal")
plt.legend()
plt.savefig("figures/relax.pdf")
#plt.xlim(-2,+2)
#plt.ylim(-2,+2)
