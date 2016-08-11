ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

#using JLD
#@load "compute/frames.jld"

function axialProjection(p)
    z = []
    r = []

    for i in 1:size(p,2)
        pos = p[:,i]
        push!(z,pos[3])
        push!(r,sqrt(pos[1]^2 + pos[2]^2))
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

session = "rotatingfield"
calculate = false
viewmeshactive = true  ### avoiding of plotting

# cases = [57,58,65,66,61,64,55,54]
# lastelement = [-1,-1,-1,-1,30,18,54,294]

cases = [66,58,55,54]
lastelement = [-1, -1, 54, 294]

for (case,n) in zip(cases,lastelement)
    include("diaryrot.jl")

    if n==-1
        n = length(memory)
    end
    
    faces = memory[n][1]
    points = memory[n][2]

    z = []
    r = []

    for i in 1:size(points,2)
        pos = points[:,i]
        push!(z,pos[3])
        push!(r,sqrt(pos[1]^2 + pos[2]^2))
    end

    println("max z is $(maximum(z)) and max r is $(maximum(r))")
    ### And now simple plotting

    s = sortperm(z)
    # plt.plot(z[s],r[s],"b-",lw=2)
    # plt.plot(z[s],-r[s],"b-",lw=2)

    Bm = H0^2 * (a*b*c)^(1/3) / gammap
    
    plt.plot([z[s];z[s[end:-1:1]]],[r[s];-r[s[end:-1:1]]],".",lw=2,label="$Bm")

end
    
plt.axis("equal")
plt.legend(title=L"Bm")
plt.savefig("figures/rotatingproj.pdf")
#plt.xlim(-2,+2)
#plt.ylim(-2,+2)
