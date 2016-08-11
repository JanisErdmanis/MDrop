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

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
plt.rc("legend",fontsize=8)

width = 3.487 
height = width / 1.618
fig = plt.figure(figsize=(width,height))
#plt.figure()

session = "rotatingfield"
calculate = false
viewmeshactive = true  ### avoiding of plotting

# cases = [57,58,65,66,61,64,55,54]
# lastelement = [-1,-1,-1,-1,30,18,54,294]

case = 55
n = 54


include("diaryrot.jl")

if n==-1
    n = length(memory)
end

faces = memory[n][3]
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

#s = sortperm(z)
# plt.plot(z[s],r[s],"b-",lw=2)
# plt.plot(z[s],-r[s],"b-",lw=2)

Bm = H0^2 * (a*b*c)^(1/3) / gammap

#plt.plot([z[s];z[s[end:-1:1]]],[r[s];-r[s[end:-1:1]]],".",lw=2,label="$Bm and $case")

plt.plot(z,r,".",lw=2)

####

a,c = maximum(r), maximum(z)
a = a*0.98

phi = collect(0:0.1:pi)
x = c*cos(phi)
y = a*sin(phi)
plt.plot(x,y,lw=2,label=L"z^2/c^2 + r^2/a^2 = 1")

phi = collect(0:0.1:pi)
x = c*sqrt(abs(cos(phi)))
y = a*sqrt(abs(sin(phi)))
plt.plot(x,y,lw=2,label=L"z^4/c^4 + r^4/a^4 = 1")

phi = collect(0:0.1:pi)
x = c*(sqrt(abs(cos(phi))) + abs(cos(phi)))/2
y = a*(sqrt(abs(sin(phi))) + abs(sin(phi)))/2
plt.plot(x,y,lw=2)

# phi = collect(0:0.1:pi)
# x = c*abs(cos(phi)).^(3/2)
# y = a*abs(sin(phi)).^(3/2)
# plt.plot(x,y,lw=2)

plt.xlabel(L"z")
plt.ylabel(L"\sqrt{x^2 + y^2}")
plt.title(L"Bm=35;\mu=10")
plt.axis("equal")
plt.legend(fontsize=7,frameon=false)

ax = plt.gca()

ax[:set_yticks]([0,0.5,1.,1.5])
ax[:set_yticklabels]([0,0.5,1.,1.5])

#plt.legend()
#plt.tight_layout()
plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.20)
plt.savefig("figures/rotatingprojspec.pdf")
#plt.xlim(-2,+2)
#plt.ylim(-2,+2)
plt.close("all")


### calculation of energy

include("fieldenergy.jl")

boxsize = [-2.1 -2 -1.; 2.1 2 1.] #[a a c; a a c] * 1.5

function elmesh2(a,c,stepsize)
    b = a
    
    fdis = """
        function f = fd(p)
            f = p(:,1).^2/$(a^2)+p(:,2).^2/$(b^2)+p(:,3).^2/$(c^2)-1;
        end
    
    """
    p,t = SurfaceMesh(fdis,DistmeshSurfaceMesher(stepsize,boxsize))

    return p,t
end

function elmesh4(a,c,step)

    fdis = """
    function f = fd(p)
        f = (p(:,1).^2 + p(:,2).^2).^2/$(a^4) + p(:,3).^4/$(c^4)-1;
    end
    
    """
    p,t = SurfaceMesh(fdis,DistmeshSurfaceMesher(step,boxsize))

    return p,t
end

isdefined(:calculateenergy) || (calculateenergy=false)
if calculateenergy==true

    Es, Em = RotatingFieldEnergy(points,faces,mup,gammap,H0)
    println("calculation: Es=$Es \t Em=$Em \t Et=$(Em+Es)")
    vol = volume(points,faces)

    points,faces = elmesh2(a,c,0.2)
    volel = volume(points,faces)
    Es, Em = RotatingFieldEnergy(points.*(vol/volel)^(1/3),faces,mup,gammap,H0)
    Emt = TheoreticalRotatingFieldEnergy(a,a,c,mup,H0)
    println("elipsoid: Es=$Es \t Em=$Em \t Et=$(Em+Es) \t Emt=$Emt")

    points,faces = elmesh4(a,c,0.2)
    volq = volume(points,faces)
    Es, Em = RotatingFieldEnergy(points.*(vol/volq)^(1/3),faces,mup,gammap,H0)
    println("higher order: Es=$Es \t Em=$Em \t Et=$(Em+Es)")
end
