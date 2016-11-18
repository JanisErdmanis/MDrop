#ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

#using SurfaceGeometry
#using Storage

include("../modules/field.jl")

using JLD
#include("../field.jl")

a,b,c = 2,1,1
@load "meshes/211-0.2matlab.jld"

# (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

mup = 10

# a,b,c = 2,1/4,1/4
# (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

# a,b,c = 2,1/6,1/6
# (points,faces)=EllipsoidMeshLoad(a,b,c,0.05)

# a,b,c = 1/4,1,1
# points,faces = EllipsoidMeshLoad(a,b,c,0.15)

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

#potential_r += abs((potential[xkey] - dot(Hvec_t,p0[:,xkey]))/dot(Hvec_t,p0[:,xkey]))/N


#points, faces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)

# psi = PotentialSimple(points,faces,10,[1,0,0])
#psi = PotentialSimpleRegularized(points,faces,10,[1,0,0])

### Does error decreases?


#points, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
#@time psi = PotentialGaussian(points,faces,rfaces,10,[1,0,0];NP=3)
#psi = PotentialGaussianTrapezodial(points,faces,rfaces,10,[1,0,0];NP=3)
#psi = PotentialGaussianPozikridis(points,faces,rfaces,10,[1,0,0];NP=3)

# psi,A = PotentialCurved(points,faces,rfaces,10,[1,0,0])


# Htheor = EllipsoidField(a,b,c,mup,[1,0,0])
# for xkey in 1:25:maximum(faces)
#     psit = dot(Htheor,points[:,xkey])
#     r = abs((psi[xkey] - psit)/psit)
#     println("psi is $(round(psit,4)) and computed $(round(psi[xkey],4)) and releative error $r")
# end


normals = Array(Float64,size(points)...)
NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

for xkey in 1:size(points,2)
    x,y,z = points[:,xkey]
    gradf = [x/a^2,y/b^2,z/c^2]
    nx = gradf/norm(gradf)
    normals[:,xkey] = nx
end


psit = Array(Float64,size(points,2))
Ht = Array(Float64,3,size(points,2))
Htheor = EllipsoidField(a,b,c,10,[1,0,0])
for xkey in 1:size(points,2)
    nx = normals[:,xkey]
    psit[xkey] = dot(Htheor,points[:,xkey])
    Ht[:,xkey] = Htheor
end

#@time Hn = NormalFieldRecalculated(points,faces,normals,psi)

#Hn = NormalFieldDomain(points,faces,psit,10,[1,0,0];eps=0.000001)
#Hn = NormalFieldHypersingular(points,faces,psi,Ht)


#points, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
#points, rfaces = subdivision(points,faces; method=:paraboloid)


### EXISTING METHOD
# @time psi = PotentialSimple(points,faces,mup,[1,0,0])
# H = HField(points,faces,psit)
# @time Hn = NormalFieldCurrent(points,faces,Ht,mup,[1,0,0],normals=normals)

### Traditional normal field way
#@time Hn = NormalField(points,faces,mup,[1,0,0],regularize=false,normals=normals)
### Regularized way
@time Hn = NormalField(points,faces,mup,[1,0,0],regularize=true,normals=normals)

#Hn = NormalFieldCurrentRegularized(points,faces,H,10,[1,0,0])
#Hn = NormalFieldCurrentGaussian(points,faces,H,10,[1,0,0],NP=16)

#@time Hn = NormalField(points,faces,mup,[1,0,0];regularize=false)
#@time Hn = NormalFieldGaussian(points,faces,rfaces,normals,10,[1,0,0];NP=16)

#@time Hn = NormalFieldTrapezodial(points,faces,rfaces,10,[1,0,0];NP=3)


normals = Array(Float64,size(points)...)
NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

rarr = Float64[]

Htheor = EllipsoidField(a,b,c,mup,[1,0,0])
for xkey in 1:maximum(faces)
    Hnt = dot(Htheor,normals[:,xkey])
    r = abs((Hn[xkey] - Hnt)/Hnt)
    term = dot([1,0,0],normals[:,xkey])/mup
    #importance = abs((Hn[xkey] - term)/term)
    importance = abs((Hnt - term)/term)
    #integral error
    # integrcal = Hn[xkey] - term
    # integralther = Ht - term
    bigerr = abs((Hn[xkey] - Hnt)/(Hnt - term))

    push!(rarr,r)
    #println("Hn is $(round(Hnt,3)) and computed $(round(Hn[xkey],3)) and r = $(round(r,3)) importance = $(round(importance,3)) ri = $(round(bigerr,3))")
end

println("Average relative error is $(mean(rarr)*100) %")

# ########## Some testing

### All field test 

# @time psi = PotentialSimple(points,faces,10,[1,0,0])
# @time rpoints, rfaces = subdivision(points,faces; method=:paraboloid)
# @time Hn = NormalFieldTrapezodial(rpoints,faces,rfaces,10,[1,0,0];NP=3)
# @time H = HField(points,faces,psi,Hn)

# Htheor = EllipsoidField(a,b,c,10,[1,0,0])
# for xkey in 1:25:maximum(faces)
#     r = norm(H[:,xkey] - Htheor)/norm(Htheor)
#     println("r = $r")
# end
