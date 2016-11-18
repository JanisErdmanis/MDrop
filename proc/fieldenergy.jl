### Loads a mesh and calculates a energy for the magnetic liquid drop
#ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using SurfaceGeometry

include("../modules/field.jl")

### Calculation of psi

function RotatingFieldEnergy(points,faces,mup,gammap,H0)

    normals = Array(Float64,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    Area = sum(vareas)

    psix = PotentialSimple(points,faces,mup,H0*[1,0,0],regularize=true)
    psiy = PotentialSimple(points,faces,mup,H0*[0,1,0],regularize=true)

    s = 0

    for xkey in 1:size(points,2)
        s += dot(H0/2*[psix[xkey],psiy[xkey],0],normals[:,xkey]) * vareas[xkey]
    end

    Es = gammap * Area
    Em = 1/8/pi * (1 - mup) * s

    return Es,Em
end

#### Theoretical magnetic energy ######

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

function TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Hx = EllipsoidField(a,b,c,mup,[H0,0,0])
    Hy = EllipsoidField(a,b,c,mup,[0,H0,0])

    Emt = 1/8/pi*(1 - mup) * (dot(Hx,H0*[1,0,0]) + dot(Hy,H0*[0,1,0]))/2 * 4/3*pi * a * b * c
    return Emt
end


import Elliptic
function EllipsoidArea(a,b,c)

    if a<b
        a,b = b
    end
    if b<c
        b,c = c,a
    end
    if a<b
        a,b = b,a
    end

    cosphi = c/a
    phi = acos(cosphi)
    sinphi = sqrt(1 - cosphi^2)
    k = sqrt(a^2*(b^2-c^2)/b^2/(a^2-c^2))

    S = 2*pi*c^2 + 2*pi*a*b/sinphi * (sinphi^2 * Elliptic.E(phi,k) + cosphi^2 * Elliptic.F(phi,k))

    return S
end

function TheoreticalDropEnergy(a,b,c,mup,Bm)
    ### I should scale a,b,c inside
    
    gammap = 1.
    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    Emt = TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Etotal = Emt + gammap*EllipsoidArea(a,b,c) # + surface area
end


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

############### Testing ##################

using JLD

#a,b,c = 2,1/4,1/4
#@load "meshes/desa0.1.jld"
#@load "meshes/mu10Bm50.jld"
@load "meshes/sphere0.2.jld"
a,b,c = getabc(points)

#(points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
#points,faces = subdivision(points,faces,x->x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)

mup = 10
#H0 = 100
gammap = 1
Bm = 1
H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))

Es, Em =  RotatingFieldEnergy(points,faces,mup,gammap,H0)
Emt = TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)
Est =  gammap*EllipsoidArea(a,b,c) 

println("Etotal=$(Es+Em) Etotalt=$(Est + Emt)")
