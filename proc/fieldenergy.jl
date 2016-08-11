### Loads a mesh and calculates a energy for the magnetic liquid drop
ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using SurfaceGeometry

include("Calculation/field.jl")

### Calculation of psi

function RotatingFieldEnergy(points,faces,psix,psiy,mup,gammap,H0)

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

    # psix = PotentialSimple(points,faces,mup,H0*[1,0,0],regularize=true)
    # psiy = PotentialSimple(points,faces,mup,H0*[0,1,0],regularize=true)

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


############### Testing ##################

# using Storage

# a,b,c = 2,1/4,1/4
# (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

# #points,faces = subdivision(points,faces,x->x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)

# mup = 10
# H0 = 10
# gammap = 1

# Es, Em =  RotatingFieldEnergy(points,faces,mup,gammap,H0)
# Emt = TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)
