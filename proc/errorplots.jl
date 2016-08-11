isdefined(:method) || (method=:potentialsimplereg)
isdefined(:session) || (session="")

using SurfaceGeometry
using Storage
#include("utils.jl")
include("Calculation/field.jl")

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

a,b,c = 1,1,1
p1,t1 = subdivision(EllipsoidMeshLoad(a,b,c,0.35)...,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
p2, t2 = subdivision(p1,t1,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
meshes = [EllipsoidMeshLoad(a,b,c,0.2),EllipsoidMeshLoad(a,b,c,0.1),EllipsoidMeshLoad(a,b,c,0.35),(p2,t2)]

Htheor = EllipsoidField(a,b,c,10,[1,0,0])

#(points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

import PyPlot; const plt = PyPlot
#plt.rc("text",usetex=true)
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
plt.rc("legend",fontsize=8)

plt.grid(which="minor",color="grey",linestyle="")
plt.grid(which="major",color="black",linestyle="")

width = 3.487 * 1.25
height = width / 1.618

fig = plt.figure(figsize=(width,height))

if method in [:potentialsimplereg,:potentialsimple,:potentialpozikridis,:potentialgaussiantrapezodial]

    for (p,t) in meshes

        N = size(p,2)
        
        if method==:potentialsimplereg
            psi = PotentialSimple(p,t,10,[1,0,0])  ### This part could change
        elseif method==:potentialsimple
            psi = PotentialSimple(p,t,10,[1,0,0]; regularize=false)
        elseif method==:potentialpozikridis
            #rp, rfaces = subdivision(p,t; method=:paraboloid)
            rp, rfaces = subdivision(p,t,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
            psi = PotentialGaussianPozikridis(rp,t,rfaces,10,[1,0,0];NP=3)
            N = size(rp,2)
        elseif method==:potentialgaussiantrapezodial
            rp, rfaces = subdivision(p,t,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
            psi = PotentialGaussianTrapezodial(rp,t,rfaces,10,[1,0,0];NP=3)
            N = size(rp,2)
        end
        
        r = Float64[]
        for xkey in 1:size(p,2)
            psit = dot(Htheor,p[:,xkey])
            push!(r,abs((psi[xkey] - psit)/psit))
        end
        r = r*100

        
        plt.hlines(N,0,mean(abs(r)))    
        plt.plot(abs(r),ones(r)*N,"+")
    end

elseif method in [:normalfield,:normalfieldreg,:normalfieldrecalculated,:normalfieldgaussian,:normalfieldtrapezodial,:normalfieldcurrent]
    for (p,t) in meshes

        normals = Array(Float64,size(p)...)
        NormalVectors!(normals,p,t,i->FaceVRing(i,t))

        if method==:normalfield
            Hn = NormalField(p,t,10,[1,0,0]; regularize=false)  ### This part could change
        elseif method==:normalfieldreg
            Hn = NormalField(p,t,10,[1,0,0]; regularize=true)  ### This part could change
        elseif method==:normalfieldrecalculated
            psit = Array(Float64,size(p,2))
            for xkey in 1:size(p,2)
                nx = normals[:,xkey]
                psit[xkey] = dot(Htheor,p[:,xkey])
            end
            Hn = NormalFieldRecalculated(p,t,normals,psit)
        elseif method==:normalfieldtrapezodial
            rp, rfaces = subdivision(p,t,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
            #rp, rfaces = subdivision(p,t; method=:paraboloid)
            Hn = NormalFieldTrapezodial(rp,t,rfaces,10,[1,0,0];NP=3)
        elseif method==:normalfieldcurrent
            Ht = Array(Float64,3,size(p,2))
            Htheor = EllipsoidField(a,b,c,10,[1,0,0])
            for xkey in 1:size(p,2)
                Ht[:,xkey] = Htheor
            end

            Hn = NormalFieldCurrent(p,t,Ht,10,[1,0,0])
        end
        
        r = Float64[]
        for xkey in 1:size(p,2)
            Hnt = dot(Htheor,normals[:,xkey])
            push!(r,abs((Hn[xkey] - Hnt)/Hnt))
        end
        r = r*100

        N = size(p,2)
        plt.hlines(N,0,mean(abs(r)))    
        plt.plot(abs(r),ones(r)*N,"+")
    end

elseif method in [:tangentialfieldpotential,:normalfieldpotential]
    for (p,t) in meshes

        normals = Array(Float64,size(p)...)
        NormalVectors!(normals,p,t,i->FaceVRing(i,t))

        psi = PotentialSimple(p,t,10,[1,0,0])

        if method==:tangentialfieldpotential
            H = HField(p,t,psi)
            
            r = Float64[]
            for xkey in 1:size(p,2)
                nx = normals[:,xkey]
                P = eye(3) - nx*nx'
                push!(r,norm(P*(H[:,xkey] - Htheor))/norm(P*Htheor))
            end
            r = r*100
            
        elseif method==:normalfieldpotential

            H = HField(p,t,psi)
            
            r = Float64[]
            for xkey in 1:size(p,2)
                nx = normals[:,xkey]
                push!(r,abs(dot(nx,H[:,xkey]-Htheor)/dot(nx,Htheor)))
            end
            r = r*100

        end
            
        N = size(p,2)
        plt.hlines(N,0,mean(abs(r)))    
        plt.plot(abs(r),ones(r)*N,"+")
    end
end
    
plt.xlabel("Releative error \%")
plt.ylabel("Number of verticies")
plt.ylim(0,3000)

# yscale("log")
# xscale("log")
plt.xlim(0,10)
plt.grid("off")

plt.tight_layout()
fig[:savefig]("figures/$(string(method))$session.pdf")
plt.close("all")
