### A bare way of making calculation
### Because of simplicity here without passive stabilisation

isdefined(:adaptive) || (adaptive=false)
isdefined(:activestabilisation) || (activestabilisation=false)

session = "calcpar"
fname = "test"

using SurfaceGeometry


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


# interface - (t,points,faces,[energy])
memory = []
t = 0.0
h = 0.1
i = 0
#par =

info("Simulation started")
push!(memory,(t,points,faces))

if !isdir(dire*"/"*session)
    mkdir(dire*"/"*session)
end
save("$dire/$session/$fname.jld","memory",memory)

info("Integrator initialised")

Energy = Inf
eps = 1e-6
#success = true

### For now we are going to assume simplicity
### A simple Euler method for testing
while true #success==true
    i +=1

    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    ### Needs a change for surfacefield function
    psix,Htx,Hnx = surfacefield(points,faces,normals,H0*[1,0,0]) 
    psiy,Hty,Hny = surfacefield(points,faces,normals,H0*[0,1,0])

    Es, Em = RotatingFieldEnergy(points,faces,psix,psiy,mup,gammap,H0)
    tensorn = mup*(mup-1)/8/pi/2 * (Hnx.^2 + Hny.^2) + (mup-1)/8/pi/2 * (Htx.^2 + Hty.^2)
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    ### Finishing with integration
    if (Energy-(Em+Es))<0 || maximum(abs(vn)) < eps
        break
    else
        Energy = Em + Es
    end
    
    v = Array(Float64,size(points)...)
    for i in 1:size(points,2)
        v[:,i] = vn[i]*normals[:,i]
    end

    actualdt,points,faces = improvemeshcol(points,faces,points + h*v, par)    
    push!(memory,(t,points,faces))

    if mod(i,5)==0
        save("$dire/$session/$fname.jld","memory",memory)
    end

    info("    E = $(Em + Es)")
    info("Step $i have been finished.")
end
