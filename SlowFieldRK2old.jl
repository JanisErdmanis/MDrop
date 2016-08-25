include("defaultpar.jl")
include("modules/field.jl")
include("modules/velocity.jl")

bname = "/SlowFieldRK2/"

function DropEnergy(points,faces,normals,psi,H0)

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
        #s += dot(H0/2*[psix[xkey],psiy[xkey],0],normals[:,xkey]) * vareas[xkey]
        s += psi[xkey]*dot(H0,normals[:,xkey]) * vareas[xkey]
    end

    Es = gammap * Area
    Em = 1/8/pi * (1 - mup) * s

    ### Here I could also do the normalisation of it    
    return Es+Em
end

if !isdir(datadir*bname)
    mkdir(datadir*bname)
end

if isdir(datadir*bname*outdir)
    info("Continuing from last simulation")
    ### Now for testing purposes

    # Reinitialise step number i
    # Gets the latest 

    run(`rm -rf $(datadir*bname*outdir)`)
    mkdir(datadir*bname*outdir)
else
    mkdir(datadir*bname*outdir)
end

memory = []
E = []
ti = 0
i = 1
push!(memory,(ti,points,faces))

volume0 = volume(points,faces)

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    psi,Ht,Hn = surfacefield(points,faces,normals,mup,H0*[1.,0,0])

    ### At this point I can calculate energy
    Ei = DropEnergy(points,faces,normals,psi,H0*[1.,0,0])
    println("E = $Ei")
    rV = volume(points,faces)/volume0
    println("V/V0 = $rV")
    push!(E,Ei)

    if mod(i,5)==0
        storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
        save(datadir*bname*outdir*"/$i.jld","memory",storage)
        memory = []
        E = []
    end
    
    tensorn = mup*(mup-1)/8/pi * Hn.^2 + (mup-1)/8/pi * Ht.^2
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    points2 = Array(Float64,size(points)...)
    for j in 1:size(points,2)
        points2[:,j] = points[:,j] + normals[:,j]*vn[j]*h/2
    end

    ### Maybe it works if I skip recalculation of normal vectors at this point!
    psi2,Ht2,Hn2 = surfacefield(points2,faces,normals,mup,H0*[1.,0,0])
    tensorn2 = mup*(mup-1)/8/pi * Hn2.^2 + (mup-1)/8/pi * Ht2.^2
    vn2 = InterfaceSpeedZinchenko(points2,faces,tensorn2,etap,gammap)
        
    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn2[j]*h
    end

    push!(memory,(ti,points,faces))
    info("Step $i has finished")
    
    ti += h
    i += 1
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
save(datadir*bname*outdir*"/$i.jld","memory",storage)

