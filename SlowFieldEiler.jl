include("defaultpar.jl")
### This is the place where one adds his own set of parameters
include("modules/field.jl")
include("modules/velocity.jl")

if !isdir(datadir*"/SlowFieldEiler")
    mkdir(datadir*"/SlowFieldEiler")
end

if isdir(datadir*"/SlowFieldEiler/"*outdir)
    info("Continuing from last simulation")
    ### Now for testing purposes

    # Reinitialise step number i
    # Gets the latest 

    rm(datadir*"/SlowFieldEiler/"*outdir)
    mkdir(datadir*"/SlowFieldEiler/"*outdir)
else
    mkdir(datadir*"/SlowFieldEiler/"*outdir)
end

memory = []
E = []
ti = 0
i = 0
push!(memory,(t,points,faces))

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    # psi = PotentialSimple(points,faces,mup,Htime)
    # H = HField(points,faces,psi)
    # Ht = Array(Float64,size(points,2))
    # for xkey in 1:size(points,2)
    #     nx = normals[:,xkey]
    #     P = eye(3) - nx*nx'
    #     Ht[xkey] = norm(P*H[:,xkey])
    # end

    psi,Ht,Hn = surfacefield(points,faces,normals,mup,H0*[1.,0,0])

    ### At this point I can calculate energy
    push!(E,1)

    if mod(i,5)==0
        storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
        save(datadir*"/SlowFieldEiler/"*outdir*"/$i.jld","memory",storage)
        memory = []
        E = []
    end
    
    tensorn = mup*(mup-1)/8/pi * Hn.^2 + (mup-1)/8/pi * Ht.^2
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    ti += h
    i += 1

    push!(memory,(t,points,faces))
    info("Step $i has finished")
    
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
save(datadir*"/SlowFieldEiler/"*outdir*"/$i.jld","memory",storage)

