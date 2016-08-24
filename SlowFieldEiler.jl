using JLD
using SurfaceGeometry

include("defaultpar.jl")
include("modules/field.jl")
include("modules/velocity.jl")

bname = "/SlowFieldEiler/"

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

memory = []
E = []

if con==true
    info("Continuing from last simulation")
    if !isdir(datadir*bname*outdir) || isempty(datadir*bname*outdir)
        error("No previous simulation found")
    end
    # 107;Bm=25.jld ### when viewing find a point or semicolon
    outfiles = readdir(datadir*bname*outdir)

    import Base.isless
    function Base.isless(x::ASCIIString,y::ASCIIString)
        nx = parse(Int,x[1:length(x)-4])
        ny = parse(Int,y[1:length(y)-4])
        return nx < ny
    end

    sort!(outfiles)
    last = outfiles[end]
    i = parse(Int,last[1:length(last)-4]) 
    data = load(datadir*bname*outdir*"/"*last)["memory"][end]
    ti,points,faces = data[1],data[2],data[3]
else
    info("Starting fresh simulation")
    run(`rm -rf $(datadir*bname*outdir)`)
    mkdir(datadir*bname*outdir)
    ti = 0
    i = 1
end

push!(memory,(ti,points,faces))
volume0 = volume(points,faces)

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    psi,Ht,Hn = surfacefield(points,faces,normals,mup,H0*[1.,0,0])

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

    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    ti += h
    i += 1

    push!(memory,(ti,points,faces))
    info("Step $i has finished")
    
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
save(datadir*bname*outdir*"/$i.jld","memory",storage)

