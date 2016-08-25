using JLD
using SurfaceGeometry

include("modules/field.jl")
include("modules/velocity.jl")

#bname = "/SlowFieldEiler/"

function DropEnergy(points,faces)

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

    Es = gammap * Area

    ### Here I could also do the normalisation of it    
    return Es
end

# if !isdir(datadir*bname)
#     mkdir(datadir*bname)
# end

memory = []
E = []

if con==true
    info("Continuing from last simulation")
    if !isdir(outdir) || isempty(outdir)
        error("No previous simulation found")
    end
    # 107;Bm=25.jld ### when viewing find a point or semicolon
    outfiles = readdir(outdir)

    import Base.isless
    function Base.isless(x::ASCIIString,y::ASCIIString)
        nx = parse(Int,x[1:length(x)-4])
        ny = parse(Int,y[1:length(y)-4])
        return nx < ny
    end

    sort!(outfiles)
    last = outfiles[end]
    i = parse(Int,last[1:length(last)-4]) 
    data = load("$outdir/$last")["memory"][end]
    ti,points,faces = data[1],data[2],data[3]
else
    info("Starting fresh simulation")
    run(`rm -rf $outdir`)
    mkdir(outdir)
    ti = 0
    i = 1
end

push!(memory,(ti,points,faces))
volume0 = volume(points,faces)

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    Ei = DropEnergy(points,faces)
    println("E = $Ei")
    rV = volume(points,faces)/volume0
    println("V/V0 = $rV")
    push!(E,Ei)

    if mod(i,5)==0
        storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
        save("$outdir/$i.jld","memory",storage)
        memory = []
        E = []
    end
    
    tensorn = zeros(points)
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    points1 = Array(Float64,size(points)...)
    for j in 1:size(points,2)
        points1[:,j] = points[:,j] + normals[:,j]*vn[j]*h/2
    end

    ### Maybe it works if I skip recalculation of normal vectors at this point!
    tensorn1 = zeros(points)
    vn1 = InterfaceSpeedZinchenko(points1,faces,tensorn1,etap,gammap)

    oldpoints = copy(points)
    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn1[j]*h
    end
    
    ### Can be commented out with ease
    actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
    ti += h
    i += 1
    push!(memory,(ti,points,faces))
    info("Step $i has finished")
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
save("$outdir/$i.jld","memory",storage)
