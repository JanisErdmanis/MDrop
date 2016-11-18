using JLD
using SurfaceGeometry

@everywhere include("field.jl")
include("velocity.jl")

#bname = "/SlowFieldEiler/"

function DropEnergy(points,faces,normals,psix,psiy,H0)

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

    ### Here I could also do the normalisation of it    
    return Es+Em
end

# if !isdir(datadir*bname)
#     mkdir(datadir*bname)
# end

memory = []
E = []
tau = []

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

### Some event variables
vn = zeros(points)
tp = 0
vi = 0
xi = 0
taui = 0
Ep = Inf
oldpoints = copy(points)
FluctatingEnergy = false

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    fieldx = @spawn surfacefield(points,faces,normals,mup,H0*[1,0,0])
    fieldy = @spawn surfacefield(points,faces,normals,mup,H0*[0,1,0])

    psix,Htx,Hnx = fetch(fieldx)
    psiy,Hty,Hny = fetch(fieldy)

    Ei = DropEnergy(points,faces,normals,psix,psiy,H0)
    println("E = $Ei")
    rV = volume(points,faces)/volume0
    println("V/V0 = $rV")
    push!(E,Ei)

    tensorn = mup*(mup-1)/8/pi/2 * (Hnx.^2 + Hny.^2) + (mup-1)/8/pi/2 * (Htx.^2 + Hty.^2)
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    oldpoints = copy(points)
    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    ti += h
    i += 1

    ### When to stop?

    xp = xi
    #pDx = xp - xi
    vp = vi
    vi = maximum(abs(vn))
    xi = vi*(ti-tp)
    taup = taui
    taui = h/log(vp/vi)
    push!(tau,taui)

    if Ei>Ep
        FluctatingEnergy = true
    end

    println("$xi < $scale; xi-xp=$(xp-xi)")
    println("$(h/log(vp/vi)*vi) < $scale")
    println("$(abs(taui/taup))")

    if ti!=tp && xp-xi>0 && h/log(vp/vi)*vi < scale/100  # 1000
        break
    end

    if mod(i,5)==0
        storage = [tuple(memory[i]...,E[i],tau[i]) for i in 1:length(memory)]
        save("$outdir/$i.jld","memory",storage)
        memory = []
        E = []
        tau = []
    end
    
    ### Can be commented out with ease
    actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
    push!(memory,(ti,points,faces))
    info("Step $i has finished")
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
push!(tau,NaN)
storage = [tuple(memory[i]...,E[i],tau[i]) for i in 1:length(memory)]
save("$outdir/$i.jld","memory",storage)
