using JLD
using SurfaceGeometry
using Distributed

@everywhere include("field.jl")
include("velocity.jl")

#bname = "/SlowFieldEiler/"

function DropEnergy(points,faces,normals,psix,psiy,H0)

    vareas = zero(points)
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
    @info "Starting fresh simulation"
    run(`rm -rf $outdir`)
    mkdir(outdir)
    ti = 0
    i = 1
end

volume0 = volume(points,faces)
sc = 0.01*(3*volume0/4/pi)^(1/3) ### charectaristic scale

@info "Proceding with simulation at Bm=$Bm"
H0 = sqrt(Bm*gammap/(volume(points,faces)*3/4/pi)^(1/3))

memory = []
vn = zero(points)
tp = 0
ip = i
vi = 0
xi = 0
taui = 0
Ei = Inf
FluctatingEnergy = false
Equilibrium = false

# WTF issue
points2 = copy(points)
faces2 = copy(faces)

while true
    @info "Starting with step $i"

    points = points2
    faces = faces2
    
    Ep = Ei
    xp = xi
    vp = vi
    taup = taui
    pointsp = copy(points)
    
    normals = Array{Float64}(undef,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    fieldx = @spawn surfacefield(points,faces,normals,mup,H0*[1,0,0])
    fieldy = @spawn surfacefield(points,faces,normals,mup,H0*[0,1,0])

    psix,Htx,Hnx = fetch(fieldx)
    psiy,Hty,Hny = fetch(fieldy)
    
    tensorn = mup*(mup-1)/8/pi/2 * (Hnx.^2 + Hny.^2) + (mup-1)/8/pi/2 * (Htx.^2 + Hty.^2)
    global vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    ### Calculating variables for this step
    #pDx = xp - xi
    global Ei = DropEnergy(points,faces,normals,psix,psiy,H0)
    rV = volume(points,faces)/volume0
    global vi = maximum(abs.(vn))
    global xi = vi*(ti-tp)
    global taui = h/log(vp/vi)

    if i==ip
        global v0max = vi
    end
    
    println("E = $Ei")
    println("$(h/log(vp/vi)*vi) < $(sc)")
    println("vi/v0max = $(vi/v0max)")
    push!(memory,(ti,copy(points),copy(faces),Ei,taui))

    if mod(i,5)==0
        save("$outdir/$i.jld","memory",memory)
        global memory = []
    end

    ### This is more cosmetic one
    if rV==NaN || abs(rV - 1)>0.5
#        break
    end

    if (i-ip)>1000
 #       break
    end

    if ti!=tp && h/log(vp/vi)*vi < sc && vi/v0max<0.01  # 1000
        global Equilibrium = true
#        break
    end
    
    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    global ti += h
    global i += 1

    # ElTopo magic.
    actualdt,points,faces = improvemeshcol(pointsp,faces,points,par)
    #points,faces = improvemesh(points,faces,par)

end

save("$outdir/$i.jld","memory",memory)
if Equilibrium==false
    @info "Simulation did not achieve equilibrium"
end

