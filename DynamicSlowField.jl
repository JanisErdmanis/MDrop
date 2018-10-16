using JLD
using SurfaceGeometry

include("field.jl")
include("velocity.jl")

#bname = "/SlowFieldEiler/"

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

# if !isdir(datadir*bname)
#     mkdir(datadir*bname)
# end

memory = []
E = []

if con==true
    @info "Continuing from last simulation"
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

#push!(memory,(ti,points,faces))
volume0 = volume(points,faces)
memory = []

while true
    @info "Starting with step $i"
    
    normals = Array{Float64}(undef,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    psi,Ht,Hn = surfacefield(points,faces,normals,mup,H0*[cos(omega*ti),sin(omega*ti),0])

    Ei = DropEnergy(points,faces,normals,psi,H0*[cos(omega*ti),sin(omega*ti),0])
    println("E = $Ei")
    rV = volume(points,faces)/volume0
    println("V/V0 = $rV")

    tensorn = mup*(mup-1)/8/pi * Hn.^2 + (mup-1)/8/pi * Ht.^2
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)
    
    push!(memory,(ti,copy(points),copy(faces),Ei))
    if mod(i,5)==0
        save("$outdir/$i.jld","memory",memory)
        global memory = []
    end

    ### Mesh stabilisation. Passive stabilisation tends to be usefull
    oldpoints = copy(points)
    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    global ti += h
    global i += 1

    ### Can be commented out with ease. ElTopo magic.
    actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
end

save("$outdir/$i.jld","memory",memory)
