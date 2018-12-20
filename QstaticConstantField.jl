#include(".juliarc.jl")

using JLD
using SurfaceGeometry


include("field.jl")
include("velocity.jl")

function DropEnergy(points,faces,normals,psi,H0)

    vareas = zero(points)
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    Area = sum(vareas)

    ### For testing
    normals = Array{Float64}(undef,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    s = 0

    for xkey in 1:size(points,2)
        s += psi[xkey]*dot(H0,normals[:,xkey]) * vareas[xkey]
    end

    Es = gammap * Area
    Em = 1/8/pi * (1 - mup) * s

    ### Here I could also do the normalisation of it    
    return Es+Em
end

if con==true
    error("Not implemented")
    
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

    if typeof(Bm)<:Real
        Bm_ = 1:Bm
    elseif typeof(Bm)<:Range
        Bm_=Bm
    elseif typeof(Bm)<:Array
        Bm_=Bm
    else
        error("Type of $(typeof(Bm)) for Bm not supported.")
    end
end

#
volume0 = volume(points,faces)
sc = 0.01*(3*volume0/4/pi)^(1/3) ### charectaristic scale

for Bmi in Bm_
    @info "Proceding with step at Bm=$Bmi"
    H0i = sqrt(Bmi*gammap/(volume(points,faces)*3/4/pi)^(1/3))

    global memory = []
    vn = zero(points)
    global tp = 0
    global ip = i
    vi = 0
    xi = 0
    taui = 0
    Ei = Inf
    FluctatingEnergy = false
    Equilibrium = false

    ### Need to rescale for volume
    factor = (volume0/volume(points,faces))^(1/3)
    global points = factor*points

    while true
        @info "Starting with step $i"
        
        Ep = Ei
        xp = xi
        vp = vi
        taup = taui
        pointsp = copy(points)
        
        normals = Array{Float64}(undef,size(points)...);
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

        psi,Ht,Hn = surfacefield(points,faces,normals,mup,H0i*[1,0,0])

        tensorn = mup*(mup-1)/8/pi * Hn.^2 + (mup-1)/8/pi * Ht.^2
        vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

        ### Calculating variables for this step
        #pDx = xp - xi
        Ei = DropEnergy(points,faces,normals,psi,H0i*[1,0,0])
        rV = volume(points,faces)/volume0
        vi = maximum(abs.(vn))
        xi = vi*(ti-tp)
        taui = h/log(vp/vi)

        if i==ip
            global v0max = vi
        end
        
        println("E = $Ei")
        println("$(h/log(vp/vi)*vi) < $(sc)")
        println("vi/v0max = $(vi/v0max)")
        push!(memory,(ti,copy(points),copy(faces),Ei,taui))

        ### This is more cosmetic one
        if rV==NaN || abs(rV - 1)>0.5
            break
        end

        if (i-ip)>1000
            break
        end

        if ti!=tp && h/log(vp/vi)*vi < sc && vi/v0max<0.01  # 1000
            Equilibrium = true
            break
        end

        ### Integration must be allowed only afterwards
        for j in 1:size(points,2)
            points[:,j] += normals[:,j]*vn[j]*h
        end
        global ti += h
        global i += 1

        # ElTopo magic
        global actualdt,points,faces = improvemeshcol(pointsp,faces,points,par)
    end

    if FluctatingEnergy==true
        @info "Step terminated since energy did fluctate"
    end

    save("$outdir/$i.jld","memory",memory)

    if Equilibrium==false
        @info "Simulation did not achieve equilibrium"
        #break
    end

end

