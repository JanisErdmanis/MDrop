#include(".juliarc.jl")

using JLD
using SurfaceGeometry

using Distributed
@everywhere include("field.jl")
include("velocity.jl")

function getabc(points)

    ar = 0
    al = 0
    br = 0
    bl = 0
    cr = 0
    cl = 0
    
    for xkey in 1:size(points,2)
        x = points[:,xkey]
        x[1]<al && (al=x[1])
        x[1]>ar && (ar=x[1])
        x[2]<bl && (bl=x[2])
        x[2]>br && (br=x[2])
        x[3]<cl && (cl=x[3])
        x[3]>cr && (cr=x[3])
    end

    a = (ar - al)/2
    b = (br - bl)/2
    c = (cr - cl)/2

    return a,b,c
end

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

    ### For testing
    normals = Array{Float64}(undef,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    s = 0

    for xkey in 1:size(points,2)
        s += dot(H0/2*[psix[xkey],psiy[xkey],0],normals[:,xkey]) * vareas[xkey]
    end

    Es = gammap * Area
    Em = 1/8/pi * (1 - mup) * s

    ### Here I could also do the normalisation of it    
    return Es+Em
end

function PointPerturbation(scale)
    DR = rand()*scale/50
    phi = rand()*2*pi
    theta = rand()*pi

    x = DR*sin(theta)*cos(phi)
    y = DR*sin(theta)*sin(phi)
    z = DR*cos(theta)
    return [x,y,z]
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

    ### Stretch a bit before every quasystep
    # a,b,c = getabc(points)
    # if a/b<1.03
    #     for j in 1:size(points,2)
    #         x,y,z = points[:,j]
    #         points[:,j] = [1.03*x,y,z]
    #     end
    # end

    a,b,c = getabc(points)
    if abs(a-b)<sc
        for j in 1:size(points,2)
            x,y,z = points[:,j]
            points[:,j] = [(1+sc/a)*x,y,z]
        end
    end
    ### For more complex perturbation
    # for j in 1:size(points,2)
    #     x,y,z = points[:,j]
    #     r = sqrt(x^2 + y^2)
    #     theta = atan2(z,r) 
    #     #phi = atan2(y,x)
    #     nx,ny,nz = normals[:,j]
    #     phi = atan2(ny,nx)
    #     Delta = scale/10*(cos(2*phi) + cos(3*phi))*cos(theta)^2
    #     points[:,j] += Delta*normals[:,j]
    # end

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

        fieldx = @spawn surfacefield(points,faces,normals,mup,H0i*[1,0,0])
        fieldy = @spawn surfacefield(points,faces,normals,mup,H0i*[0,1,0])

        psix,Htx,Hnx = fetch(fieldx)
        psiy,Hty,Hny = fetch(fieldy)
        
        tensorn = mup*(mup-1)/8/pi/2 * (Hnx.^2 + Hny.^2) + (mup-1)/8/pi/2 * (Htx.^2 + Hty.^2)
        vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

        ### Calculating variables for this step
        #pDx = xp - xi
        Ei = DropEnergy(points,faces,normals,psix,psiy,H0i)
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

        # if Ei>Ep
        #     FluctatingEnergy = true 
        #     Equilibrium = true
        #     break
        # end
                
        #if ti!=tp && xp-xi>0 && h/log(vp/vi)*vi < scale/100  # 1000
        # abs((taui-taup)/taui)<1e-5 &&
        #if ti!=tp && h/log(vp/vi)*vi < sc/10 && vi/v0max<0.01  # 1000
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
        actualdt,points,faces = improvemeshcol(pointsp,faces,points,par)
        

        # when xp-xi becomes positive
        # what if oscillates??
        # if pDx<0 && xp-xi>0
        #     for j in 1:size(points,2)
        #         points[:,j] += PointPerturbation(scale)
        #     end
        # end

        # ### Can be commented out with ease
        # if !(maximum(abs(vn))*(ti-tp) < scale && ti!=tp && xp-xi>0)
        #     actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
        # end
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

