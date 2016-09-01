using JLD
using SurfaceGeometry

include("modules/field.jl")
include("modules/velocity.jl")

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

function PointPerturbation(scale)
    DR = rand()*scale/20
    phi = rand()*2*pi
    theta = rand()*pi

    x = DR*sin(theta)*cos(phi)
    y = DR*sin(theta)*sin(phi)
    z = DR*cos(theta)
    return [x,y,z]
end

# if !isdir(datadir*bname)
#     mkdir(datadir*bname)
# end

memory = []
E = []

if con==true
    error("Not implemented")
    
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

    Bmi = 15.
    H0i = sqrt(Bmi*gammap/(volume(points,faces)*3/4/pi)^(1/3))
end

push!(memory,(ti,points,faces))
volume0 = volume(points,faces)

vn = zeros(points)
tp = 0
xp = 0
xi = 0
oldpoints = copy(points)

while true
    normals = Array(Float64,size(points)...);
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    psix,Htx,Hnx = surfacefield(points,faces,normals,mup,H0i*[1,0,0])
    psiy,Hty,Hny = surfacefield(points,faces,normals,mup,H0i*[0,1,0])

    Ei = DropEnergy(points,faces,normals,psix,psiy,H0)
    println("E = $Ei")
    rV = volume(points,faces)/volume0
    println("V/V0 = $rV")
    push!(E,Ei)

    tensorn = mup*(mup-1)/8/pi/2 * (Hnx.^2 + Hny.^2) + (mup-1)/8/pi/2 * (Htx.^2 + Hty.^2)
    vn = InterfaceSpeedZinchenko(points,faces,tensorn,etap,gammap)

    for j in 1:size(points,2)
        points[:,j] += normals[:,j]*vn[j]*h
    end
    ti += h
    i += 1

    xp = xi
    xi = maximum(abs(vn))*(ti-tp)
    println("$xi < $scale; xi-xp=$(xp-xi)")
    if maximum(abs(vn))*(ti-tp) < scale/5 && ti!=tp && xp-xi>0
        storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
        save("$outdir/$i.jld","memory",storage)
        memory = []
        E = []
        if Bmi>Bm
            break
        else
            Bmi += 1.
        end
        H0i = sqrt(Bmi*gammap/(volume(points,faces)*3/4/pi)^(1/3))
        tp = ti
        info("Proceeding with next quasistep")
        xi = 0

        actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
        oldpoints = copy(points)

        for j in 1:size(points,2)
            points[:,j] += PointPerturbation(scale)
        end
    end
    
    # ### Can be commented out with ease
    # if !(maximum(abs(vn))*(ti-tp) < scale && ti!=tp && xp-xi>0)
    #     actualdt,points,faces = improvemeshcol(oldpoints,faces,points,par)
    # end

    push!(memory,(ti,copy(points),faces))
    info("Step $i has finished")
end

### Storage for last steps
### Makes sense only if simulation exits while loop by itself
push!(E,NaN)
storage = [tuple(memory[i]...,E[i]) for i in 1:length(memory)]
save("$outdir/$i.jld","memory",storage)
