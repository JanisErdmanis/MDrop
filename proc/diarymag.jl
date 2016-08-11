# This file represents a wrapper for a dimensionless units

ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"
using Storage
using SurfaceGeometry

include("Calculation/drivers.jl")
include("Calculation/integrator.jl")

session = "rotatingmag"

# Some default parameters
gammap = 1
etap = 1
mup = 10
N = 100
tmax = 4

elparameters(scale) = Elparameters(
 m_use_fraction = false,
 m_min_edge_length = 0.7*scale,
 m_max_edge_length = 1.5*scale,
 m_max_volume_change = 0.1*scale^3,
 m_min_curvature_multiplier = 1,
 m_max_curvature_multiplier = 1,
 m_merge_proximity_epsilon = 0.5*scale,
 m_proximity_epsilon = 0.00001,
 m_perform_improvement = true, 
 m_collision_safety = false,
 m_min_triangle_angle = 15,
 m_max_triangle_angle = 120,
 m_allow_vertex_movement = true, ### 
 m_use_curvature_when_collapsing = false,
 m_use_curvature_when_splitting = false,
 m_dt = 1
)

zc = nothing
par = nothing

isdefined(:calculate) || (calculate=true)

if case==1

    Bm = 30
    omega = 0

    N = 100
    tmax = 10
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)

    par = elparameters(0.2)
    zc = nothing
    
elseif case==2

    Bm = 30
    omega = 1

    N = 100
    tmax = 10
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau

    par = elparameters(0.2)
    zc = nothing
    
elseif case==3

    Bm = 0

    N = 100
    tmax = 9*2

    a,b,c = 2,1/4,1/4
    points,faces=EllipsoidMeshLoad(a,b,c,0.1)
    
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *=tau

    par = elparameters(0.1)
    zc = nothing

### relaxationtest
elseif case==4

    Bm = 0

    N = 480
    tmax = 9

    a,b,c = 1.1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)
    
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau

    curvatureless = false
    par = nothing
    zc = nothing

elseif case==5

    Bm = 0

    N = 480
    tmax = 9

    a,b,c = 1.1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)
    
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *=tau

    curvatureless = true
    par = nothing
    zc = nothing

elseif case==6

    Bm = 0

    N = 480
    tmax = 9/2

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)
    
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau

    curvatureless = true 
    par = elparameters(0.1)
    zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

### constantfieldtest
elseif case==7

    Bm = 20
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing

elseif case==8

    Bm = 15
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing

elseif case==9

    Bm = 10
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing
    

elseif case==10

    Bm = 5
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    
    par = elparameters(0.1)
    zc = nothing

elseif case==11

    Bm = 24
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing

elseif case==12

    Bm = 2.5
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing

elseif case==13

    Bm = 7.5
    omega = 0

    N = 400
    tmax = 10

    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
        
    par = elparameters(0.1)
    zc = nothing

### rotating field dynamics    
elseif case==14

    Bm = 35
    omega = 0.175

    N = 600
    tmax = 90

    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau
        
    par = elparameters(0.2)
    zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

elseif case==15

    Bm = 35
    omega = 0.5

    N = 600*2
    tmax = 80*2

    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau
        
    par = elparameters(0.2)
    zc = nothing

elseif case==16

    Bm = 35
    omega = 3

    N = 600*2
    tmax = 80*2

    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau
    
    par = elparameters(0.2)
    zc = nothing

elseif case==17

    Bm = 35
    omega = 1

    N = 600*2
    tmax = 80*2

    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau
    
    par = elparameters(0.2)
    zc = nothing

elseif case==18

    Bm = 35
    omega = 0.05

    N = 600*2
    tmax = 80*2

    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau
    
    par = elparameters(0.2)
    zc = nothing
    

# pancake instability    
elseif case==19

    Bm = 35
    omega = Inf

    N = 900
    tmax = 250

    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau

    par = elparameters(0.2)
    zc = nothing

elseif case==20

    Bm = 22.5
    omega = Inf

    N = 900
    tmax = 250

    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau

    par = elparameters(0.2)
    zc = nothing

### histeresis
# a worm
elseif case==21

    Bm = 72
    omega = Inf

    N = 1200
    tmax = 20

    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    points,faces = subdivision(points,faces)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau

    par = elparameters(0.05)
    zc = nothing

# a star
elseif case==22

    Bm = 50
    omega = Inf

    N = 2400
    tmax = 240

    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)

    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    tau = gammap/etap/(a*b*c)^(1/3)
    tmax *= tau
    omega /= tau

    par = elparameters(0.1)
    zc = nothing

end



using JLD
dire = homedir()*"/SimulationData" #Storage.SimulationData

if calculate==true

    if Bm==0
        info("relaxation")
        velocityrelax(t,points,faces) = velocity(t,points,faces;Htime=nothing)
        memory = integrate(points,faces,velocityrelax;N=N,tmax=tmax,zc=zc,par=par)
    else
        if omega==0
            info("stretching")
            if typeof(H0) <: Real
                velocityconst(t,points,faces) = velocity(t,points,faces;Htime=H0*[1.,0,0])
            else
                velocityconst(t,points,faces) = velocity(t,points,faces;Htime=H0)
            end
            memory = integrate(points,faces,velocityconst;N=N,tmax=tmax,zc=zc,par=par)
        elseif omega==Inf
            info("infinitely fast rotating field")
            memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=zc,par=par)
        else
            info("dynamics of rotating field")            
            memory = integrate(points,faces,velocityrot;N=N,tmax=tmax,zc=zc,par=par)
        end
    end

    # Saving calculation
    if !isdir(dire*"/"*session)
        mkdir(dire*"/"*session)
    end
    save("$dire/$session/$case.jld","memory",memory)
else
    data = load(dire*"/$session/$case.jld")
    memory = data["memory"]
end

    
    
    
# 
