ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using Storage
using SurfaceGeometry

include("Calculation/drivers.jl")
include("Calculation/integrator.jl")

isdefined(:calculate) || (calculate=true)
session = "constantfield4"

gammap = 1
etap = 1
mup = 10
N = 500
tmax = 20

isdefined(:case) || (case=nothing)

fname = "$case.jld"

parr(scale) = Elparameters(
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
m_dt = 1)

par = parr(0.1)


par = nothing
zc = nothing

if case==1
    
    H0 = [1,0,0]
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(1,1,1,0.2)

elseif case==2
    
    H0 = [2,0,0]
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(1,1,1,0.2)

elseif case==3
    
    H0 = [2.5,0,0]
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(1,1,1,0.2)

elseif case==4

    H0 = [2,0,0]
    a,b,c = 2,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

elseif case==5
    
    H0 = [2.5,0,0]
    a,b,c = 2,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

elseif case==6    

    H0 = [3,0,0]
    a,b,c = 2,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

elseif case==7
    
    H0 = [3.5,0,0]
    a,b,c = 2,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

elseif case==9    

    N = 1600
    tmax=30
    H0 = [6,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

elseif case==10
    
    H0 = [6.25,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

elseif case==11
    
    H0 = [6.25,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

elseif case==12    

    H0 = [6.25,0,0]
    a,b,c = 2,1/6,1/6
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.05)

elseif case==13    

    H0 = [7,0,0]
    a,b,c = 2,1/6,1/6
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.05)

elseif case==14    

    H0 = [8,0,0]
    a,b,c = 2,1/6,1/6
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.05)

elseif case==15

    tmax = 5
    N = 100
    H0 = [6.25,0,0]

    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    points, faces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)

### A new tests
elseif case==16
    
    H0 = [6.25,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = (2/4/4)^(1/3) * 20
        
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N

elseif case==17
    
    H0 = [5,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
        
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N
    
elseif case==18
    
    H0 = [5.7,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   

elseif case==19
    
    H0 = [2,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   

elseif case==20
    
    H0 = [2.8,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   

elseif case==21
    
    H0 = [7,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   

elseif case==22
    
    H0 = [8,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   


elseif case==23
    
    H0 = [9,0,0]
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    N = 400
    tmax = 20
    
    #zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
    par = parr(0.1)
    par.m_dt = tmax/N   
    
end

using JLD
if calculate==true
    velocityconst(t,points,faces) = velocity(t,points,faces;Htime=H0)
    memory = integrate(points,faces,velocityconst;N=N,tmax=tmax,zc=zc,par=par)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
else
    # dire = Pkg.dir("Storage","simulations")
    # data = load(dire*"/"*"/$session/$case.jld")
    dire = Storage.SimulationData
    data = load(dire*"/$session/$case.jld")
    memory = data["memory"]
end

# using JLD
# if !isdefined(:viewmeshactive) || viewmeshactive==false
#     eval(:(using Escher))
#     include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
#     @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
#     viewmeshactive = true
# end

