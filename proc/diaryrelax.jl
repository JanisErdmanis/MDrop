ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using Storage
using SurfaceGeometry

include("Calculation/drivers.jl")
include("Calculation/integrator.jl")

isdefined(:calculate) || (calculate=true)
session = "relaxation"

gammap = 1
etap = 1
mup = 1
N = 100
tmax = 9

#scale = 0.1
parf(scale) = Elparameters(
                   m_use_fraction = false,
                   m_min_edge_length = 0.7*scale,
                   m_max_edge_length = 1.5*scale,
                   m_max_volume_change = 0.1*scale^3,
                   m_min_curvature_multiplier = 1.0,
                   m_max_curvature_multiplier = 1.0,
                   m_merge_proximity_epsilon = 0.5*scale,
                   m_proximity_epsilon = 0.00001,
m_perform_improvement = true, 
m_collision_safety = false,
m_min_triangle_angle = 15,
m_max_triangle_angle = 120,
m_allow_vertex_movement = true, ### 
m_use_curvature_when_collapsing = false,
m_use_curvature_when_splitting = false,
m_dt = h)

zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

fname = "$case.jld"
isdefined(:case) || (case=1)

if case==1
    
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    zc = nothing
    par = nothing

elseif case==2

    N = 100
    tmax = 9
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    par = parf(0.1)

elseif case==3

    N=480
    a,b,c = 1.1,1,1
    zc = nothing
    par = nothing
    tmax = 9
    curvatureless = false
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

elseif case==4

    N=480
    a,b,c = 1.1,1,1
    zc = nothing
    par = nothing
    tmax = 9
    curvatureless = true
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

elseif case==5
    
    N = 100
    tmax = 9
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    par = parf(0.1)
    
end

using JLD
if calculate==true
    velocityrelax(t,points,faces) = velocity(t,points,faces;Htime=nothing)
    memory = integrate(points,faces,velocityrelax;N=N,tmax=tmax,zc=zc,par=par)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
else
    # dire = Pkg.dir("Storage","simulations")
    # data = load(dire*"/"*"/$session/$case.jld")
    dire = Storage.SimulationData
    data = load(dire*"/$session/$case.jld")
    memory = data["memory"]
end

using JLD
if !isdefined(:viewmeshactive) || viewmeshactive==false
    eval(:(using Escher))
    include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
    @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
    viewmeshactive = true
end


