#!/usr/bin/julia
# This file represents a wrapper for a dimensionless units

if length(ARGS)>0
    calcpar = ARGS[1] #parse(ASCIIString,ARGS[1])
elseif !isdefined(:calcpar)
    calcpar = nothing
end

fname = basename(calcpar)[1:end-3]

ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"
using JLD
#using Storage
using SurfaceGeometry

include("Calculation/meshes.jl")
include("Calculation/drivers.jl")
include("Calculation/integrator.jl")

# In some far future
# session = dirname(calcpar)
session = "calcpar"

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

# When it will be tested
#!(calcpar==nothing) || include(calcpar)
include(calcpar)

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
    save("$dire/$session/$fname.jld","memory",memory)
else
    data = load(dire*"/$session/$fname.jld")
    memory = data["memory"]
end

