ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using Storage
using SurfaceGeometry

include("constantfield.jl")

# addprocs(4)
# @everywhere module Calc
# include("constantfield.jl")
# export constantfieldcalc
# end
# @everywhere using Storage
# @everywhere using Calc

session = "stabilisationtest"

gammap = 1
etap = 1
mup = 10
N = 50
tmax = 10
PasiveStabilisation = true

H0 = [3,0,0]
a,b,c = 1,1,1
(points,faces)=EllipsoidMeshLoad(1,1,1,0.2)

isdefined(:case) || (case=nothing)

if case==0
    
    fname = "0.jld"

    PasiveStabilisation = false
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-3,C=0.4,gamma=0.25)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)

elseif case==1
    
    fname = "1.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.4,gamma=0.25)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    
elseif case==2
    
    fname = "2.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.5,gamma=1.1)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    
elseif case==3
    
    fname = "3.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.,gamma=1.)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)

elseif case==4
    
    fname = "4.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.,gamma=2.)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    
elseif case==5
    
    fname = "5.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.,gamma=4.)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    
elseif case==6
    
    fname = "6.jld"
    
    memory = constantfieldcalc(points,faces,N=N,tmax=tmax,mup=mup,gammap=gammap,etap=etap,H0=H0,PasiveStabilisation=PasiveStabilisation,ftol=1e-6,C=0.,gamma=0.)
    StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    
end


### For viewing results
using JLD
if !isdefined(:memory)
    dire = Storage.SimulationData
    data = load(dire*"/$session/$case.jld")
    memory = data["memory"]
    # dire = Pkg.dir("Storage","simulations",session)
    # @load dire*"/"*"1.jld"
end

if isdefined(:viewmesh) && viewmesh==true 
    isdefined(:ViewerActive) || (ViewerActive=false)
    using Escher
    if !ViewerActive
        include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
        @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
        ViewerActive = true
    end
end
