ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using Storage
using SurfaceGeometry

include("Calculation/drivers.jl")
include("Calculation/integrator.jl")

function rotatingfieldcalc(points,faces;N=10,tmax=3,zc=nothing, par=nothing)
    integrate(points,faces,velocityrot;N=N,tmax=tmax,zc=zc,par=par)
end

session = "rotatingfield"

gammap = 1
etap = 1
mup = 10
N = 1000
tmax = 40
h = tmax/N

isdefined(:calculate) || (calculate=true)

PasiveStabilisation = true

isdefined(:case) || (case=0)

fname = "$case.jld"

zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

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
m_dt = h)

par = parr(0.1)

if case==0
    
    H0 = 6
    omega = 0.1
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=10,tmax=3/10,par=nothing)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==1
    
    H0 = 6
    omega = 0.1
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end
elseif case==2
    
    H0 = 6
    omega = 0.3
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end
elseif case==3

    H0 = 6
    omega = 0.2
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==4
    
    H0 = 6
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end
elseif case==5

    H0 = 6
    omega = 1
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### case 4 but with stabilisation    
elseif case==6
    
    H0 = 6
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

    
elseif case==7
    
    H0 = 8
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end
    

elseif case==8
    
    H0 = 8
    omega = 1
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==9
    
    H0 = 8
    omega = 1
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==10
    
    H0 = 8
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==11
    
    H0 = 10
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

# failed
elseif case==12
    
    H0 = 14
    omega = 0.5
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==13
    
    H0 = 10
    gammap = 1/0.35
    omega = 0.5
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==14
    
    H0 = 10
    gammap = 1/0.35
    omega = 1
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==15
    
    H0 = 10
    gammap = 1/0.35
    omega = 4
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==16
    
    H0 = 6
    gammap = 1/0.35
    omega = 0.5
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==17
    
    H0 = 6
    gammap = 1/0.35
    omega = 1
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==18
    
    H0 = 6
    gammap = 1/0.35
    omega = 4
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==20
    
    H0 = 6
    gammap = 1/0.35
    omega = 4
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==21
    
    H0 = 10
    gammap = 1/0.35
    omega = 4
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 15
    N = 200

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==22
    
    H0 = 6
    gammap = 1/0.35
    omega = 4
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 15
    N = 200

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==23
    
    H0 = 6
    gammap = 1/0.35
    omega = 1
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 15
    N = 200

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### A new set of calculations

elseif case==24
    
    H0 = 10
    gammap = 1/0.35
    omega = 8
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 15
    N = 300

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==25

    H0 = 14
    gammap = 1/0.35
    omega = 8
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 15
    N = 300

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==26
    
    H0 = 6
    gammap = 1/0.35
    omega = 15
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30
    N = 600

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==27
    
    H0 = 14
    gammap = 1/0.35
    omega = 15
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 4.5*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==28
    
    H0 = 14
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 4*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==29
    
    H0 = 10
    gammap = 1/0.35
    omega = 15
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 7*4
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==30
    
    H0 = 10
    gammap = 1/0.35
    omega = 15
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 5*4
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==31
    
    H0 = 8
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==32
    
    H0 = 7.5
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==33
    
    H0 = 8.5
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==34
    
    H0 = 9.0
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==35
    
    H0 = 7.0
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### The next session

elseif case==36
    
    H0 = 10
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==37
    
    H0 = 11
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==38
    
    H0 = 12
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==39
    
    H0 = 13
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### Modified
elseif case==40
    
    H0 = 14
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*5
    N = 300*5

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==41
    
    H0 = 16
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*5
    N = 300*5

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==42
    
    H0 = 14
    gammap = 1/0.35
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*5/2
    N = 300*5

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==43
    
    H0 = 16
    gammap = 1/0.35
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*5
    N = 300*5*2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### Testing
elseif case==44
    
    H0 = 16
    gammap = 1/0.35
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 0.5*9*2
    N = 900*2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==45
    
    H0 = 11
    gammap = 1/0.35
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 5*2
    N = 600*2*2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==46
    
    H0 = 11
    gammap = 1/0.35
    omega = 60*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 5
    N = 600*2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### Real parameters
elseif case==47
    
    H0 = 1.25
    mup = 25
    normalfield = :traditional
    gammap = 1
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 0.5*9
    N = 900

    par = parr(0.2)
    par.m_dt = tmax/N


    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==48
    
    H0 = 1.8
    mup = 25
    normalfield = :traditional
    gammap = 1
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 0.5*9
    N = 900

    par = parr(0.2)
    par.m_dt = tmax/N


    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==49
    
    H0 = 3
    mup = 25
    normalfield = :traditional
    gammap = 1
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 0.5*9
    N = 900

    par = parr(0.2)
    par.m_dt = tmax/N


    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==50
    
    H0 = 10
    omega = 30*2
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 2*9
    N = 900

par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==51
    
    H0 = 10
    omega = 30*2
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 0.5*9/2
    N = 900

par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==52
    
    H0 = 12
    omega = 30*2
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 0.5*9/2
    N = 900

par = parr(0.1)
    par.m_dt = tmax/N


    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==53
    
    H0 = 12
    omega = 30*2
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    tmax = 2*9
    N = 900

par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=zc,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==54
    
    H0 = 16
    gammap = 1/0.35
    omega = 30*2
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 0.5*9
    N = 900

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==55
    
    H0 = 10
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==56
    
    H0 = 14
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 20
    N = 200

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==57
    
    H0 = 7
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==58
    
    H0 = 8
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==59
    
    H0 = 9
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==60
    
    H0 = 8.5
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==61   # Updated
    
    H0 = 12
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3 * 2
    N = 300*3  * 2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==62
    
    H0 = 13
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


elseif case==64
    
    H0 = 14
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==65
    
    H0 = 6
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==66
    
    H0 = 5
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==67

    H0 = 12
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 5
    N = 400

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==68 # Improved
    
    H0 = 16
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
tmax = 5*2 * 2
    N = 400 * 2

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==69
    
    H0 = 16
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
tmax = 30
    N = 400

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==70
    
    H0 = 14
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
tmax = 5 
    N = 400

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==71
    
    H0 = 14
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
tmax = 10 * 2
    N = 400 * 2

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==72

H0 = 12
    gammap = 1/0.35
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 10 * 6
    N = 400 * 6

    par = parr(0.1)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==73 ### Improved
    
    H0 = 12
    gammap = 1/0.35
    omega = 30
    a,b,c = 1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 3*3*2
    N = 300*3*2

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==74 
    
    H0 = 12
    gammap = 1.42
    omega = 30
    
a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

par = parr(0.1)

    tmax = 3*3*2
    N = 300*3*2

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==75
    
    H0 = 12
    gammap = 1.42
    omega = 15
    
a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

par = parr(0.1)

    tmax = 3*3*2
    N = 300*3

    par.m_dt = tmax/N

    if calculate==true
        memory = rotatingfieldcalc(points,faces,N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

### Raimonds
elseif case==76
    
    H0 = 10
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    points,faces = subdivision(points,faces)
    tmax = 60
    N = 400

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==77
    
    H0 = 9
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==78
    
    H0 = 8
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==79
    
    H0 = 7
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==80
    
    H0 = 7.5
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==81
    
    H0 = 6
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==82
    
    H0 = 7
mup = 25
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==83
    
    H0 = 5
mup = 25
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


### A new set today 
elseif case==84
    
    H0 = 7.6
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==85
    
    H0 = 7.7
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==86
    
    H0 = 7.85
    gammap = 1/0.35
    a,b,c = 1.1,1,1
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
    tmax = 30*3
    N = 300*3

    par = parr(0.2)
    par.m_dt = tmax/N

    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==87 ### Longer time
    
    H0 = 12
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    points,faces = subdivision(points,faces)
    tmax = 5*4
    N = 300*4

par = parr(0.05)
    par.m_dt = tmax/N


    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end

elseif case==88
    
    H0 = 14
    a,b,c = 2,1/4,1/4
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)
    points,faces = subdivision(points,faces)
    tmax = 5
    N = 300

par = parr(0.05)
    par.m_dt = tmax/N


    if calculate==true
        memory = integrate(points,faces,velocityrotfast;N=N,tmax=tmax,zc=nothing,par=par)
        StoreData(fname,session=session,memory=memory,a=a,b=b,c=c,etap=etap,gammap=gammap,H0=H0,mup=mup)
    end


end

using JLD
if calculate==true
    ### An option is to put calculation here
else
    #dire = Pkg.dir("Storage","simulations")
    dire = Storage.SimulationData
    data = load(dire*"/$session/$case.jld")
    memory = data["memory"]
end

# using JLD
# isdefined(:viewmeshactive) || (viewmeshactive=calculate) # For avoiding unnecesar
# if viewmeshactive==false
#     eval(:(using Escher))
#     include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
#     @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
#     viewmeshactive = true
# end




# using JLD
# if calculate==false
#     @load "/home/janiserdmanis/.julia/v0.4/Storage/simulations/$session/$case.jld"

#     if !isdefined(:viewmeshactive) || viewmeshactive==false
#         eval(:(using Escher))
#         include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
#         @spawn escher_serve(5555,Pkg.dir("SurfaceGeometry","examples","viewers"))
#         viewmeshactive = true
#     end
    
# end

