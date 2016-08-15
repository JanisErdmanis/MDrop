using SurfaceGeometry
include("modules/meshes.jl")

datadir = homedir()*"/SimulationData/"
outdir = "test"
### Again I would need to have a wa of knowing if code is being included
points,faces = EllipsoidMeshLoad(1,1,1,0.2)

mup = 10
Bm = 10
omega = 10
tau = 1 ### Unit of time
h = 0.1

function dimensionless()
    #tau = gammap/etap/(volume(points,faces)*3/4/pi)^(1/3)

    global gammap = 1
    global etap = gammap/tau/(volume(points,faces)*3/4/pi)^(1/3)
    global H0 = sqrt(Bm*gammap/(volume(points,faces)*3/4/pi)^(1/3))
    #global h *= tau
    #global omega /= tau
end

dimensionless()

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

par = elparameters(0.2)
zc = nothing
