using SurfaceGeometry
using ElTopo
using JLD

mup = 10.
Bm = 10.
omega = 0.
h = 0.1
mesh = "sphere0.2"
tau = 1
con = false ### do not continue simulation

scale = 0.2

if !isdir("$datadir/config")
    mkdir("$datadir/config")
end
outdir = "$datadir/config/$mesh:mu=$mup;Bm=$Bm;omega=$omega"

### loading mesh
data = load("meshes/$mesh.jld")
points,faces = data["points"],data["faces"]

### Dimensional parameters
gammap = 1.
etap = gammap/tau/(volume(points,faces)*3/4/pi)^(1/3)
H0 = sqrt(Bm*gammap/(volume(points,faces)*3/4/pi)^(1/3))

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
 m_allow_vertex_movement = false, #true, ### This is where is a bug
 m_use_curvature_when_collapsing = false,
 m_use_curvature_when_splitting = false,
 m_dt = 1
)

par = elparameters(scale)
zc = nothing

