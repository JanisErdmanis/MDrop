using SurfaceGeometry
#include("modules/meshes.jl")

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--mu"
    help = "Magnetic permeability"
    arg_type = Float64
    default = 10.
    "--Bm"
    help = "Magnetic bond number"
    arg_type = AbstractString #Float64
    default = "10."
    "--omega"
    help = "Rotation frequency in dimensionless units"
    arg_type = Float64
    default = 0.
    "--Dt"
    help = "Stepsize for advancing smulation"
    arg_type = Float64
    default = 0.1
    "--mesh"
    help = "Initial mesh picked from mesh foleder."
    default = "sphere0.2"
    "--config" # The filename could be derived from this one if given!
    help = "Configuration file for changing other less important parameters." # Stored in calcpar directory
    default = nothing
    "--continue"
    help = "Continue by using data from previous simulation."
    #default = false
    action = :store_true
    #arg_type = Bool
    "--simulation"
    help = "The simulation file which is going to be used."
    default = "FastField"
end

parsed_args = parse_args(ARGS, s)

con = parsed_args["continue"]
sim = parsed_args["simulation"]

if !isdir("$datadir/$sim")
    mkdir("$datadir/$sim")
end

mup = parsed_args["mu"]
Bm = eval(parse(parsed_args["Bm"]))
omega = parsed_args["omega"]
h = parsed_args["Dt"]
mesh = parsed_args["mesh"]
tau = 1

config = parsed_args["config"]

using JLD
data = load("meshes/$mesh.jld")
points,faces = data["points"],data["faces"]

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

### Here now I will write a edge length estimator

vareas = zeros(Float64,size(points,2))
for i in 1:size(faces,2)
    v1,v2,v3 = faces[:,i]
    area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
    vareas[v1] += area/3
    vareas[v2] += area/3
    vareas[v3] += area/3
end

Area = sum(vareas)
N = size(faces,2)
scale = sqrt(Area/N*4/sqrt(3))

par = elparameters(scale)
zc = nothing

if config==nothing
    outdir = "$datadir/$sim/$mesh:mu=$mup;Bm=$Bm;omega=$omega"
else
    include("calcpar/$config.jl")
    Bm = H0^2*(volume(points,faces)*3/4/pi)^(1/3)/gammap
    outdir = "$datadir/$sim/$mesh:mu=$mup;Bm=$Bm;omega=$omega:$config"
end

println("######################################")
println("mu=$mup")
println("Bm=$Bm")
println("omega=$omega")
println("gamma=$gammap")
println("eta=$etap")
println("H0=$H0")
println("mesh=$mesh")
println("Dt=$h")
println("scale=$scale")
println("N=$(size(faces,2))")
println("######################################")

if !isinteractive()
    if isfile("$sim.jl")
        if !isdir("$datadir/$sim") 
            mkdir("$datadir/$sim")
        end
        include("$sim.jl")
    else
        error("No simulation $sim found.")
    end
else
    outdir = "$datadir/$sim/test"
end



