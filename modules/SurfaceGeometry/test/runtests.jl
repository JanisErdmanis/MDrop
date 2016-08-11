using SurfaceGeometry
using Base.Test

# Pkg.test("SurfaceGeometry")
# write your own tests here
#@test 1 == 2

### Loading of spherical mesh

function SphereError(points)
    error = 0
    for xkey in 1:size(points,2)
        pos = points[:,xkey]
        n = pos/norm(pos)
        err = norm(n - pos)
        if err>error
            error = err
        end
    end
    return error
end

### Loading of mesh for other tests
println("Mesh loading test")
using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

@test SphereError(points) < 0.01

println("Mesh generation with CGAL")

### Interface will change to pass signed distance function
### CGAL mesh generator

mesher = CGALSurfaceMesher()
fdis(x,y,z) = x^2 + y^2 + z^2 - 1
points, faces = SurfaceMesh(fdis,mesher)
@test SphereError(points) < 0.01

### Matlab mesh generator
println("Mesh generation with distmesh")
try
    mesher = DistmeshSurfaceMesher(0.35)
    points, faces = EllipsoidMesh(1,1,1,mesher)
    @test SphereError(points) < 0.01
catch
    warn("does not work")
end

### Testing a pushback
# sdist(x) = x[1]^2 + x[2]^2 + x[3]^2 - 1
# @test isapprox(norm(pushback(sdist,[1.,1.,1.])),1)

### Testing topology

println("Topology function tests")

t = [5,7,10,7,5,6,4,0,3,0,4,6,4,7,6,4,9,10,7,4,10,0,2,1,2,0,6,2,5,1,5,2,6,8,4,3,4,11,9,8,11,4,9,11,3,11,8,3]
faces = reshape(t,(3,div(length(t),3))) + 1 

triangles = []
for i in FaceVRing(5,faces)
    #println("i")
    push!(triangles,i)
end

@test sort(triangles)==[3,4,5,6,7,12,13,14]

triverticies = []
for i in DoubleVertexVRing(5,faces)
    #println("i")
    push!(triverticies,i)
end

@test (1,4) in triverticies

verticies = []
for i in VertexVRing(5,faces)
    push!(verticies,i)
end

@test sort(verticies)==[1,4,7,8,9,10,11,12] #[1,2,6,7]

println("Testing Complex DS")

points = zeros(size(faces)...)
fb = FaceBasedDS(faces)

triangles = []
for i in FaceVRing(5,fb)
    push!(triangles,i)
end
@test sort(triangles)==[3,4,5,6,7,12,13,14]

triverticies = []
for i in DoubleVertexVRing(5,fb)
    push!(triverticies,i)
end

@test (1,4) in triverticies

verticies = []
for i in VertexVRing(5,fb)
     push!(verticies,i)
end

@test sort(verticies)==[1,4,7,8,9,10,11,12]

println("Topology tests for Connectivity table")

using JLD
data = load("sphere.jld")
faces = data["faces"]

con = ConnectivityDS(faces,10)

verticies = []
for i in VertexVRing(1,con)
    push!(verticies,i)
end

@test sort(verticies)==[2,4,5,14,15]

triverticies = []
for i in DoubleVertexVRing(1,con)
    push!(triverticies,i)
end

@test (2,14) in triverticies

# Here I will also have tests for connectivity table

##### Surface Properties ######

println("Tests for surface properties")

using JLD
data = load("sphere.jld")
points = data["points"]
faces = data["faces"]

n = Array(Float64,size(points)...)
NormalVectors!(n,points,faces,i->FaceVRing(i,faces))
curvatures = Array(Float64,size(points,2))
MeanCurvatures!(curvatures,points,faces,n,i->FaceVRing(i,faces))

angle = 0
curvaturer = 0

for xkey in 1:size(points,2)
    iter = FaceVRing(xkey,faces)
    ncalc = n[:,xkey]
    ccalc = curvatures[xkey]

    nT = points[:,xkey]/norm(points[:,xkey])
    angl = acos(dot(ncalc,nT))
    if angl>angle
        angle = angl
    end

    cT = 1.
    N = size(points,2)
    curvaturer += abs((ccalc - cT)/cT)/N
end

@test angle*180/pi < 1.3
@test curvaturer < 0.01

@test isapprox(volume(points,faces),4/3*pi,atol=0.05)
@test isapprox(*(FitEllipsoid(points)...),1)

### And now integrators

println("With Eigene test checking surface integrators")

N = 100
PasiveStabilisation = false
include("integrator.jl")
@test SphereError(p2) < 0.0008

### ElTopo mesh stabilisation

mesher = CGALSurfaceMesher()
fdis(x,y,z) = x^2 + y^2 + z^2 - 1
p,t = SurfaceMesh(fdis,mesher)

par = Elparameters()
inmsh_verticies = 0.15*p .+ [0.35,0.35,0.35]
inmsh_triangles = map(Int32,t) - 1

points, faces = improvemesh(inmsh_verticies,inmsh_triangles,par)

