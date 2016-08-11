Bm = 22
omega = Inf

N = 600
tmax = 100

using JLD
@load "meshes/criticalmesh.jld"
points[1,:] = 1.1*points[1,:]
points[2,:] = 0.9*points[2,:]

vol = volume(points,faces)
R0 = (3*vol/4/pi)^(1/3)

H0 = sqrt(Bm*gammap/R0)
tau = etap*R0/gammap
tmax *= tau
omega /= tau

par = elparameters(0.2)
zc = nothing


