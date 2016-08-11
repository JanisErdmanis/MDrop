mup = 6

Bm = 100
omega = Inf

N = 600
tmax = 100

#a,b,c = 1,1,1
using JLD
(points,faces)=EllipsoidMeshLoad(1.1,1,1,0.2)

vol = volume(points,faces)
R0 = (3*vol/4/pi)^(1/3)

# (points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
# points,faces = subdivision(points,faces)

H0 = sqrt(Bm*gammap/R0)
tau = gammap/etap/R0
tmax *= tau
omega /= tau

par = elparameters(0.1)
zc = nothing
