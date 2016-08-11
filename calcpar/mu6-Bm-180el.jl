mup = 6
Bm = 180
omega = Inf

N = 1000
tmax = 10

a,b,c = 2, 1/4, 1/4
points,faces = EllipsoidMeshLoad(a,b,c,0.1)
points,faces = subdivision(points,faces)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau
omega /= tau

par = elparameters(0.05)
zc = nothing
