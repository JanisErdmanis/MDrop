Bm = 50
omega = Inf

N = 2400
tmax = 240

a,b,c = 1,1,1
(points,faces)=EllipsoidMeshLoad(a,b,c,0.2)
points,faces = subdivision(points,faces)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau
omega /= tau

par = elparameters(0.1)
zc = nothing

