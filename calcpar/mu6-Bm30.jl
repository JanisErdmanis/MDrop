mup = 6
Bm = 30
omega = Inf

N = 500
tmax = 50

a,b,c = 1,1,1
(points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau
omega /= tau

par = elparameters(0.2)
zc = nothing
