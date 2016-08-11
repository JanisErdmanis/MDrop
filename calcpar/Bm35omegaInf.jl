Bm = 35
omega = Inf

gammap = 1/0.35
etap = 1

N = 900
tmax = 90/0.35

a,b,c = 1,1,1
(points,faces)=EllipsoidMeshLoad(a,b,c,0.2)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = etap*(a*b*c)^(1/3)/gammap
tmax *= tau
omega /= tau

par = elparameters(0.2)
zc = nothing

