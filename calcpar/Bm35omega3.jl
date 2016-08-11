Bm = 35
omega = 3

N = 600*2
tmax = 20

a,b,c = 1,1,1
points,faces = EllipsoidMeshLoad(a,b,c,0.2)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau
omega /= tau
        
par = elparameters(0.2)
zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)
