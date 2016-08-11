Bm = 0

N = 480
tmax = 9/2

a,b,c = 2, 1/4, 1/4
points,faces = EllipsoidMeshLoad(a,b,c,0.1)
    
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau

adaptive = true
curvatureless = true 
par = elparameters(0.1)
zc = SurfaceGeometry.Erdmanis2016(C=0.01,ftol=1e-3)

