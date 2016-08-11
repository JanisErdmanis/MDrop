Bm = 0

N = 480
tmax = 9

a,b,c = 1.1,1,1
points,faces = EllipsoidMeshLoad(a,b,c,0.2)
    
tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau

curvatureless = true
par = nothing
zc = nothing

