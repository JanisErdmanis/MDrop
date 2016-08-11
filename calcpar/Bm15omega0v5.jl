Bm = 15
omega = 0

gammap = 1
etap = 3

N = 400
tmax = 10

a,b,c = 2,1/4,1/4
points,faces = EllipsoidMeshLoad(a,b,c,0.1)

H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
tau = etap*(a*b*c)^(1/3)/gammap
#tau = gammap/etap/(a*b*c)^(1/3)
tmax *= tau
        
par = elparameters(0.1)
zc = nothing

