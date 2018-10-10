using SurfaceGeometry
using LinearAlgebra

eye(A::AbstractMatrix{T}) where T = Matrix{eltype(A)}(I,size(A))
eye(m::Integer) = Matrix(1.0I,m,m)
diagm(x) = Matrix(Diagonal(x))


### The Best way of making calculation now
### One which needs signifficant improvement for speed
function surfacefield(points,faces,normals,mup,Htime)

    psi = PotentialSimple(points,faces,mup,Htime)
    H = HField(points,faces,psi)
    Ht = Array{Float64}(undef,size(points,2))
    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        P = eye(3) - nx*nx'
        Ht[xkey] = norm(P*H[:,xkey])
    end
    
    # rpoints, rfaces = subdivision(points,faces; method=:paraboloid)
    # Hn = NormalFieldTrapezodial(rpoints,faces,rfaces,mup,H0;NP=3)
    #if normalfield==:current
    Hn = NormalFieldCurrent(points,faces,H,mup,Htime)
    # elseif normalfield==:traditional
    #     Hn = NormalField(points,faces,mup,Htime;regularize=false)
    # end

    return psi,Ht,Hn
end    

### Old way
# using Cubature
# function strquad(q::Function,x1,x2,x3;abstol=0)
    
#     x(xi,mu) = x1*(1 - xi - mu) + x2*xi + x3*mu

#     xi(rho,hi) = rho*cos(hi)
#     mu(rho,hi) = rho*sin(hi)
#     R(hi) = 1/(cos(hi) + sin(hi))
    
#     B = dot(x3-x1,x2-x1)/norm(x2-x1)^2
#     C = norm(x3-x1)^2/norm(x2-x1)^2

#     ff(t,hi) = R(hi)*q(xi(t*R(hi),hi),mu(t*R(hi),hi)) / sqrt(cos(hi)^2 + B*sin(2*hi) + C*sin(hi)^2)

#     hS = norm(cross(x2-x1,x3-x1))
#     hcubature(x -> ff(x[1],x[2])*hS/norm(x2-x1),[0,0],[1,pi/2],abstol=abstol)    
# end

### For a better performance and simplicity
using FastGaussQuadrature
const NP = 10
const t, w = gausslegendre(NP)
function strquad(q::Function,x1,x2,x3)

    B = dot(x3-x1,x2-x1)/norm(x2-x1)^2
    C = norm(x3-x1)^2/norm(x2-x1)^2
    hS = norm(cross(x2-x1,x3-x1))

    s = 0 
    for i in 1:NP
        Chi = pi/4*(1 + t[i])

        R = 1/(cos(Chi) + sin(Chi))
        si = 0
        for j in 1:NP
            rho = R/2*(1 + t[j])
            si += q(rho*cos(Chi),rho*sin(Chi))*w[j]
        end

        s += si*R/2 / sqrt(cos(Chi)^2 + B*sin(2*Chi) + C*sin(Chi)^2) * w[i]
    end
    s *= pi/4

    return s*hS/norm(x2 - x1) 
end

function NormalFieldCurrent(points,faces,Ht,hmag,H0; eps=0.0001, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        Ht[:,xkey] = (eye(3) - nx*nx')*Ht[:,xkey]
    end
    
    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end


    function qs(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        s = - dot(nx,cross(Hty - Htx,cross(ny,-(y-x)/norm(y-x)^2)))
    end

    Hn = Array{Float64}(undef,size(points,2))
    
    for xkey in 1:size(points,2)

        nx = normals[:,xkey]
        x = points[:,xkey] + eps*nx
        Htx = Ht[:,xkey]

        s = 0 
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            ny = normals[:,ykey]
            Hty = Ht[:,ykey]

            #s += dot(nx,-(Hty-Htx)*dot((y-x)/norm(y-x)^3,ny)) * vareas[ykey]
            s += dot(nx,-(Hty-Htx)*dot((y-x)/norm(y-x)^3,ny)) * vareas[ykey]
            s += -dot(nx,cross(Hty-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * vareas[ykey]
            #s += dot(nx,cross(cross(ny,),-(y-x)/norm(y-x)^3)) * vareas[ykey]
        end

        ### Making a proper hole
        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2

            ny = normals[:,v2]
            y = points[:,v2]
            s -= -dot(nx,cross(Ht[:,v2]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

            ny = normals[:,v3]
            y = points[:,v3]
            s -= -dot(nx,cross(Ht[:,v3]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

            ### Singular triangle integration

            #s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,points[:,v2],points[:,v3],abstol=abs(s/100))[1]
            s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,points[:,v2],points[:,v3])
        end

        #println("xkey is $xkey")

        Hn[xkey] = dot(H0,nx)/hmag + 1/4/pi * (1-hmag)/hmag * s
    end

    return Hn
end

function NormalFieldCurrentGaussian(points,faces,Ht,hmag,H0; eps=0.0001, NP=3, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        Ht[:,xkey] = (eye(3) - nx*nx')*Ht[:,xkey]
    end

    xw = gaussianpoints(NP)

    function fs(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        #area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        s = dot(nx,-(Hty - Htx)*dot((y-x)/norm(y-x)^3,ny)) - dot(nx,cross(Hty - Htx,cross(ny,-(y-x)/norm(y-x)^3)))

    end

    # function qs(xi,eta,v1,v2,v3,x,nx,Htx)
    #     y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
    #     Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
    #     ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
    #     #area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
    #     s = dot(nx,-(Hty - Htx)*dot((y-x)/norm(y-x)^2,ny)) - dot(nx,cross(Hty - Htx,cross(ny,-(y-x)/norm(y-x)^2)))

    # end

    function fss(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        #area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        s = dot(nx,-(Hty - Htx)*dot((y-x)/norm(y-x)^3,ny)) 

    end

    function fs2(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        # ny = cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])
        # ny = ny/norm(ny)
        #area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        s = - dot(nx,cross(Hty,cross(ny,-(y-x)/norm(y-x)^3)))
    end


    function qs(xi,eta,v1,v2,v3,x,nx,Htx)
        y = (1 - xi - eta)*points[:,v1] + xi*points[:,v2] + eta*points[:,v3]
        Hty = (1 - xi - eta)*Ht[:,v1] + xi*Ht[:,v2] + eta*Ht[:,v3]
        ny =  (1 - xi - eta)*normals[:,v1] + xi*normals[:,v2] + eta*normals[:,v3]
        #area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        s = - dot(nx,cross(Hty - Htx,cross(ny,-(y-x)/norm(y-x)^2)))

    end

    Hn = Array(Float64,size(points,2))
    
    for xkey in 1:size(points,2)

        nx = normals[:,xkey]
        x = points[:,xkey] + eps*nx
        Htx = (eye(3) - nx*nx')*Ht[:,xkey]

        s = 0

        # f1(vi) = !(xkey==vi) ? dot(nx,-(Ht[:,vi]-Htx)*dot((points[:,vi]-x)/norm(points[:,vi]-x)^3,normals[:,vi])) : 0
        # f2(vi) = !(xkey==vi) ? -dot(nx,cross(Ht[:,vi]-Htx,cross(normals[:,vi],-(points[:,vi]-x)/norm(points[:,vi]-x)^3))) : 0 
        
        for i in 1:size(faces,2)
            !(xkey in faces[:,i]) || continue
            v1,v2,v3 = faces[:,i]
            area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2

            #function fs(xi,eta,ti,x,nx,Htx)
            
            s += (fs(0,0,v1,v2,v3,x,nx,Htx) + fs(1,0,v1,v2,v3,x,nx,Htx) + fs(0,1,v1,v2,v3,x,nx,Htx)) * area/3
            
        end

        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            v1 = xkey
            area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
            
            x2 = points[:,v2]
            x3 = points[:,v3]

            ### Gaussian quadrature

            # for j in 1:NP
            #     xi,eta,w = xw[j,1],xw[j,2],xw[j,3]
                
            #     s += w/2*fs2(xi,eta,v1,v2,v3,x,nx,Htx)*area
            #     # f1 += w/2*f*(1-xi-eta)
            #     # f2 += w/2*f*xi
            #     # f3 += w/2*f*eta
            # end

            ### Adaptive Gaussian quadrature

            # [ ] Need to fix variable s

            ##
                        
            # ds = Inf
            # ksi0 = 1
            # eta0 = 1

            # i = 1
            # ss = 0
            
            # for i in 1:50  ### Now for simplicity
            #     ksi0 /=2
            #     eta0 /=2

            #     s1 = ss
                
            #     for j in 1:NP
            #         xi,eta = ksi0 + ksi0*xw[j,1],0. + eta0*xw[j,2]
            #         local w = xw[j,3]
            #         f = fs2(xi,eta,v1,v2,v3,x,nx,Htx)
            #         ss += w/2*ksi0*eta0*f
            #     end

            #     for j in 1:NP
            #         xi,eta = 0. + ksi0*xw[j,1],eta0 + eta0*xw[j,2]
            #         local w = xw[j,3]
            #         f = fs2(xi,eta,v1,v2,v3,x,nx,Htx)
            #         ss += w/2*ksi0*eta0*f
            #     end

            #     for j in 1:NP
            #         xi,eta = ksi0 - ksi0*xw[j,1],eta0 - eta0*xw[j,2]
            #         local w = xw[j,3]
            #         f = fs2(xi,eta,v1,v2,v3,x,nx,Htx)
            #         ss += w/2*ksi0*eta0*f
            #     end

            #     s2 = ss

            #     ds_new = norm(s2-s1)
            #     # if ds_new<ds
            #     #     ds = ds_new
            #     # else
            #     #     break
            #     # end
                    
            # end

            #println("i is $i")

            # s += ss*area
            ### Ordinary quadrature            
            s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,x2,x3)

            s += (fss(1,0,v1,v2,v3,x,nx,Htx) + fss(0,1,v1,v2,v3,x,nx,Htx)) * area/3
        end
        
        Hn[xkey] = dot(H0,nx)/hmag + 1/4/pi * (1-hmag)/hmag * s
    end

    return Hn
end


function NormalFieldCurrentOld(points,faces,Ht,hmag,H0; eps=0.0001, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    normals = Array{Float64}(undef,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros{Float64}(undef,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    Hn = Array{Float64}(undef,size(points,2))
    
    for xkey in 1:size(points,2)

        nx = normals[:,xkey]
        x = points[:,xkey] + eps*nx

        s = 0 
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            ny = normals[:,ykey]

            s += dot(nx,cross(cross(ny,Ht[:,ykey]),-(y-x)/norm(y-x)^3)) * vareas[ykey]
        end

        ### Making a proper hole
        # for (v2,v3) in DoubleVertexVRing(xkey,faces)
        #     area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2

        #     ny = normals[:,v2]
        #     y = points[:,v2]
        #     s -= dot(nx,cross(cross(ny,Ht[:,v2]),-(y-x)/norm(y-x)^3))*area/3

        #     ny = normals[:,v3]
        #     y = points[:,v3]
        #     s -= dot(nx,cross(cross(ny,Ht[:,v3]),-(y-x)/norm(y-x)^3))*area/3
        # end

        Hn[xkey] = dot(H0,nx)/hmag + 1/4/pi * (1-hmag)/hmag * s
    end

    return Hn
end

function NormalFieldDomain(points,faces,psi,hmag,H0; eps=0.0001, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    Hn = Array{Float64}(undef,size(points,2))
    
    for xkey in 1:size(points,2)

        nx = normals[:,xkey]
        x = points[:,xkey] + eps*nx
        psix = psi[xkey]

        s = 0 
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            ny = normals[:,ykey]
            psiy = psi[ykey]

            #s += (psiy-psix)*dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
            s += psiy*dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
        end

        psixp = 2*dot(H0,x)/(hmag+1) + 1/2/pi * (hmag-1)/(hmag+1) * s

        Hn[xkey] = (psixp - psix)/eps
    end

    return Hn
end

function NormalFieldHypersingular(points,faces,psi,Ht; normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros{Float64}(undef,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    A = zeros(Float64,size(points,2),size(points,2))

    for xkey in 1:size(points,2)
        x = points[:,xkey]
        nx = points[:,xkey]
        Htx = Ht[:,xkey]
        psix = psi[xkey]

        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            
            y= points[:,ykey]
            ny = points[:,ykey]
            Hty = Ht[:,xkey]
            psiy = psi[xkey]

            A[ykey,xkey] = -1/4/pi*dot(y-x,nx)/norm(y-x)^3 * vareas[ykey]
        end
    end

    ### diognal elements

    for xkey in 1:size(points,2)
        x = points[:,xkey]
        nx = points[:,xkey]
        Htx = Ht[:,xkey]
        psix = psi[xkey]

        s = 0 
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            
            y= points[:,ykey]
            ny = points[:,ykey]
            Hty = Ht[:,xkey]
            psiy = psi[xkey]

            G = 1/4/pi*( dot(nx,ny)/norm(y-x)^3 - 3 * dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^5 )
            s += - G * dot(y-x,nx) * vareas[ykey]
            
            s += 1/4/pi*dot(y-x,nx)/norm(y-x)^3 * dot(nx,ny) * vareas[ykey]
            
        end

        A[xkey,xkey] = 1 + s
    end

    B = zeros(Float64,size(points,2))

    for xkey in 1:size(points,2)
        x = points[:,xkey]
        nx = points[:,xkey]
        Htx = Ht[:,xkey]
        psix = psi[xkey]

        s = 0 
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            
            y= points[:,ykey]
            ny = points[:,ykey]
            Hty = Ht[:,xkey]
            psiy = psi[xkey]

            G = 1/4/pi*( dot(nx,ny)/norm(y-x)^3 - 3 * dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^5 )

            s += -G * (psiy - psix - dot(Htx,y-x) ) * vareas[ykey]
            s += -1/4/pi * dot(y-x,nx)/norm(y-x)^3 * dot(ny,Htx) * vareas[ykey]
        end

        B[xkey] = s
        
    end

    A = A'
    
    return A\B
end
    



function NormalFieldRecalculated(points,faces,psi; normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    
    B = Array{Float64}(undef,size(points,2))

    for xkey in 1:size(points,2)
        x = points[:,xkey]
        psix = psi[xkey]

        s = 0
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            
            y = points[:,ykey]
            ny = normals[:,ykey]
            psiy = psi[ykey]
            
            s += -dot(y-x,ny)/norm(y-x)^3 * (psiy - psix) * vareas[ykey]
        end

        B[xkey] = s
    end

    A = zeros(Float64,size(points,2),size(points,2))
    for xkey in 1:size(points,2)
        x = points[:,xkey]

        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            
            A[ykey,xkey] = 1/norm(y-x) * vareas[ykey]
        end
    end

    ### and now the singular quadratures
    for xkey in 1:size(points,2)
        x1 = points[:,xkey]

        s = 0
        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            x2 = points[:,v2]
            x3 = points[:,v3]

            s += strquad((xi,eta) -> (1 - xi - eta),x1,x2,x3)
            #s += strquad((xi,eta) -> 1,x1,x2,x3)[1]
        end

        A[xkey,xkey] = s
        
    end


    # for (v2,v3) in DoubleVertexVRing(xkey,faces)
    #     area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2

    #     ny = normals[:,v2]
    #     y = points[:,v2]
    #     s -= -dot(nx,cross(Ht[:,v2]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

    #     ny = normals[:,v3]
    #     y = points[:,v3]
    #     s -= -dot(nx,cross(Ht[:,v3]-Htx,cross(ny,-(y-x)/norm(y-x)^3))) * area/3

    #     ### Singular triangle integration

    #     s += strquad((xi,eta) -> qs(xi,eta,xkey,v2,v3,x,nx,Htx),x,points[:,v2],points[:,v3],abstol=abs(s/100))[1]
    # end


    
    A = A'

    return A \ B
end

function HtField(points,faces,psi;normals=nothing)
    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    H = HField(points,faces,psi)
    for xkey in 1:size(points,2)
        nx = normals[:,xkey]
        H[:,xkey] = (eye(3) - nx*nx')*H[:,xkey]
    end

    return H
end

function HField(points,faces,psi)

    # if normals==nothing
    #     normals = Array(Float64,size(points)...)
    #     NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    # end

    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    H = Array{Float64}(undef,size(points)...)

    for xkey in 1:size(points,2)

        x = points[:,xkey]
        #nx = normals[:,xkey]
        psix = psi[xkey]
        #Hnx = Hn[xkey]

        vvec = Float64[]
        dphi = Float64[]

        for ykey in VertexVRing(xkey,faces)
            y = points[:,ykey]

            vvec = [vvec; y-x]
            dphi = [dphi; psi[ykey]-psix]
            # println(ykey)
        end

        # vvec = [vvec; nx]
        # dphi = [dphi; Hnx]    
        
        A, B = vvec, dphi

        A = transpose(reshape(A,3,div(length(A),3))) ### This looks unecesarry
        H[:,xkey] = inv(transpose(A)*A)*transpose(A)*B

    end

    return H
end

function NormalField(points,faces,hmag,H0; regularize=false, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    A = zeros(Float64,size(points,2),size(points,2))
    B = zeros(Float64,size(points,2))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]
        
        for ykey in 1:size(points,2)
            !(xkey==ykey) || continue
            y = points[:,ykey]
            A[ykey,xkey] = dot(y-x,nx)/norm(y-x)^3 * vareas[ykey]
        end
        B[xkey] = dot(H0,nx)
    end


    ### Regularized version
    if regularize==true
        # C = Array(Float64,size(points,2))
        # for xkey in 1:size(points,2)
        #     x = points[:,xkey]
        #     s = 0
        #     for ykey in 1:size(points,2)
        #         !(xkey==ykey) || continue
        #         y = points[:,ykey]
        #         ny = normals[:,ykey]
        #         s += - dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
        #     end
        #     C[xkey] = s
        # end
         
        # B = B*2/(hmag+1)
        # A = A*(hmag-1)/(hmag+1)/2/pi
        # A = transpose(A)

        # Hn = (eye(A)*(1 + (hmag-1)/(hmag+1)) + 1/2/pi * (hmag-1)/(hmag+1)*diagm(C) + A)\B
        A = A'
        B = B*2/(hmag+1)
        A = A - diagm(Float64[sum(A[i,:]) for i in 1:size(A,2)])
        #Hn = (eye(A)*(1 - (hmag-1)/(hmag+1)) + 1/2/pi * (hmag-1)/(hmag+1) * A ) \ B
        Hn = (eye(A)*(1 + (hmag-1)/(hmag+1)/2) + 1/2/pi * (hmag-1)/(hmag+1) * A ) \ B
    else
        A = A'
        B = B*2/(hmag+1)
        A = A*(hmag-1)/(hmag+1)/2/pi
        A = transpose(A)
        Hn = (eye(A) + A)\B
    end

    return Hn
end


function PotentialSimple(points,faces,hmag,H0;regularize=true,normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    A = zeros(Float64,size(points,2),size(points,2))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
        
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]

        for ykey in 1:size(points,2)
            if xkey==ykey
                continue
            end

            y = points[:,ykey]
            ny = normals[:,ykey]
            
            A[ykey,xkey] = dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
        end
    end

    B = zeros(Float64,size(points,2))
    for xkey in 1:size(points,2)
        B[xkey] = 2*dot(H0,points[:,xkey])/(hmag+1)
    end


    if regularize==true
        A = A'
        psi = (eye(A)*(1- (hmag-1)/(hmag+1)) - 1/2/pi * (hmag-1)/(hmag+1) * (A - diagm(Float64[sum(A[i,:]) for i in 1:size(A,2)]))) \ B
    else
        A = A*(hmag-1)/(hmag+1)/2/pi
        A = A'
        psi = (eye(A) - A)\B
    end
    
    return psi
end



#"Skip singular, integrate others with Gaussian. Start with objective to get 0.05 % error going over triangles"


function gaussianpoints(NP)

    if NP==3
        [0.16666666666667    0.16666666666667    0.33333333333333
              0.16666666666667    0.66666666666667    0.33333333333333
              0.66666666666667    0.16666666666667    0.33333333333333]
    elseif NP==12
        [0.24928674517091 0.24928674517091 0.11678627572638
         0.24928674517091 0.50142650965818 0.11678627572638
         0.50142650965818 0.24928674517091 0.11678627572638
         0.06308901449150 0.06308901449150 0.05084490637021
         0.06308901449150 0.87382197101700 0.05084490637021
         0.87382197101700 0.06308901449150 0.05084490637021
         0.31035245103378 0.63650249912140 0.08285107561837
        0.63650249912140 0.05314504984482 0.08285107561837
        0.05314504984482 0.31035245103378 0.08285107561837
        0.63650249912140 0.31035245103378 0.08285107561837
        0.31035245103378 0.05314504984482 0.08285107561837
        0.05314504984482 0.63650249912140 0.08285107561837]
    elseif NP==16
        [0.33333333333333    0.33333333333333    0.14431560767779
            0.45929258829272    0.45929258829272    0.09509163426728
            0.45929258829272    0.08141482341455    0.09509163426728
            0.08141482341455    0.45929258829272    0.09509163426728
            0.17056930775176    0.17056930775176    0.10321737053472
            0.17056930775176    0.65886138449648    0.10321737053472
        0.65886138449648    0.17056930775176    0.10321737053472
        0.05054722831703    0.05054722831703    0.03245849762320
        0.05054722831703    0.89890554336594    0.03245849762320
        0.89890554336594    0.05054722831703    0.03245849762320
        0.26311282963464    0.72849239295540    0.02723031417443
        0.72849239295540    0.00839477740996    0.02723031417443
        0.00839477740996    0.26311282963464    0.02723031417443
        0.72849239295540    0.26311282963464    0.02723031417443
        0.26311282963464    0.00839477740996    0.02723031417443
        0.00839477740996    0.72849239295540    0.02723031417443]
    else
        error("Gaussian quadrature with $NP points not implemented")
    end
end


macro calculationf()
    return quote
        zeta = 1-xi-eta

        y = zeta*(2*zeta-1)*x1 + xi*(2*xi-1)*x2 + eta*(2*eta-1)*x3 + 4*xi*zeta*x4 + 4*xi*eta*x5 + 4*eta*zeta*x6
        exi = (1-4*zeta)*x1 + (4*xi-1)*x2 + 0*x3 + 4*(zeta-xi)*x4 + 4*eta*x5 - 4*eta*x6
        eeta = (1-4*zeta)*x1 + 0*x2 + (4*eta-1)*x3 - 4*xi*x4 + 4*xi*x5 + 4*(zeta-eta)*x6
        
        crr = cross(exi,eeta)
        J=norm(crr)
        global f = -dot(y-x,nx)/norm(y-x)^3 * J
    end
end


#using Debug
function NormalFieldGaussian(points,faces,rfaces,hmag,H0; NP=3, normals=nothing)
    
    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    
    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)
    B = zeros(Float64,N)

    xw = gaussianpoints(NP)
        
    ### This one is must to know
        
    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]
        
        for ti in 1:size(faces,2)
            v1,v2,v3 = faces[:,ti]
            !(v1==xkey || v2==xkey || v3==xkey) || continue
            
            x1 = points[:,v1]
            x2 = points[:,v2]
            x3 = points[:,v3]

            x4 = points[:,rfaces[3,4*ti]]
            x5 = points[:,rfaces[1,4*ti]]
            x6 = points[:,rfaces[2,4*ti]]

            f1,f2,f3 = 0,0,0
            
            for j in 1:NP
                xi,eta,w = xw[j,1],xw[j,2],xw[j,3]
                @calculationf
                
                f1 += w/2*f*(1-xi-eta)
                f2 += w/2*f*xi
                f3 += w/2*f*eta
            end

            A[v1,xkey] += f1
            A[v2,xkey] += f2 
            A[v3,xkey] += f3 
        end

    end

    ### And here comes a singular integrals

    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]

        #println("xkey is $xkey")
        
        for ti in FaceVRing(xkey,faces)

            w = faces[:,ti] .== xkey
            cw = w[[3,1,2]]
            ccw = w[[2,3,1]]
                        
            x1 = points[:,faces[w,ti]...]
            x2 = points[:,faces[cw,ti]...]
            x3 = points[:,faces[ccw,ti]...]

            x4 = points[:,rfaces[ccw,4*ti]...]
            x5 = points[:,rfaces[w,4*ti]...]
            x6 = points[:,rfaces[cw,4*ti]...]
                        
            s = zeros(3)

            
            ds = Inf
            ksi0 = 1
            eta0 = 1
            i = 1
            for i in 1:10  ### Now for simplicity
                ksi0 /=2
                eta0 /=2

                s1 = s
                
                for j in 1:NP
                    xi,eta = ksi0 + ksi0*xw[j,1],0. + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = 0. + ksi0*xw[j,1],eta0 + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = ksi0 - ksi0*xw[j,1],eta0 - eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                s2 = s

                ds_new = norm(s2-s1)
                if ds_new<ds
                    ds = ds_new
                else
                    break
                end
                    
            end

            A[faces[w,ti]...,xkey] += s[1]
            A[faces[cw,ti]...,xkey] += s[2]
            A[faces[cw,ti]...,xkey] += s[3]
        end
    end

    B = Float64[dot(H0,normals[:,xkey]) for xkey in 1:N]
    
    B = B*2/(hmag+1)
    A = A*(hmag-1)/(hmag+1)/2/pi
    A = transpose(A)
    Hn = (eye(A) - A)\B

    return Hn
    
end



function NormalFieldTrapezodial(points,faces,rfaces,hmag,H0; NP=3, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)
    B = zeros(Float64,N)

    xw = gaussianpoints(NP)

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]
        for ykey in 1:N
            !(xkey==ykey) || continue
            
            y = points[:,ykey]
            
            A[ykey,xkey] = -dot(y-x,nx)/norm(y-x)^3 * vareas[ykey]
        end

        ### Making a proper hole
        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2
            f2 = -dot(points[:,v2]-x,nx)/norm(points[:,v2]-x)^3
            f3 = -dot(points[:,v3]-x,nx)/norm(points[:,v3]-x)^3
            A[v2,xkey] -= f2*area/3
            A[v3,xkey] -= f3*area/3
        end
    end

    ### And here comes a singular integrals

    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]

        #println("xkey is $xkey")
        
        for ti in FaceVRing(xkey,faces)

            w = faces[:,ti] .== xkey
            cw = w[[3,1,2]]
            ccw = w[[2,3,1]]
                        
            x1 = points[:,faces[w,ti]...]
            x2 = points[:,faces[cw,ti]...]
            x3 = points[:,faces[ccw,ti]...]

            x4 = points[:,rfaces[ccw,4*ti]...]
            x5 = points[:,rfaces[w,4*ti]...]
            x6 = points[:,rfaces[cw,4*ti]...]
                        
            s = zeros(3)

            
            ds = Inf
            ksi0 = 1
            eta0 = 1
            i = 1
            for i in 1:10  ### Now for simplicity
                ksi0 /=2
                eta0 /=2

                s1 = s
                
                for j in 1:NP
                    xi,eta = ksi0 + ksi0*xw[j,1],0. + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = 0. + ksi0*xw[j,1],eta0 + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = ksi0 - ksi0*xw[j,1],eta0 - eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationf
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                s2 = s

                ds_new = norm(s2-s1)
                if ds_new<ds
                    ds = ds_new
                else
                    break
                end
                    
            end

            A[faces[w,ti]...,xkey] += s[1]
            A[faces[cw,ti]...,xkey] += s[2]
            A[faces[cw,ti]...,xkey] += s[3]
        end
    end

    B = Float64[dot(H0,normals[:,xkey]) for xkey in 1:N]
    
    B = B*2/(hmag+1)
    A = A*(hmag-1)/(hmag+1)/2/pi
    A = transpose(A)
    Hn = (eye(A) - A)\B

    return Hn
    
end


macro calculationfpot()
    return quote
        zeta = 1-xi-eta

        y = zeta*(2*zeta-1)*x1 + xi*(2*xi-1)*x2 + eta*(2*eta-1)*x3 + 4*xi*zeta*x4 + 4*xi*eta*x5 + 4*eta*zeta*x6
        exi = (1-4*zeta)*x1 + (4*xi-1)*x2 + 0*x3 + 4*(zeta-xi)*x4 + 4*eta*x5 - 4*eta*x6
        eeta = (1-4*zeta)*x1 + 0*x2 + (4*eta-1)*x3 - 4*xi*x4 + 4*xi*x5 + 4*(zeta-eta)*x6
        
        crr = cross(exi,eeta)
        J=norm(crr)
        ny = crr/J   ### I can also check signatures
        global f = dot(y-x,ny)/norm(y-x)^3 * J
    end
end

function PotentialGaussian(points,faces,rfaces,hmag,H0; NP=3, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)
    B = zeros(Float64,N)

    xw = gaussianpoints(NP)
        
    ### This one is must to know
        
    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]
        
        for ti in 1:size(faces,2)
            v1,v2,v3 = faces[:,ti]
            !(v1==xkey || v2==xkey || v3==xkey) || continue
            
            x1 = points[:,v1]
            x2 = points[:,v2]
            x3 = points[:,v3]

            x4 = points[:,rfaces[3,4*ti]]
            x5 = points[:,rfaces[1,4*ti]]
            x6 = points[:,rfaces[2,4*ti]]

            f1,f2,f3 = 0,0,0
            
            for j in 1:NP
                xi,eta,w = xw[j,1],xw[j,2],xw[j,3]
                @calculationfpot
                
                f1 += w/2*f*(1-xi-eta)
                f2 += w/2*f*xi
                f3 += w/2*f*eta
            end

            A[v1,xkey] += f1
            A[v2,xkey] += f2 
            A[v3,xkey] += f3 
        end

    end

    ### And here comes a singular integrals

    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]

        #println("xkey is $xkey")
        
        for ti in FaceVRing(xkey,faces)

            w = faces[:,ti] .== xkey
            cw = w[[3,1,2]]
            ccw = w[[2,3,1]]
                        
            x1 = points[:,faces[w,ti]...]
            x2 = points[:,faces[cw,ti]...]
            x3 = points[:,faces[ccw,ti]...]

            x4 = points[:,rfaces[ccw,4*ti]...]
            x5 = points[:,rfaces[w,4*ti]...]
            x6 = points[:,rfaces[cw,4*ti]...]
                        
            s = zeros(3)

            
            ds = Inf
            ksi0 = 1
            eta0 = 1
            i = 1
            for i in 1:10  ### Now for simplicity
                ksi0 /=2
                eta0 /=2

                s1 = s
                
                for j in 1:NP
                    xi,eta = ksi0 + ksi0*xw[j,1],0. + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = 0. + ksi0*xw[j,1],eta0 + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = ksi0 - ksi0*xw[j,1],eta0 - eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                s2 = s

                ds_new = norm(s2-s1)
                if ds_new<ds
                    ds = ds_new
                else
                    break
                end
                    
            end

            A[faces[w,ti]...,xkey] += s[1]
            A[faces[cw,ti]...,xkey] += s[2]
            A[faces[cw,ti]...,xkey] += s[3]
        end
    end

    B = Float64[dot(H0,points[:,xkey]) for xkey in 1:N]
    
    B = B*2/(hmag+1)
    A = A*(hmag-1)/(hmag+1)/2/pi
    A = transpose(A)
    psi = (eye(A) - A)\B

    return psi
end


function PotentialGaussianPozikridis(points,faces,rfaces,hmag,H0; NP=3, regularize=true, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)
    B = zeros(Float64,N)

    xw = gaussianpoints(NP)
        
    ### This one is must to know
        
    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]
        
        for ti in 1:size(faces,2)
            v1,v2,v3 = faces[:,ti]
            #!(v1==xkey || v2==xkey || v3==xkey) || continue
            
            x1 = points[:,v1]
            x2 = points[:,v2]
            x3 = points[:,v3]

            x4 = points[:,rfaces[3,4*ti]]
            x5 = points[:,rfaces[1,4*ti]]
            x6 = points[:,rfaces[2,4*ti]]

            f1,f2,f3 = 0,0,0
            
            for j in 1:NP
                xi,eta,w = xw[j,1],xw[j,2],xw[j,3]
                @calculationfpot
                
                f1 += w/2*f*(1-xi-eta)
                f2 += w/2*f*xi
                f3 += w/2*f*eta
            end

            A[v1,xkey] += f1
            A[v2,xkey] += f2 
            A[v3,xkey] += f3 
        end

    end
    
    ### And here comes a singular integrals

    B = Float64[dot(H0,points[:,xkey]) for xkey in 1:N]

    if regularize==true
        B = B*2/(hmag+1)
        A = A - diagm(diag(A))
        A = transpose(A)
        psi = (eye(A)*(1- (hmag-1)/(hmag+1)) - 1/2/pi * (hmag-1)/(hmag+1) * (A - diagm(Float64[sum(A[i,:]) for i in 1:size(A,2)]))) \ B
    else
        B = B*2/(hmag+1)
        A = A*(hmag-1)/(hmag+1)/2/pi
        A = transpose(A)
        psi = (eye(A) - A)\B
    end
    
    return psi
    
end


function PotentialGaussianTrapezodial(points,faces,rfaces,hmag,H0; NP=3, normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end
    
    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)
    B = zeros(Float64,N)

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end

    for xkey in 1:N
        x = points[:,xkey]
        for ykey in 1:N
            !(xkey==ykey) || continue
            
            y = points[:,ykey]
            ny = normals[:,ykey]
            A[ykey,xkey] = dot(y-x,ny)/norm(y-x)^3 * vareas[ykey]
        end

        ### Making a proper hole
        for (v2,v3) in DoubleVertexVRing(xkey,faces)
            area = norm(cross(points[:,v2]-x,points[:,v3]-x))/2
            f2 = dot(points[:,v2]-x,normals[:,v2])/norm(points[:,v2]-x)^3
            f3 = dot(points[:,v3]-x,normals[:,v3])/norm(points[:,v3]-x)^3
            A[v2,xkey] -= f2*area/3
            A[v3,xkey] -= f3*area/3
        end
    end
    
    xw = gaussianpoints(NP)
    
    for xkey in 1:N
        x = points[:,xkey]
        nx = normals[:,xkey]

        #println("xkey is $xkey")
        
        for ti in FaceVRing(xkey,faces)

            w = faces[:,ti] .== xkey
            cw = w[[3,1,2]]
            ccw = w[[2,3,1]]
                        
            x1 = points[:,faces[w,ti]...]
            x2 = points[:,faces[cw,ti]...]
            x3 = points[:,faces[ccw,ti]...]

            x4 = points[:,rfaces[ccw,4*ti]...]
            x5 = points[:,rfaces[w,4*ti]...]
            x6 = points[:,rfaces[cw,4*ti]...]
                        
            s = zeros(3)

            
            ds = Inf
            ksi0 = 1
            eta0 = 1
            i = 1
            for i in 1:10  ### Now for simplicity
                ksi0 /=2
                eta0 /=2

                s1 = s
                
                for j in 1:NP
                    xi,eta = ksi0 + ksi0*xw[j,1],0. + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = 0. + ksi0*xw[j,1],eta0 + eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                for j in 1:NP
                    xi,eta = ksi0 - ksi0*xw[j,1],eta0 - eta0*xw[j,2]
                    local w = xw[j,3]
                    @calculationfpot
                    s += w/2*ksi0*eta0*f*[1-xi-eta,xi,eta]                    
                end

                s2 = s

                ds_new = norm(s2-s1)
                if ds_new<ds
                    ds = ds_new
                else
                    break
                end
                
            end

            A[faces[w,ti]...,xkey] += s[1]
            A[faces[cw,ti]...,xkey] += s[2]
            A[faces[cw,ti]...,xkey] += s[3]
        end
    end

    B = Float64[dot(H0,points[:,xkey]) for xkey in 1:N]
    
    B = B*2/(hmag+1)
    A = A*(hmag-1)/(hmag+1)/2/pi
    A = transpose(A)
    psi = (eye(A) - A)\B

    return psi
    
end


function NormalFieldCurved(points,faces,rfaces,hmag,H0; normals=nothing)

    if normals==nothing
        normals = Array{Float64}(undef,size(points)...)
        NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    end

    N = maximum(faces)
    
    # normals = Array(Float64,size(points)...)
    # NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))    
    
    A = zeros(Float64,N,N)

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end


    for xkey in 1:N

        x = points[:,xkey]
        nx = normals[:,xkey]

        for ykey in 1:N
            if xkey==ykey   ### This is where all magic will begin

                for ti in FaceVRing(xkey,faces)

                    w = faces[:,ti].==xkey
                    cw = w[[3,1,2]]
                    ccw = w[[2,3,1]]
                    
                    x1 = points[:,faces[w,ti]...]
                    x2 = points[:,faces[cw,ti]...]
                    x3 = points[:,faces[ccw,ti]...]

                    x4 = points[:,rfaces[ccw,4*ti]...]
                    x5 = points[:,rfaces[w,4*ti]...]
                    x6 = points[:,rfaces[cw,4*ti]...]

                    #x4,x5 = x5,x4
                    #x1,x2 = x2,x1

                    xx(xi,eta) = (zeta=1-xi-eta; zeta*(2*zeta-1)*x1 + xi*(2*xi-1)*x2 + eta*(2*eta-1)*x3 + 4*xi*zeta*x4 + 4*xi*eta*x5 + 4*eta*zeta*x6)
                    
                    exi(xi,eta) = (zeta=1-xi-eta; (1-4*zeta)*x1 + (4*xi-1)*x2 + 0*x3 + 4*(zeta-xi)*x4 + 4*eta*x5 - 4*eta*x6)
                    eeta(xi,eta) = (zeta=1-xi-eta; (1-4*zeta)*x1 + 0*x2 + (4*eta-1)*x3 - 4*xi*x4 + 4*xi*x5 + 4*(zeta-eta)*x6)

                    
                    f(xi,eta) = (crr = cross(exi(xi,eta),eeta(xi,eta));J=norm(crr);ny=crr/J;y=xx(xi,eta);val = -dot(y-x,nx)/norm(y-x)^3 * J ; [val*(1-xi-eta),val*xi,val*eta])
                    #f(xi,eta) = (crr = cross(exi(xi,eta),eeta(xi,eta));J=norm(crr);ny=crr/J;y=xx(xi,eta); -dot(y-x,ny)/norm(y-x)^3 * J * (1-xi-eta))
                    #f(xi,eta) = (crr = cross(exi(xi,eta),eeta(xi,eta));J=norm(crr);ny=crr/J;y=xx(xi,eta); -dot(y-x,ny)/norm(y-x)^3 * J )
#dot(y-x,nx)/norm(y-x)^3 * vareas[ykey]

                    
                    s = zeros(3)
                    
                    ksi0 = 1
                    eta0 = 1
                    for i in 1:10  ### Now for simplicity
                        ksi0 /=2
                        eta0 /=2
                        
                        s += trquad((xi,eta) -> ksi0*eta0*f(ksi0 + ksi0*xi,0 + eta0*eta))
                        s += trquad((xi,eta) -> ksi0*eta0*f(0. + ksi0*xi,eta0 + eta0*eta))
                        s += trquad((xi,eta) -> ksi0*eta0*f(ksi0 - ksi0*xi,eta0 - eta0*eta))
                    end
                    
                    A[ykey,xkey] += s[1]
                    A[faces[cw,ti]...,xkey] += s[2]
                    A[faces[ccw,ti]...,xkey] += s[3]

                    # # ### I also need to remove

                    # vcw, = faces[cw,ti]
                    # vccw, = faces[ccw,ti]
                    
                    # S = norm(cross(points[:,vcw]-x,points[:,vccw]-x))/2
                    
                    # A[vcw,xkey] -= -dot(points[:,vcw]-x,normals[:,vcw])/norm(points[:,vcw]-x)^3 * S/3
                    # A[vccw,xkey] -= -dot(points[:,vccw]-x,normals[:,vccw])/norm(points[:,vccw]-x)^3 * S/3
                end

                #A[ykey,xkey] = s
                continue
            end

            y = points[:,ykey]
            #ny = normals[:,ykey]
            
            A[ykey,xkey] = -dot(y-x,nx)/norm(y-x)^3 * vareas[ykey]
        end
    end

    B = zeros(Float64,N)
    for xkey in 1:N
        B[xkey] = 2*dot(H0,normals[:,xkey])/(hmag+1)
    end

    A = A'

    #psi = (eye(A)*(1-1/4/pi*(hmag-1)/(hmag+1)) - 1/2/pi * (hmag-1)/(hmag+1) * (A - diagm(Float64[sum(A[i,:]) for i in 1:size(A,2)]))) \ B
    Hn = (eye(A) - (hmag-1)/(hmag+1)/2/pi*A)\B

    #psi = (eye(A) - (hmag-1)/(hmag+1)/2/pi * A) \ B

    return Hn,A
end
