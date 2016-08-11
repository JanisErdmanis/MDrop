#using Cubature
# function InterfaceSpeed(points,faces,tensorn,etaP)

#     ### Adapting directly the algorithm from strquad
#     q(t,hi,x,normalx,fn,y1,y2,y3) = begin

#         normaly = cross(y2-y1,y3-y1)
#         normaly /= norm(normaly)

#         R = 1/(cos(hi) + sin(hi))
#         xi = t*R*cos(hi)
#         eta = t*R*sin(hi)

#         y = y1*(1 - xi - eta) + y2*xi + y3*eta

#         qq = 1./8/pi/etaP*(dot(normalx,normaly) + dot(normalx,x-y)*dot(normaly,x-y)/norm(x-y)^2)

#         B = dot(y3-y1,y2-y1)/norm(y2-y1)^2
#         C = norm(y3-y1)^2/norm(y2-y1)^2
#         hS = norm(cross(y2-y1,y3-y1))
#         R*qq / sqrt(cos(hi)^2 + B*sin(2*hi) + C*sin(hi)^2)*hS/norm(y2-y1) * fn
#     end

#     ### A possible improvement would be the interpolation of curvature and normal
#     f(xi,eta,x,normalx,fn,y1,y2,y3) = begin

#         normaly = cross(y2-y1,y3-y1)
#         crr = norm(normaly)
#         normaly /= crr
         
#         y = y1*(1 - xi - eta) + y2*xi + y3*eta

#         1./8/pi/etaP*(dot(normalx,normaly)/norm(x-y) + dot(normalx,x -y)*dot(normaly,x-y)/norm(x-y)^3) * fn * crr
#     end

#     normals = Array(Float64,size(points)...)
#     NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))
    
#     velocityn = zeros(Float64,size(points,2))
#     for xkey in 1:size(points,2)
                        
#         normalx = normals[:,xkey]  #CollocactionNormal(xkey,faces,enodes,points) ###
#         x = points[:,xkey]

#         s = 0 
#         for ykey in 1:size(faces,2)

#             ytri = faces[:,ykey]

#             if xkey in ytri ### This one is questionable
#                 ### Singular integration
#                 w = xkey.==ytri

#                 # So I need checking here

#                 vy1, = ytri[w[[1,2,3]]]
#                 vy2, = ytri[w[[3,1,2]]]
#                 vy3, = ytri[w[[2,3,1]]]
                
#                 y1 = points[:,vy1]
#                 y2 = points[:,vy2]
#                 y3 = points[:,vy3]


#                 ds = hcubature(par-> q(par[1],par[2],x,normalx,tensorn[ykey],y1,y2,y3),[0,0],[1,pi/2],reltol=1e-3)
                                
#                 s += ds[1]
                
#             else
#                 ### A simple one

#                 y1 = points[:,ytri[1]]
#                 y2 = points[:,ytri[2]]
#                 y3 = points[:,ytri[3]]

#                 ftilde = (xi,eta) -> f(xi,eta,x,normalx,tensorn[ykey],y1,y2,y3)
#                 ds = UnitTriangleIntegration(ftilde,NP=3) #

#                 s += ds[1]
#             end
#         end

#         velocityn[xkey] = s
#     end

#     return velocityn
# end


function InterfaceSpeedPozikridis(points,faces,forcen,etaP)

    ### A possible improvement would be the interpolation of curvature and normal
    normals = Array(Float64,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    
    #phi = Array(Float64,size(points,2))
    velocityn = zeros(Float64,size(points,2))
    
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]
        fx = forcen[xkey]
                
        s = 0
        for ykey in 1:size(points,2)
            if ykey==xkey
                continue
            end

            y = points[:,ykey]
            ny = normals[:,ykey]
            fy = forcen[ykey]
            
            s += vareas[ykey]*1./8/pi/etaP* ( dot(nx,ny)/norm(x-y) + dot(nx,x -y)*dot(ny,x-y)/norm(x-y)^3 )*(fy - fx)
        end

        velocityn[xkey] = s
    end

    return velocityn
end


function InterfaceSpeedZinchenko(points,faces,forcen,etaP,gammap)

    ### A possible improvement would be the interpolation of curvature and normal
    normals = Array(Float64,size(points)...)
    NormalVectors!(normals,points,faces,i->FaceVRing(i,faces))

    vareas = zeros(Float64,size(points,2))
    for i in 1:size(faces,2)
        v1,v2,v3 = faces[:,i]
        area = norm(cross(points[:,v2]-points[:,v1],points[:,v3]-points[:,v1]))/2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    
    #phi = Array(Float64,size(points,2))
    velocityn = zeros(Float64,size(points,2))
    
    for xkey in 1:size(points,2)

        x = points[:,xkey]
        nx = normals[:,xkey]
        fx = forcen[xkey]
                
        s = 0
        for ykey in 1:size(points,2)
            if ykey==xkey
                continue
            end

            y = points[:,ykey]
            ny = normals[:,ykey]
            fy = forcen[ykey]

            ### I will need to check a missing 2
            s += vareas[ykey]*1./8/pi/etaP* dot(y-x,nx+ny)/norm(y-x)^3*(1-3*dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^2) * gammap

            s += vareas[ykey]*1./8/pi/etaP* ( dot(nx,ny)/norm(x-y) + dot(nx,x -y)*dot(ny,x-y)/norm(x-y)^3 )*(fy - fx)
        end

        velocityn[xkey] = s
    end

    return velocityn
end


function Velocity3D(cmsh,etaP,gammaP)

    curvaturep = Array(Float64,size(cmsh.points,2))
    tensorn = Array(Float64,size(cmsh.faces,2))
    
    for xkey in 1:size(cmsh.points,2)
        curvaturep[xkey] = vcurvature(xkey,cmsh)
    end
    
    for xkey in 1:size(cmsh.faces,2)
        v1,v2,v3 = cmsh.faces[:,xkey]
        tensorn[xkey] = -(curvaturep[v1] + curvaturep[v2] + curvaturep[v3])/3*2*gammaP
    end

    velocityn = InterfaceSpeed(cmsh,tensorn,etaP)
    velocity = Array(Float64,size(cmsh.points))
    for vi in 1:size(cmsh.points,2)

        normall = vnormal(vi,cmsh)
        velocity[:,vi] = vnormal(vi,cmsh)*velocityn[vi]
    end

    return velocity
end

function MagneticVelocity3D(cmsh,hmag,H0,etaP,gammaP)

    curvaturep = Array(Float64,size(cmsh.points,2))
    tensorn = Array(Float64,size(cmsh.faces,2))

    println("Curvature...")
    for xkey in 1:size(cmsh.points,2)
        curvaturep[xkey] = vcurvature(xkey,cmsh)
    end
    ### One should tjke in account of usual multiplicaction by 2

    ### Now I also need calculation for field
    #smooth_surface_edge_nodes!(cmsh) ### The one for which I fighted for

    ### Probably for this I could write a macro
    points = cmsh.points
    faces = cmsh.faces
    points = cmsh.points
    neighs = cmsh.neighs
    enodes = cmsh.enodes
    epoints = cmsh.epoints
    


    println("Calculation of both field components...")
    Hn, Ht = SurfaceFields(cmsh,hmag,H0)

    #################### Now construction for a field tensor ############
    println("and velocity...")
    
    for xkey in 1:size(cmsh.faces,2)
        v1,v2,v3 = cmsh.faces[:,xkey]
        tensorn[xkey] = -(curvaturep[v1] + curvaturep[v2] + curvaturep[v3])/3*2*gammaP + hmag*(hmag-1)/8/pi*Hn[xkey]^2 + (hmag-1)/8/pi*Ht[xkey]^2
    end

    velocityn = InterfaceSpeed(cmsh,tensorn,etaP)
    velocity = Array(Float64,size(cmsh.points))
    for vi in 1:size(cmsh.points,2)

        normall = vnormal(vi,cmsh)
        
        velocity[:,vi] = vnormal(vi,cmsh)*velocityn[vi]

    end

    return velocity
end
