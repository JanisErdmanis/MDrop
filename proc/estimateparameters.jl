ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"
using Storage
using SurfaceGeometry

function rotatetoaxis(points)

    r = 0
    v1 = nothing
    v2 = nothing
    for xkey in 1:size(points,2)
        x = points[:,xkey]
        for ykey in 1:size(points,2)
            y = points[:,ykey]
            rr = norm(x-y)
            if rr>r
                v1 = xkey
                v2 = ykey
                r = rr
            end
        end
    end

    l1 = points[:,v1] - points[:,v2]
    nl1 = l1/norm(l1)

    n = [0,0,1] + nl1
    n /= norm(n)

    Lx = [0 0 0; 0 0 -1; 0 1 0]
    Ly = [0 0 1; 0 0 0; -1 0 0]
    Lz = [0 -1 0; 1 0 0; 0 0 0]

    Ln = n[1]*Lx + n[2]*Ly + n[3]*Lz
    R = expm(pi*Ln)

    pointsr = zeros(points)

    for xkey in 1:size(points,2)
        pointsr[:,xkey] = R*points[:,xkey]
    end

    r2 = 0
    v3 = nothing
    for xkey in 1:size(pointsr,2)
        x = pointsr[:,xkey]
        rr = sqrt(x[1]^2 + x[2]^2)
        if rr>r2
            v3 = xkey
            r2 = rr
        end
    end

    l2 = pointsr[:,v3]
    nl2 = l2/norm(l2)
    alpha = atan2(nl2[2],nl2[1])

    R2 = expm(-alpha*Lz)

    return R2*R
end
    #pointsrr = zeros(points)


# a,b,c = 2,1/4,1/4
# (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

# case = 55
# session = "rotatingfield"
# calculate = false

# include("diaryrot.jl")

# tmax,points,faces = memory[end]
function getabc(points)
    # R = rotatetoaxis(points)
    # for xkey in 1:size(points,2)
    #     points[:,xkey] = R*points[:,xkey]
    # end

    ar = 0
    al = 0
    br = 0
    bl = 0
    cr = 0
    cl = 0
    
    for xkey in 1:size(points,2)
        x = points[:,xkey]
        x[1]<al || (al=x[1])
        x[1]>ar || (ar=x[1])
        x[2]<bl || (bl=x[2])
        x[2]>br || (br=x[2])
        x[3]<cl || (cl=x[3])
        x[3]>cr || (cr=x[3])
    end

    a = (ar - al)/2
    b = (br - bl)/2
    c = (cr - cl)/2

    return a,b,c
end




#### Old try of getting a projection ####
# x = points[1,:][:]
# y = points[2,:][:]
# z = points[3,:][:]

# using PyPlot


# tri = Any[]

# for ti in 1:size(faces,2)
#     push!(tri,faces[:,ti]-1)
# end
# #push!(tri,[0,1,2])
# #push!(tri,[90,2,200])

# # x = x[1:100:end]
# # y = y[1:100:end]

# #triplot(x,y,tri)
# #triplot(z,x,tri)
# #colors = 
# tripcolor(z,x,tri)
