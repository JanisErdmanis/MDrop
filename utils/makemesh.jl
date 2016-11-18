# immutable DistmeshSurfaceMesher
#     step::Float64
# end

# function DistmeshSurfaceMesher(;step=0.2)
#     DistmeshSurfaceMesher(step)
# end

# function EllipsoidMesh(a::Real,b::Real,c::Real,ds::DistmeshSurfaceMesher)

#     olddir = pwd()
#     try
#         cd(Pkg.dir("SurfaceGeometry","src","libraries","distmesh"))
#         Base.syntax_deprecation_warnings(false)
#         eval(:(using MATLAB))

#         a = float(a)
#         b = float(b)
#         c = float(c)
#         step = float(ds.step)
#         p, t = mxcall(:elipsoid,2,a,b,c,step)
#     finally
#         cd(olddir)
#     end
        
    
#     t = transpose(map(Int,t))
#     p = transpose(p)

#     ### Converting to tuples
#     ### Checking if there is need to change order for triangles    


#     if !isoriented(p,t)
#         t[1,:], t[2,:] = t[2,:], t[1,:]
#     end

#     p,t = map(Float64,p), map(Int,t)
    
#     for vi in 1:size(p)[end]
#         p[:,vi] = pushback(x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1,p[:,vi])
#     end
    
#     return p,t
# end


# function SurfaceMesh(fdis::AbstractString,ds::DistmeshSurfaceMesher;boxsize=[-1. -1. -1.; 1. 1. 1.])


#     eval(:(using MATLAB))

#     step = float(ds.step)

#     ### A function definition here

#     tmpdir = mktempdir()

#     eval_string("""addpath('$(Pkg.dir("SurfaceGeometry","src","libraries","distmesh"))')""")
#     eval_string("""addpath('$tmpdir')""")

#     #cd(tmpdir)
    
#     file = open("meshgen.m","w")

#     #boxsize = [-2.1 -1.1 -1.6; 2.1 1.1 1.6]
    
#     code = """
#     function [p,t] = meshgen() %,boxsize)
#         $fdis
        
#         [p,t]=distmeshsurface(@fd,@huniform,$step,$boxsize);
#     end

#     function h=huniform(p,varargin)
#         h=ones(size(p,1),1);
#     end
    
#     """
#     write(file,code)
#     close(file)
    
#     p, t = mxcall(:meshgen,2)
#     # finally
#     #     cd(olddir)
#     # end
    
    
#     t = transpose(map(Int,t))
#     p = transpose(p)

#     ### Converting to tuples
#     ### Checking if there is need to change order for triangles    


#     if !isoriented(p,t)
#         t[1,:], t[2,:] = t[2,:], t[1,:]
#     end

#     p,t = map(Float64,p), map(Int,t)
    
#     for vi in 1:size(p)[end]
#         p[:,vi] = pushback(x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1,p[:,vi])
#     end
    
#     return p,t
# end

using SurfaceGeometry


function EllipsoidMesh(a::Real,b::Real,c::Real,ds::DistmeshSurfaceMesher)

    
    fdis = """
        function f = fd(p)
            f = p(:,1).^2/$(a^2)+p(:,2).^2/$(b^2)+p(:,3).^2/$(c^2)-1;
        end
    
    """

    SurfaceMesh(fdis,ds;boxsize=[-2.1 -1.1 -1.6; 2.1 1.1 1.6])    
end

function elmesh4(a,c,step;boxsize = [-2.1 -2 -1.; 2.1 2 1.])

    fdis = """
    function f = fd(p)
        f = (p(:,1).^2 + p(:,2).^2).^2/$(a^4) + p(:,3).^4/$(c^4)-1;
    end
    
    """
    p,t = SurfaceMesh(fdis,DistmeshSurfaceMesher(step,boxsize))

    return p,t
end

function starmesh(a,c,step;boxsize = [-2.1 -2 -1.; 2.1 2 1.],k=3)

    ###  
    fdis = """
    function f = fd(p)

        phi = atan2(p(:,2),p(:,1));
        r0 = cos(3*phi);
        r2 = p(:,1).^2 + p(:,2).^2;
        theta = atan(p(:,3)./sqrt(r2));
        f = (r2 - r0.^2).*cos(theta).^2  + p(:,3).^2/$(c^2) - 1;
    end
    """
    p,t = SurfaceMesh(fdis,DistmeshSurfaceMesher(step,boxsize))

    return p,t
end


using DataFrames
using JLD

function EllipsoidMeshBm(Bm,step)
    
    df1 = readtable("$(dirname(@__FILE__))/cebers/ca_9_short_en_1_200.dat",separator='\t',header=false)
    df2 = readtable("$(dirname(@__FILE__))/cebers/cb_9_short_en_1_200.dat",separator='\t',header=false)

    Bm_ = df1[:x1]
    c_ = (df1[:x2].*df2[:x2]).^(1/3)
    a_ = (c_ ./ df1[:x2])
    b_ = (c_ ./ df2[:x2])

    i = 1
    for i in 1:length(Bm_)
        if Bm<Bm_[i]
            break
        end
    end

    a,b,c = a_[i], b_[i], c_[i]
    info("Generating mesh $i with a=$a b=$b c=$c")

    fdis = """
    function f = fd(p)
        f = p(:,1).^2/$(a^2)+p(:,2).^2/$(b^2)+p(:,3).^2/$(c^2)-1;
    end
    """
    
    #boxsize = [-1.1*a -1.1*b -1.1*c; 1.1*a 1.1*b 1.1*c]
    boxsize = [-1.2*a -1.2*b -1.2*c; 1.2*a 1.2*b 1.2*c]

    points,faces = SurfaceMesh(fdis,DistmeshSurfaceMesher(step,boxsize))
    println("Number of faces N=$(size(faces,2))")

    save("meshes/elips$step-Bm$Bm.jld","points",points,"faces",faces)
end

### For CGAL
# a,b,c=2.9985,0.68425,0.48738
# fdis(x,y,z) = x^2/Float32(a^2) + y^2/Float32(b^2) + z^2/Float32(c^2) - 1
# mesher = CGALSurfaceMesher(10,0.1,0.1,3)
# points,faces = SurfaceMesh(fdis,mesher)
# save("meshes/elipsCGAL-Bm25.jld","points",points,"faces",faces)
