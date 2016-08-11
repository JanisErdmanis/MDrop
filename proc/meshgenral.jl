ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

immutable DistmeshSurfaceMesher
    step::Float64
end

function DistmeshSurfaceMesher(;step=0.2)
    DistmeshSurfaceMesher(step)
end

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


function SurfaceMesh(fdis::AbstractString,ds::DistmeshSurfaceMesher)


    eval(:(using MATLAB))

    step = float(ds.step)

    ### A function definition here

    tmpdir = mktempdir()

    eval_string("""addpath('$(Pkg.dir("SurfaceGeometry","src","libraries","distmesh"))')""")
    eval_string("""addpath('$tmpdir')""")

    #cd(tmpdir)
    
    file = open("meshgen.m","w")

    boxsize = [-2.1 -1.1 -1.6; 2.1 1.1 1.6]
    
    code = """
    function [p,t] = meshgen() %,boxsize)
        $fdis
        
        [p,t]=distmeshsurface(@fd,@huniform,$step,$boxsize);
    end

    function h=huniform(p,varargin)
        h=ones(size(p,1),1);
    end
    
    """
    write(file,code)
    close(file)
    
    p, t = mxcall(:meshgen,2)
    # finally
    #     cd(olddir)
    # end
    
    
    t = transpose(map(Int,t))
    p = transpose(p)

    ### Converting to tuples
    ### Checking if there is need to change order for triangles    


    if !isoriented(p,t)
        t[1,:], t[2,:] = t[2,:], t[1,:]
    end

    p,t = map(Float64,p), map(Int,t)
    
    for vi in 1:size(p)[end]
        p[:,vi] = pushback(x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1,p[:,vi])
    end
    
    return p,t
end



function EllipsoidMesh(a::Real,b::Real,c::Real,ds::DistmeshSurfaceMesher)

    
    fdis = """
        function f = fd(p)
            f = p(:,1).^2/$(a^2)+p(:,2).^2/$(b^2)+p(:,3).^2/$(c^2)-1;
        end
    
    """

    SurfaceMesh(fdis,ds)    
end


# function EllipsoidMeshOld(a::Real,b::Real,c::Real,ds::DistmeshSurfaceMesher)


#     eval(:(using MATLAB))

#     a = float(a)
#     b = float(b)
#     c = float(c)
#     step = float(ds.step)

#     ### A function definition here

#     tmpdir = mktempdir()

#     mat"""addpath($(Pkg.dir("SurfaceGeometry","src","libraries","distmesh")))"""
#     mat"""addpath($tmpdir)"""

#     #cd(tmpdir)
    
#     file = open("meshgen.m","w")
    
#     code = """
#     function [p,t] = meshgen(a,b,c,step) %,boxsize)
#         function f = fd(p)
#             f = p(:,1).^2/a^2+p(:,2).^2/b^2+p(:,3).^2/c^2-1;
#         end
        
#         [p,t]=distmeshsurface(@fd,@huniform,step,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
#     end

#     function h=huniform(p,varargin)
#         h=ones(size(p,1),1);
#     end
    
#     """
#     write(file,code)
#     close(file)
    
#     p, t = mxcall(:meshgen,2,a,b,c,step)
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





EllipsoidMesh(1.,1.,1.,DistmeshSurfaceMesher(0.2))
