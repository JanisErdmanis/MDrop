MeshData = dirname(dirname(@__FILE__)) * "/meshes/"

using JLD
using SurfaceGeometry

function loadmesh(fname)
    try
        d = load(MeshData*fname)
        points = d["points"]
        faces = d["faces"]
        return (points,faces)
    catch
        error("No valid mesh file found with name [preproc/$fname]")
    end
end

ellipsoidhash_matlab(a,b,c,step) = "$a$b$c-$step"*"-matlab.jld"
ellipsoidhash_cgal(a,b,c,args...) = "$a$b$c-$(args...)"*"-cgal.jld"

function EllipsoidMeshLoad(a,b,c,args...;method="matlab")
    if method=="matlab"
        fname = ellipsoidhash_matlab(a,b,c,args...)
        try
            p,t = loadmesh(fname)
        catch
            mesher = DistmeshSurfaceMesher(args...,[-1. -1 -1.; 1. 1. 1.]*1.5)
            points, faces = EllipsoidMesh(a,b,c,mesher)
            save(MeshData*fname,"points",points,"faces",faces)
            #save("ellipsoid.jld",hash*":p",p,hash*"points",t)
        end
    elseif method=="cgal"
        fname = ellipsoidhash_cgal(a,b,c,args...)
        p,t = loadmesh(fname)
    else
        println("Method not valid")
        throw(KeyError)
        #return 0 ### Something better should be implemented
    end

    return loadmesh(fname)
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

#SimulationData = homedir()*"/SimulationData" 

# using JLD
# points,faces = starmesh(1.5,0.5,0.1)
# @save "meshes/star.jld" points faces


