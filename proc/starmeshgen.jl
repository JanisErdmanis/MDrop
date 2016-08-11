ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using SurfaceGeometry

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

using JLD
points,faces = starmesh(1.5,0.5,0.1)
@save "meshes/star.jld" points faces


