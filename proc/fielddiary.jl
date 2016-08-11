ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

using SurfaceGeometry
using Storage

include("Calculation/field.jl")

#mup = 25
isdefined(:mup) || (mup = 10)
H0 = [1,0,0]

function ellipsoid_demagnetization_coefficients(a,b,c)

    UP_BOUND = 1000

    Ru2(u) = (u+a^2)*(u+b^2)*(u+c^2)

    nx = 1/2 * a*b*c * quadgk(s -> 1/(s+a^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    ny = 1/2 * a*b*c * quadgk(s -> 1/(s+b^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]
    nz = 1/2 * a*b*c * quadgk(s -> 1/(s+c^2)/sqrt(Ru2(s)), 0, UP_BOUND)[1]

    return [nx, ny, nz]
end

function EllipsoidField(a,b,c,mu,H0)

    H0x, H0y, H0z = H0
    nx, ny, nz = ellipsoid_demagnetization_coefficients(a,b,c)

    Hix = H0x/(1 + (mu-1)*nx)
    Hiy = H0y/(1 + (mu-1)*ny)
    Hiz = H0z/(1 + (mu-1)*nz)

    return [Hix,Hiy,Hiz]
end

function EllipsoidField(a,b,c)
    EllipsoidField(a,b,c,mup,H0)
    # H0x, H0y, H0z = H0
    # nx, ny, nz = ellipsoid_demagnetization_coefficients(a,b,c)

    # Hix = H0x/(1 + (mup-1)*nx)
    # Hiy = H0y/(1 + (mup-1)*ny)
    # Hiz = H0z/(1 + (mup-1)*nz)

    # return [Hix,Hiy,Hiz]
end

function ellipsoid_normals(a,b,c,points)
    normals = zeros(Float64,size(points))
    for xkey in 1:size(points,2)
        x,y,z = points[:,xkey]
        nx = [x/a^2,y/b^2,z/c^2]
        nx /= norm(nx)
        normals[:,xkey] = nx
    end
    return normals
end

function theoretical_potential(a,b,c,p)
    normals = ellipsoid_normals(a,b,c,p)
    Htheor = EllipsoidField(a,b,c)

    psit = zeros(Float64,size(p,2))
    for xkey in 1:size(p,2)
        nx = normals[:,xkey]
        psit[xkey] = dot(Htheor,p[:,xkey])
    end
    return psit
end

function theoretical_tfield(a,b,c,p)
    normals = ellipsoid_normals(a,b,c,p)
    Htheor = EllipsoidField(a,b,c)

    Ht = zeros(Float64,size(p))
    for xkey in 1:size(p,2)
        nx = normals[:,xkey]
        Ht[:,xkey] = (eye(3) - nx*nx')*Htheor#dot(Htheor,p[:,xkey])
    end
    return Ht
end

function theoretical_nfield(a,b,c,p)
    normals = ellipsoid_normals(a,b,c,p)
    Htheor = EllipsoidField(a,b,c)

    Hn = zeros(Float64,size(p,2))
    for xkey in 1:size(p,2)
        nx = normals[:,xkey]
        Hn[xkey] = dot(Htheor,nx) #dot(Htheor,p[:,xkey])
    end
    return Hn
end

if !isdefined(:meshcase)
    if isdefined(:ARGS)
        meshcase = parse(Int,ARGS[2])
    else
        meshcase = 5
    end
end

if !isdefined(:case)
    if isdefined(:ARGS)
        case = parse(Int,ARGS[1])
    else
        case = 10
    end
end
#isdefined(:case) || (case=10)

if meshcase==1
    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)

elseif meshcase==2
    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)

elseif meshcase==3
    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.35)

elseif meshcase==4
    a,b,c = 1,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.2)
    points, faces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)

elseif meshcase==5
    a,b,c = 2,1/4,1/4
    # a,b,c = 2,1/6,1/6
    (points,faces)=EllipsoidMeshLoad(a,b,c,0.1)

elseif meshcase==6
    a,b,c = 1/4,1,1
    points,faces = EllipsoidMeshLoad(a,b,c,0.15)
    
elseif meshcase==7
    a,b,c = 2,1/4,1/4
    points,faces = EllipsoidMeshLoad(a,b,c,0.1)
    points, faces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
end

normals = ellipsoid_normals(a,b,c,points)


if case==1
    psi = PotentialSimple(points,faces,mup,H0,regularize=false, normals = normals)

elseif case==2
    psi = PotentialSimple(points,faces,mup,H0,regularize=true, normals=normals)
      
elseif case==3
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    normals = ellipsoid_normals(a,b,c,points)
    psi = PotentialGaussianPozikridis(rp,faces,rfaces,mup,H0;NP=3, normals=normals)

elseif case==4
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    psi = PotentialGaussianTrapezodial(rp,faces,rfaces,mup,H0;NP=3, normals=normals)

elseif case==5
    Hn = NormalField(points,faces,mup,H0; regularize=false, normals=normals)

elseif case==6
    Hn = NormalField(points,faces,mup,H0; regularize=true, normals=normals)

elseif case==7
    psit = theoretical_potential(a,b,c,points)
    Hn = NormalFieldRecalculated(points,faces,psit,normals=normals)

elseif case==8
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    Hn = NormalFieldTrapezodial(rp,faces,rfaces,mup,H0;NP=3,normals=normals)

elseif case==9
    Htt = theoretical_tfield(a,b,c,points)
    Hn = NormalFieldCurrent(points,faces,Htt,mup,H0,normals=normals)

### Combined approach    
elseif case==10
    psi = PotentialSimple(points,faces,mup,H0,regularize=true)
    Ht = HtField(points,faces,psi)
    Hn = NormalFieldCurrent(points,faces,Ht,mup,H0)

elseif case==11
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    Hn = NormalFieldTrapezodial(rp,faces,rfaces,mup,H0;NP=3)

elseif case==12
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    psi = PotentialGaussianPozikridis(rp,faces,rfaces,mup,H0;NP=12, normals=normals)

elseif case==13
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    Hn = NormalFieldGaussian(rp,faces,rfaces,mup,H0; NP=16, normals=normals)

elseif case==14
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    psi = PotentialGaussian(rp,faces,rfaces,mup,H0; NP=3,normals=normals)

### Another combination    
elseif case==15
    psi = PotentialSimple(points,faces,mup,H0,regularize=true)
    rp, rfaces = subdivision(points,faces,x -> x[1]^2/a^2 + x[2]^2/b^2 + x[3]^2/c^2 - 1)
    Hn = NormalFieldTrapezodial(rp,faces,rfaces,mup,H0;NP=3,normals=normals)
end


if isdefined(:psi)
    Ht = HtField(points,faces,psi,normals=normals)

    rpsi = Float64[]
    rHt = Float64[]

    psit = theoretical_potential(a,b,c,points)
    Htt = theoretical_tfield(a,b,c,points)

    Htheor = EllipsoidField(a,b,c)
        
    for xkey in 1:size(points,2)
        push!(rpsi,abs((psi[xkey] - psit[xkey])/psit[xkey]))        
    end

    for xkey in 1:size(points,2)
        push!(rHt,norm(Ht[:,xkey] - Htt[:,xkey])/norm(Htt[:,xkey]))
    end

    rpsi *= 100
    rHt *= 100

    println("mean=$(round(mean(rpsi),4)) mean2=$(round(sum(abs(psi - psit))/sum(abs(psit))*100,4)) q1/4=$(round(quantile(rpsi,1/4),4)) q1/2=$(round(quantile(rpsi,1/2),4)) q3/4=$(round(quantile(rpsi,3/4),4))")
    println("mean=$(round(mean(rHt),4)) mean2=$(round(sum(abs(Ht-Htt))/sum(abs(Htt))*100,4)) q1/4=$(round(quantile(rHt,1/4),4)) q1/2=$(round(quantile(rHt,1/2),4)) q3/4=$(round(quantile(rHt,3/4),4))")
    
end

if isdefined(:Hn)
    rHn = Float64[]
    Htn = theoretical_nfield(a,b,c,points)

    for xkey in 1:size(points,2)
        push!(rHn,abs((Hn[xkey] - Htn[xkey])/Htn[xkey]))
    end

    rHn *= 100
    println("mean=$(round(mean(rHn),4)) mean2=$(round(sum(abs(Hn-Htn))/sum(abs(Htn))*100,4)) q1/4=$(round(quantile(rHn,1/4),4)) q1/2=$(round(quantile(rHn,1/2),4)) q3/4=$(round(quantile(rHn,3/4),4))")
end


# for meshcase in [3,1,2,5,6]
#     println("meshcase=$meshcase")
#     include("fielddiary.jl")
#     println("N = $(size(points,2))")
# end

