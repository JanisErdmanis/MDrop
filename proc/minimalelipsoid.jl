### Implementation of mathematica thing

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

function TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Hx = EllipsoidField(a,b,c,mup,[H0,0,0])
    Hy = EllipsoidField(a,b,c,mup,[0,H0,0])

    Emt = 1/8/pi*(1 - mup) * (dot(Hx,H0*[1,0,0]) + dot(Hy,H0*[0,1,0]))/2 * 4/3*pi * a * b * c
    return Emt
end

import Elliptic
function EllipsoidArea(a,b,c)

    if a<b
        a,b = b,a
    end
    if b<c
        b,c = c,b
    end
    if a<b
        a,b = b,a
    end

    cosphi = c/a
    phi = acos(cosphi)
    sinphi = sqrt(1 - cosphi^2)
    k = sqrt(a^2*(b^2-c^2)/b^2/(a^2-c^2))

    S = 2*pi*c^2 + 2*pi*a*b/sinphi * (sinphi^2 * Elliptic.E(phi,k) + cosphi^2 * Elliptic.F(phi,k))

    return S
end

function TheoreticalDropEnergy(a,b,c,mup,Bm)
    ### I should scale a,b,c inside
    a,b,c = abs(a),abs(b),abs(c)
    # if a<0 || b<0 || c<0
    #     return a^2 + b^2 + c^2
    # end
        
    #println("a=$a b=$b c=$c")
    gammap = 1.
    H0 = sqrt(Bm*gammap/(a*b*c)^(1/3))
    Emt = TheoreticalRotatingFieldEnergy(a,b,c,mup,H0)

    Etotal = Emt + gammap*EllipsoidArea(a,b,c) # + surface area
end



using Optim

mu = 4.8

ab = Float64[]
E = Float64[]
Bm = Float64[]

ai = 1.
bi = 1.1
Bmi = 40.

while Bmi>1.

    f(x) = TheoreticalDropEnergy(x[1],x[2],1./x[1]/x[2],mu,Bmi)

    lower = [1., 1.]
    upper = [Inf, Inf]
    initial_x = [ai, bi]
    #results = optimize(DifferentiableFunction(f), initial_x, lower, upper, Fminbox(), optimizer = GradientDescent)
    results = optimize(f,[ai,bi])
    #results = optimize(f,[1.1,2.])

    ai,bi = Optim.minimizer(results)
    Ei = Optim.minimum(results)

    Einorm = Ei/4/pi

    push!(ab,ai/bi)
    push!(E,Einorm)
    push!(Bm,Bmi)

    Bmi -= 1.

    println("Bmi=$Bmi")    

end

using DataFrames

df = DataFrame(Bm=Bm,E=E,ab=ab)
writetable("$(dirname(@__FILE__))/raimonds/myatempt.dat",df,separator='\t')
