ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"

### Calculation of big axis

function mainaxis(points)

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

    return nl1
end

function angledep(memory)
    tarr = []
    
    phiarr = []
    prevphi = 0

    jump = 0
    
    for tkey in 1:10:length(memory)
        nl = mainaxis(memory[tkey][2])

        x = dot(nl,[1,0,0])
        y = dot(nl,[0,1,0])
        #phi = acos(dot(nl,[1,0,0]))
        #phi = atan2(y,x)
        phi = atan(y/x)
        #prevphi = phiarr[end]

        if (prevphi-phi)>pi*0.7
            jump += pi
        end
        
        push!(tarr,memory[tkey][1])
        push!(phiarr,phi+jump)
        #println("t=$(memory[tkey][1]) phi=$phi")
        prevphi = phi
    end

return tarr,phiarr
end

function rotatingspeed(tarr)


    tarr,phiarr = angledep(memory)
    len = length(tarr)
    
    dt = tarr[len] - tarr[len-window]
    dphi = phiarr[len] - phiarr[len-window]

    return dphi/dt
end

#case = 17

# #calcpar = "calcpar/Bm35omega0.175.jl"
# calcpar = "calcpar/Bm35omega0.35.jl"

parameters = ["calcpar/dynamics/Bm35omega0.175.jl","calcpar/dynamics/Bm35omega0.20.jl","calcpar/dynamics/Bm35omega0.25.jl","calcpar/dynamics/Bm35omega0.35.jl","calcpar/dynamics/Bm35omega0.4.jl","calcpar/dynamics/Bm35omega0.3.jl","calcpar/dynamics/Bm35omega0.325.jl","calcpar/dynamics/Bm35omega0.275.jl"] # ,"calcpar/dynamics/Bm35omega0.5.jl",
windows = [40,40,40,40,40,25,28,40]

# parameters = ["calcpar/dynamics/Bm35omega0.325.jl",]
# windows = [30,]

omegaField = []
omegaDrop = []

using PyPlot; const plt = PyPlot
fig = plt.figure()

for (calcpar,window) in zip(parameters,windows)

    println("$calcpar")

    calculate = false
    include("manager.jl");

    #println("angular frequency is $(rotatingspeed(memory,20))")

    tarr,phiarr = angledep(memory)
    #plt.plot(tarr,phiarr)
    rrr = 1:length(tarr)
    plt.plot(rrr,phiarr,label="$calcpar")
    plt.plot(rrr[length(tarr) - window],phiarr[length(tarr) - window],"ko")

    len = length(tarr)
    dt = tarr[len] - tarr[len-window]
    dphi = phiarr[len] - phiarr[len-window]

    tau = etap*(volume(memory[end][2],memory[end][3])/ 4/pi*3)^(1/3)/gammap
        
    push!(omegaDrop,dphi/dt*tau)
    push!(omegaField,omega*tau)
end

plt.legend()
fig[:savefig]("figures/rotatingfreqcalc.pdf")
plt.close("all")

plt.rc("text",usetex=true)
plt.rc("font", family="serif", serif="Times")
plt.rc("xtick",labelsize=10)
plt.rc("ytick",labelsize=10)
plt.rc("axes",labelsize=10)
plt.rc("legend",fontsize=8)

width = 3.487 
height = width / 1.618
fig = plt.figure(figsize=(width,height))

plt.plot(omegaField,omegaDrop./omegaField,".")

plt.title(L"Bm=35;\mu=10")
plt.xlabel(L"\omega_{field} \cdot \gamma/\eta R_0")
plt.ylabel(L"\omega_{drop}/\omega_{field}")
plt.ylim(0,1.1)
plt.xlim(0,0.5)

plt.grid("off")
plt.tight_layout()
fig[:savefig]("figures/rotatingfrq.pdf")
plt.close("all")


#isinteractive()
