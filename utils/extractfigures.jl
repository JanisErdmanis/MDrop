using JLD

import Base.isless
function Base.isless(x::ASCIIString,y::ASCIIString)
    nx = parse(Int,x[1:length(x)-4])
    ny = parse(Int,y[1:length(y)-4])
    return nx < ny
end

#simulation = "sphere0.1:mu=4.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=6.0;Bm=200.0;omega=0.0"
#simulation = "sphere0.1:mu=10.0;Bm=170.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=172.0;omega=0.0"
#simulation = "Psphere0.2:mu=10.0;Bm=158.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=50.0;omega=0.0"
#simulation = "sphere0.2:mu=100.0;Bm=50.0;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:36;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:81;omega=0.0"
#simulation = "sphere0.2:mu=10.0;Bm=1:84;omega=0.0"
simulation = "sphere0.1:mu=10.0;Bm=1:84;omega=0.0"

input = "$datadir/QstaticFastFieldEiler/$simulation/"
output = "data/outdir/QstaticFastFieldEiler:$simulation"

# simulation = "mu10Bm50:mu=10.0;Bm=50.0;omega=0.0"
# input = "$datadir/QstaticFastFieldEilerDecreasing/$simulation/"
# output = "$datadir/ExtractedEqFigures/QstaticFastFieldEilerDecreasing:$simulation"

if !isdir(dirname(output))
    mkdir(dirname(output))
end

if isdir(output)
    run(`rm -rf $output`)
end
mkdir(output)

infiles = readdir(input)
sort!(infiles)

Bmi = 1
#Bmi = 50
for i in 1:(length(infiles) - 1)
    memi = load("$input/$(infiles[i])","memory")
    if length(memi[end])>=4
        Ei = memi[end][4]
    else
        Ei = nothing
    end

    if length(memi[end])>=5
        taui = memi[end][5]
    else
        taui = nothing
    end
    
    save("$output/$Bmi.jld","points",memi[end][2],"faces",memi[end][3],"Ei",Ei,"taui",taui)
    Bmi += 1
    #Bmi += -1
end
