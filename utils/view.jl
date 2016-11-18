#!/usr/bin/julia

using SurfaceGeometry

using JLD
SimulationData = homedir() * "/SimulationData/"

using Escher
include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
using ThreeJS
using Compat

info("Modules are loaded")

@spawn run(`chromium-browser --app=http://0.0.0.0:5555/NeatFrameViewer.jl`)

viewers = dirname(@__FILE__) *"/Viewers"
@spawn escher_serve(5555,viewers)

import Base.isless
function Base.isless(x::ASCIIString,y::ASCIIString)
    nx = parse(Int,x[1:length(x)-4])
    ny = parse(Int,y[1:length(y)-4])
    return nx < ny
end

memory = nothing

while true
    try
        println("Which integrator data do you want to check")

        content = readdir(SimulationData)
        for i in 1:length(content)
            println("$i \t $(content[i])")
        end

        N = parse(Int,readline(STDIN))

        directory = content[N]
        
        content = readdir(SimulationData*"/"*directory)
        for i in 1:length(content)
            println("\t$i \t $(content[i])")
        end

        N = parse(Int,readline(STDIN))

        if isfile(SimulationData*"/"*directory*"/"*content[N])
            memory = load(SimulationData*"/"*directory*"/"*content[N],"memory")
        else
            info("Loading a set of files")
            outfiles = readdir(SimulationData*"/"*directory*"/"*content[N])
            sort!(outfiles)

            memory = []
            for i in 1:length(outfiles)
                memi = load(SimulationData*"/"*directory*"/"*content[N]*"/"*outfiles[i],"memory")
                memory = [memory;memi]
            end
        end
        println("File $(directory*"/"*content[N]) loaded. It has $(length(memory)) frames.")

        ### Here I now also can extract meaningfull parameters from configuration file
        
    catch
        info("Error")
    end
end
