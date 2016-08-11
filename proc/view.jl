#!/usr/bin/julia

ENV["JULIA_PKGDIR"] = dirname(@__FILE__) * "/Packages"
using JLD
#using Storage
using SurfaceGeometry


### A mesh viewers

calculate = false
include("manager.jl")
#include("diarymag.jl")

using Escher
include(Pkg.dir("Escher", "src", "cli", "serve.jl"))
using ThreeJS
using Compat

info("Modules are loaded")

if isinteractive()
    isdefined(:viewmeshactive) || (viewmeshactive=false)
    if viewmeshactive==false
        @spawn escher_serve(5555,"Viewers")
        viewmeshactive = true
    end
else
    escher_serve(5555,"Viewers")
end

