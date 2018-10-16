using Pkg
pkg"activate ."
pkg"resolve"

BaseDir = dirname(@__FILE__)

ENV["JULIA_PKGDIR"] = BaseDir * "/packages"
ENV["LOAD_CACHE_PATH"] = homedir() * "/.julia/.cache/"
#push!(LOAD_PATH,BaseDir * "/libs")

datadir = homedir()*"/SimulationData/"
if !isdir(datadir)
    mkdir(datadir)
end

println("\n")

