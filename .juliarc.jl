BaseDir = dirname(@__FILE__)
#push!(LOAD_PATH,BaseDir * "/modules")

ENV["JULIA_PKGDIR"] = BaseDir * "/packages"
ENV["LOAD_CACHE_PATH"] = homedir() * "/.julia/.cache/"

println("\n")

