BaseDir = dirname(@__FILE__)
push!(LOAD_PATH,BaseDir * "/modules")

ENV["JULIA_PKGDIR"] = BaseDir * "/packages"

