# add to ~/.juliarc.jl
# push!(LOAD_PATH, "/pathTo/MathIntSeq/src")
# push!(LOAD_PATH, "/pathTo/MathIntSeq/src/modules")

devdir = realpath(joinpath(dirname(@__FILE__)))
srcdir = dirname(devdir)
moddir = joinpath(srcdir, "modules")

push!(LOAD_PATH, srcdir)
push!(LOAD_PATH, moddir)

info(LOAD_PATH)
