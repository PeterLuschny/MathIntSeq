# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

#=
Executing SeqBuild.jl produces MathIntSeq.jl which is the
  main module. It also produces
  docs/src/index and
  docs/src/modules which are the starting points for the docs.

! While building watch for duplicates (i.e. implementations in different
! modules of the same sequence).
=#

firsttime = false
if firsttime
    Pkg.add("Nemo")
    Pkg.add("Memoize")
    Pkg.add("OffsetArrays")
    Pkg.add("HTTP")
    Pkg.add("IterTools")
#   Pkg.add("Documenter")
#   Pkg.add("Lint")
    Pkg.update()
end

ROOTDIR = Pkg.dir("MathIntSeq")
srcdir = joinpath(ROOTDIR, "src")
tstdir = joinpath(ROOTDIR, "test")

if srcdir ∉ LOAD_PATH
    push!(LOAD_PATH, srcdir)
end
info(LOAD_PATH)

using MathIntSeqBuild #, Lint

cd(srcdir)

info(" *** building MathIntSeq")
run(`julia MathIntSeqBuild.jl`)

#info(" *** running MathIntSeq")
#run(`julia MathIntSeq.jl`)

#info(" *** lint MathIntSeq running ...")
#lintfile("MathIntSeq.jl")

cd(tstdir)

#info(" *** lint runtests running ...")
#lintfile("runtests.jl")

info(" *** runtests")
run(`julia runtests.jl`)

#info(" *** documenting performance")
#run(`julia perftests.jl`)
#info(" results are in performance.txt")

cd(ROOTDIR)

info(" *** build finished")

#info("  -  TODO: check for performance regressions")
#info("  -  TODO: update git")
#info("  -  TODO: docs/makedocs.jl")
#info("  -  TODO: publish docs/site")
#info("  -  TODO: deploy MathIntSeq.jl")
