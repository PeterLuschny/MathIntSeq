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
    Pkg.add("Combinatorics")
    Pkg.add("OffsetArrays")
    Pkg.add("Requests")
    Pkg.add("URIParser")
    Pkg.add("Documenter")
    Pkg.add("Lint")
    Pkg.update()
end
using Lint

srcdir = realpath(joinpath(dirname(@__FILE__)))
pkgdir = dirname(srcdir)
tstdir = joinpath(pkgdir, "test")

if srcdir ∉ LOAD_PATH
    push!(LOAD_PATH, srcdir)
end

info(LOAD_PATH)

cd(srcdir)

info(" *** building MathIntSeq")
run(`julia MathIntSeqBuild.jl`)

info(" *** running MathIntSeq")
run(`julia MathIntSeq.jl`)

#info(" *** lint MathIntSeq running ...")
#lintfile("MathIntSeq.jl")

pkg = joinpath(Pkg.dir(), "MathIntSeq")
src = joinpath(pkg, "src")
mis = joinpath(src, "MathIntSeq.jl")
cp("MathIntSeq.jl", mis; remove_destination=true)

cd(tstdir)

#info(" *** lint runtests running ...")
#lintfile("runtests.jl")

info(" *** runtests")
run(`julia runtests.jl`)

info(" *** documenting performance")
run(`julia perftests.jl`)
# info(" results are in performance.txt")

cd(srcdir)

info(" *** build finished successfully")

#info("  -  TODO: check for performance regressions")
#info("  -  TODO: update git")
#info("  -  TODO: docs/makedocs.jl")
#info("  -  TODO: publish docs/site")
#info("  -  TODO: deploy MathIntSeq.jl")
