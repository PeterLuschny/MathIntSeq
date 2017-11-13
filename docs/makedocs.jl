# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

using Documenter

#=
https://juliadocs.github.io/Documenter.jl/latest/man/guide.html

(1) Executing makedocs() and gives
    MathIntSeq/docs/build/index.md.
(2) In MathIntSeq/docs running "mkdocs build" gives
    MathIntSeq/docs/site/index.html.
(3) Transform MathIntSeq/docs/site/index.html with HtmlTrans.jl
    to a file with the same name in the same dir.
=#

docdir = realpath(joinpath(dirname(@__FILE__)))
cd(docdir)

info(" *** running makedocs")
Documenter.makedocs(root = docdir)

info(" *** running MdTrans")
run(`julia MdTrans.jl`)

info(" *** running mkdocs")
run(`mkdocs build`)

info(" *** running HtmlTrans")
run(`julia HtmlTrans.jl`)

# Documenter: copying assets to build directory.
# !! Overwriting 'build\assets\mathjaxhelper.js'.
# However we want our own mathjax setup!

file = "mathjaxhelper.js"
src = joinpath(joinpath(joinpath(docdir, "src"), "assets"), file)
dst = joinpath(joinpath(joinpath(docdir, "site"), "assets"), file)
cp(src, dst; remove_destination = true)

# We don't want this big html file on github.
file = "SeqNotebook.html"
src = joinpath(joinpath(dirname(docdir), "local"), file)
dst = joinpath(joinpath(joinpath(docdir, "site"), "notebook"), file)
cp(src, dst; remove_destination = true)

info(" *** docs done")

function open_file(filename)
    if is_apple()
        run(`open $(filename)`)
    elseif is_linux() || is_bsd()
        run(`xdg-open $(filename)`)
    elseif is_windows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    else
        warn("This file is not supported on OS $(string(Compat.KERNEL))")
    end
end

info(" *** running mkdocs as server")
run(`mkdocs serve`)
open_file("http://127.0.0.1:8000/")
