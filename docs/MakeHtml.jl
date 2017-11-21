# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

docdir = realpath(joinpath(dirname(@__FILE__)))
srcdir = joinpath(dirname(docdir), "src")
if srcdir ∉ LOAD_PATH
    push!(LOAD_PATH, srcdir)
end
info(LOAD_PATH)

module HtmlDoc
using MathIntSeqBuild

#gitdir = "https://github.com/OpenLibMathSeq/MathIntSeq/blob/master/src/modules/"
gitdir = "https://github.com/PeterLuschny/MathIntSeq/blob/master/src/modules/"

function compressfile(_head, _art, _foot, _dest)

    docdir = realpath(joinpath(dirname(@__FILE__)))
    head = open(joinpath(docdir, _head), "r")
    art  = open(joinpath(docdir, _art), "r")
    foot = open(joinpath(docdir, _foot), "r")
    dest = open(joinpath(docdir, _dest), "w")

    for l in eachline(head) print(dest, l) end
    for l in eachline(art)  print(dest, l) end
    for l in eachline(foot) print(dest, l) end

    head = close(head)
    art  = close(art)
    tail = close(foot)
    dest = close(dest)
end

function escape_chars(i::AbstractString)
    o = replace(i, "&", "&amp;")
    o = replace(o, "\"", "&quot;")
    o = replace(o, "'", "&#39;")
    o = replace(o, "<", "&lt;")
    o = replace(o, ">", "&gt;")
    o = replace(o, "``", "\$")
    return o
end

function writeseqarticle()

    docdir = realpath(joinpath(dirname(@__FILE__)))
    pkgdir = dirname(docdir)
    srcdir = joinpath(pkgdir, "src")
    mis = open(joinpath(srcdir, "MathIntSeq.jl"), "r")
    doc = open(joinpath(docdir, "seqarticle.html"), "w")

    docline = false
    example = false
    line = 0
    linenumber = 0

    buffer0 = IOBuffer()
    buffer1 = IOBuffer()
    def = ""
    tit = ""

    for l in eachline(mis)
        line += 1
        n = lstrip(l)
        g = length(n)
        n == "" && continue

        if startswith(n, "doc\"\"\"")
            print(buffer1, "<p class=\"b\">")
            docline = true
            example = false
            linenumber = line
            line -= 1
            continue
        elseif startswith(n, "\"\"\"")
            println(doc, tit)
            println(buffer1, "</p>")
            print(doc, String(take!(buffer1)))
            println(doc, def)
            if example
                println(buffer0, "</code></p>")
                println(doc, String(take!(buffer0)))
            end
            docline = false
            continue
        elseif startswith(n, "julia")
            example = true
            println(buffer0, "<p style=\"margin-left: 30px;\"><code>")
        end

        if docline
            s = escape_chars(n)
            if startswith(s, "*")
                s = replace(s, "* ", "")
                s = replace(s, "function ", "")
                i = searchindex(s, "(")
                name = s[1:i - 1]
                tit = "<p class=\"c\"><a id=\"$(name * string(linenumber))\" href=\"#L$(linenumber)\">$(name)</a></p>"
                def = "<ul class=\"def\"><li>" * s * "</li></ul>"
            else
                if example println(buffer0, s, "<br>") else print(buffer1, " ", s) end
            end
        end
    end

    close(buffer0)
    close(buffer1)
    close(mis)
    close(doc)
end

function writemodarticle()

    docdir = realpath(joinpath(dirname(@__FILE__)))
    pkgdir = dirname(docdir)
    srcdir = joinpath(pkgdir, "src")
    moddir = joinpath(srcdir, "modules")
    modart = open(joinpath(docdir, "modarticle.html"), "w")

    exclude = ["OLMS.jl", "SeqTests.jl"]
    seq_modules = filter!(r"\.jl$", readdir(moddir))
    for filename in seq_modules
        filename in exclude && continue
        name = split(filename, ".")
        println(modart, "<a href=\"" * filename * "\">" * name[1] * "</a>  ")
    end
    close(modart)
end

function makedoc()
    MathIntSeqBuild.build_all(true)
    writeseqarticle()
    compressfile("seqhead.html", "seqarticle.html", "seqfoot.html", "MathIntSeq.html")
    MathIntSeqBuild.build_all(false)

    writemodarticle()
    compressfile("modhead.html", "modarticle.html", "modfoot.html", "modules.html")

    docdir = realpath(joinpath(dirname(@__FILE__)))
    src = joinpath(docdir, "MathIntSeq.html")
    dst = joinpath(joinpath(docdir, "site"), "index.html")
    cp(src, dst; remove_destination=true)

    src = joinpath(docdir, "modules.html")
    dst = joinpath(joinpath(joinpath(docdir, "site"), "modules"), "index.html")
    cp(src, dst; remove_destination=true)
end

makedoc()

end # module
