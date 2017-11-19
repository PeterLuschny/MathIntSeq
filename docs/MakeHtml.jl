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

function compressfile()

    docdir = realpath(joinpath(dirname(@__FILE__)))
    head = open(joinpath(docdir, "head.html"), "r")
    art  = open(joinpath(docdir, "article.html"), "r")
    tail = open(joinpath(docdir, "tail.html"), "r")
    dest = open(joinpath(docdir, "MathIntSeq.html"), "w")

    for l in eachline(head) print(dest, l) end
    for l in eachline(art)  print(dest, l) end
    for l in eachline(tail) print(dest, l) end

    head = close(head)
    art  = close(art)
    tail = close(tail)
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

function writearticle()

    docdir = realpath(joinpath(dirname(@__FILE__)))
    pkgdir = dirname(docdir)
    srcdir = joinpath(pkgdir, "src")
    mis = open(joinpath(srcdir, "MathIntSeq.jl"), "r")
    doc = open(joinpath(docdir, "article.html"), "w")

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

function makedoc()
    MathIntSeqBuild.build_all(true)
    writearticle()
    compressfile()
    MathIntSeqBuild.build_all(false)

    docdir = realpath(joinpath(dirname(@__FILE__)))
    src = joinpath(docdir, "MathIntSeq.html")
    dst = joinpath(joinpath(docdir, "site"), "index.html")
    cp(src, dst; remove_destination=true)
end

makedoc()

end # module
