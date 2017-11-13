# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module MdTrans

function md_trans()
    docsdir = realpath(joinpath(dirname(@__FILE__)))
    builddir = joinpath(docsdir, "build")
    path = joinpath(builddir, "index.md")
    npath = joinpath(builddir, "nindex.md")

    src = open(path, "r")
    lines = readlines(src)
    close(src)

    sline = line = 0
    t = ["", ""]
    for l in lines
        line += 1
        n = lstrip(l)
        n == "" && continue
        if startswith(n, "<a id")
            t = split(n, '#')
            sline = line
            lines[line + 1] = ""
            continue
        end
        if startswith(n, "<a target")
            s = split(n, '#')
            L = split(s[2], ' ')
            M = split(t[2], '\'')
            lines[sline] = t[1] * "#" * L[1] * ">" * M[1] * "</a>"
            lines[line] = ""
            continue
        end
        if startswith(n, "* function")
            lines[line] = "  *" * n[11:end]
            continue
        end
        if startswith(n, "* ")
            lines[line] = "  * " * n[3:end]
            continue
        end
        lines[line] = "** " * n
    end

    dest = open(npath, "w")
    for l in lines
        print(dest, l)
    end
    close(dest)

    cp(npath, path, remove_destination = true)
    rm(npath)
end

md_trans()

end #  module
