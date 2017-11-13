# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module HtmlTrans

bas = "<base href=\"https://github.com/OpenLibMathSeq/MathIntSeq/blob/master/src/MathIntSeq.jl\" target=\"_blank\" />"
loc = "http://olms.onl/julia/MathIntSeq/"

footer = ["</footer>", "</div>",
"<script src=\"http://olms.onl/julia/IntSeq/assets/javascripts/application-7aa26ad9ec.js\"></script>",
"<script type=\"text/x-mathjax-config\">",
"MathJax.Hub.Config({",
"    \"HTML-CSS\": {",
"     availableFonts: [],",
"     preferredFont: null,",
"     webFont: \"Neo-Euler\",",
"     scale: 90",
"    },",
"    displayAlign: \"left\",",
"    displayIndent: \"1.5em\",",
"    TeX: { equationNumbers: { autoNumber: \"AMS\" } },",
"    tex2jax: {inlineMath: [[\"*\",\"*\"]]}",
"});",
"</script>",
"<script type=\"text/javascript\" async",
"src=\"https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">",
"</script>",
"</body>",
"</html>"]

function trans()
    docsdir = realpath(joinpath(dirname(@__FILE__)))
    sitedir = joinpath(docsdir, "site")
    path = joinpath(sitedir, "index.html")
    npath = joinpath(sitedir, "nindex.html")

    src = open(path, "r")
    dest = open(npath, "w")

    sty = "<style>pre.a{margin-left:10px;font-size:10pt}p.c{color:#3F51B5;}p.b{margin-left:14px;color:brown;}"

    for l in eachline(src)
        n = lstrip(l)
        g = length(n)
        n == "" && continue
        startswith(n, "</footer>") && break

        if startswith(n, "<style>")
            print(dest, sty * n[8:end])
            continue
        end

        n = replace(n, "<p>*<em> ", "<p>**")
        n = replace(n, "</head>", bas * "</head>")
        n = replace(n, "<p><a id", "<p class=\"c\"><a id")
        n = replace(n, "<li>", "<li><code>")
        n = replace(n, "</li>", "</code></li>")
        n = replace(n, "<p>**", "<p class='b'>")
        n = replace(n, "</code>**", "<p class='b'>")
        n = replace(n, ">MathIntSeq.", ">")
        ##
        n = replace(n, "about", "./about")
        n = replace(n, "license", "./license")
        n = replace(n, "modules", "./modules")
        n = replace(n, "repository", "./repository")
        n = replace(n, "href=\".\"", "href=\"./\"")
        n = replace(n, "notebook", "./notebook")
        ##
        n = replace(n, "./", loc)
        if startswith(n, "<p class='b'>") && contains(n, "*")
            z = split(n, '*')
            println(n)
            println(z)
            a = z[1] * "</p>"
            b = "<code>" * z[2]
            b = replace("</p>", "</code>")
        end
        if startswith(n, "<li><code>")
            n = replace(n, "<em>", "×")
            n = replace(n, "</em>", "×")
        end
        print(dest, n)
    end
    for l in footer
        n = replace(l, '*', '$')
        println(dest, n)
    end
    close(src)
    close(dest)
    cp(npath, path, remove_destination = true)
    rm(npath)
end

function shrink(name)
    docsdir = realpath(joinpath(dirname(@__FILE__)))
    sitedir = joinpath(docsdir, "site")
    destdir = joinpath(sitedir, name)
    path = joinpath(destdir, "index.html")
    npath = joinpath(destdir, "nindex.html")

    src = open(path, "r")
    dest = open(npath, "w")

    for l in eachline(src)
        n = lstrip(l)
        n == "" && continue
        n = replace(n, "../", loc)
        print(dest, n)
    end
    close(src)
    close(dest)
    cp(npath, path, remove_destination = true)
    rm(npath)
end

trans()
shrink("about")
shrink("license")
shrink("modules")
shrink("notebook")

end #  module
