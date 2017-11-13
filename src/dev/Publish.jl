# https://docs.julialang.org/en/release-0.6/manual/packages/#Generating-the-package-1
# https://github.com/JuliaLang/PkgDev.jl/blob/master/README.md
# https://docs.julialang.org/en/stable/manual/packages/

# Pkg.add("PkgDev")
using PkgDev

PkgDev.generate("MathIntSeq", "MIT")

# PkgDev.config()
    # PkgDev.jl configuration:
    # User name: OpenLibMathSeq
    # User email: seq@olms.onl
    # Enter GitHub user [OpenLibMathSeq]:

Pkg.update()
Pkg.resolve()

###
# https://github.com/attobot/attobot
###

# Pkg.status()
#PkgDev.register("MathIntSeq")
#PkgDev.tag("MathIntSeq", v"0.0.1")
#PkgDev.publish()
