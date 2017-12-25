srcdir = dirname(realpath(joinpath(dirname(@__FILE__))))

if srcdir ∉ LOAD_PATH
    push!(LOAD_PATH, srcdir)
    push!(LOAD_PATH, joinpath(srcdir, "modules"))
    info(LOAD_PATH)
end

module CompileCache

using SeqBase, OEISUtils, SeqTests, PrimeSieve, Products, NumberTheory
using GeneralBinomial, Deleham, Abundant, StirlingLahNumbers
using Andre, Clausen, BernoulliNumbers, ZumkellerNumbers, BinaryQF, DedekindEta
using Fibonacci, GaussFactorial, Hyper1F1, JacobiTheta, Kolakoski, Maxima
using Partitions, Recursive2, BellNumbers, OrthoPolynomials
using SwingingFactorial, SelfConvolutive, FigurativeNumbers, BinaryInteger

end # module
