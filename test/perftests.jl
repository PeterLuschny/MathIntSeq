# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.
# This file includes portions from Combinatorics.jl and Primes.jl
# in modified form. License is MIT, http://julialang.org/license

# Version of: UTC 2018-01-31 13:17:47
# bfc8ab60-0680-11e8-046d-eb5cdf504ac0

# Do not edit this file, it is generated from the modules and will be overwritten!
# Edit the modules in the modules directory and build this file with MathIntSeqBuild.jl!

# Reporting a bug please give date, uuid and line number in this file, not of modules.

module perftests
using Nemo, MathIntSeq
versioninfo()
start = now()
# +++ Abundant.jl +++
L005101(10)
println("\nTEST: L005101(80000)")
@time L005101(80000)
gc()
L002093(10)
println("\nTEST: L002093(200)")
@time L002093(200)
# +++ Andre.jl +++
gc()
println("\nTEST: (for m in 1:20, n in 0:100 André(m, n) end)")
@time (for m in 1:20, n in 0:100 André(m, n) end)
# +++ BellNumbers.jl +++
gc()
BellNumberList(5)
println("\nTEST: BellNumberList(1000)")
@time BellNumberList(1000)
println("\nTEST: BellTriangle(100, [1 for _ in 1:100])")
@time BellTriangle(100, [1 for _ in 1:100])
# +++ BernoulliNumbers.jl +++
gc()
for m in 1:2, n in 0:1 BernoulliInt(m, n) end
println("\nTEST: (for m in 1:10, n in 0:499 BernoulliInt(m, n) end)")
@time (for m in 1:10, n in 0:499 BernoulliInt(m, n) end)
for m in 1:2, len in 1:2 BernoulliIntList(m, len) end
println("\nTEST: (for m in 1:10 BernoulliIntList(m, 500) end)")
@time (for m in 1:10 BernoulliIntList(m, 500) end)
gc()
for n in 0:10 Bernoulli(n) end
println("\nTEST: (for n in 0:1000 Bernoulli(n) end)")
@time (for n in 0:1000 Bernoulli(n) end)
BernoulliList(10)
println("\nTEST: BernoulliList(1000)")
@time BernoulliList(1000)
for n in 0:10 A027641(n) end
println("\nTEST: (for n in 0:1000 A027641(n) end)")
@time (for n in 0:1000 A027641(n) end)
gc()
A195441(10)
println("\nTEST: (for n in 0:10000 A195441(n) end)")
@time (for n in 0:10000 A195441(n) end)
# +++ BinaryInteger.jl +++
println("\nTEST: [A295513(k) for k in 0:100000]")
@time [A295513(k) for k in 0:100000]
# +++ BinaryQF.jl +++
gc()
println("\nTEST: B034017(10000)")
@time B034017(10000)
# +++ Clausen.jl +++
gc()
println("\nTEST: ClausenNumberList(10000)")
@time ClausenNumberList(10000)
gc()
println("\nTEST: (for n in 0:10000 ClausenNumber(n) end)")
@time (for n in 0:10000 ClausenNumber(n) end)
# +++ CombinationsIterator.jl +++
# +++ Counts.jl +++
println("\nTEST: [A257993(n) for n in 1:10000]")
@time [A257993(n) for n in 1:10000]
gc()
println("\nTEST: PrimePiList(10000)")
@time PrimePiList(10000)
# +++ DedekindEta.jl +++
println("\nTEST: PartitionNumberList(10000)")
@time PartitionNumberList(10000)
println("\nTEST: L000731(10000)")
@time L000731(10000)
println("\nTEST: RamanujanTauList(10000)")
@time RamanujanTauList(10000)
# +++ Deleham.jl +++
gc()
println("\nTEST: T225478(100)")
@time T225478(100)
println("\nTEST: T055883(100)")
@time T055883(100)
# +++ Derangements.jl +++
gc()
println("\nTEST: (for n in 0:450 Subfactorial(n) end)")
@time (for n in 0:450 Subfactorial(n) end)
# +++ Fibonacci.jl +++
gc()
println("\nTEST: (for n in 1:10000 FibonacciGeneralized(n, 1) end)")
@time (for n in 1:10000 FibonacciGeneralized(n, 1) end)
println("\nTEST: FibonacciGeneralizedList(10000, 1)")
@time FibonacciGeneralizedList(10000, 1)
# +++ FigurativeNumbers.jl +++
# +++ GaussFactorial.jl +++
gc()
A001783(10)
println("\nTEST: L001783(1000)")
@time L001783(1000)
gc()
A055634(10)
println("\nTEST: [A055634(n) for n in 1:1000]")
@time [A055634(n) for n in 1:1000]
gc()
A103131(10)
println("\nTEST: [A103131(n) for n in 1:1000]")
@time [A103131(n) for n in 1:1000]
# +++ GeneralBinomial.jl +++
println("\nTEST: (for n in 0:10000 Binomial(2n, n) end)")
@time (for n in 0:10000 Binomial(2n, n) end)
println("\nTEST: (for n in -100:100, k in -100:100 Binomial(n, k) end)")
@time (for n in -100:100, k in -100:100 Binomial(n, k) end)
println("\nTEST: (for k in -10000:10000 Binomial(-5, k) end)")
@time (for k in -10000:10000 Binomial(-5, k) end)
# +++ Hyper1F1.jl +++
gc()
println("\nTEST: (for n in 0:150 A000262(n) end)")
@time (for n in 0:150 A000262(n) end)
# +++ JacobiTheta.jl +++
println("\nTEST: L000143(100000)")
@time L000143(100000)
println("\nTEST: L035016(100000)")
@time L035016(100000)
# +++ Kolakoski.jl +++
println("\nTEST: KolakoskiList(10000)")
@time KolakoskiList(10000)
# +++ Maxima.jl +++
# +++ NumberTheory.jl +++
println("\nTEST: [Divisors(n) for n in 1:10000]")
@time [Divisors(n) for n in 1:10000]
println("\nTEST: [Radical(n)  for n in 1:10000]")
@time [Radical(n)  for n in 1:10000]
# +++ OEISUtils.jl +++
# +++ OrthoPolynomials.jl +++
T111062(2)
T066325(2)
T053120(2)
gc()
println("\nTEST: T111062(1000)")
@time T111062(1000)
gc()
println("\nTEST: T066325(1000)")
@time T066325(1000)
gc()
println("\nTEST: T053120(1000)")
@time T053120(1000)
# +++ Partitions.jl +++
# +++ PrimeSieve.jl +++
# +++ Products.jl +++
gc()
println("\nTEST: (for n in 1:10000 F!(n) end)")
@time (for n in 1:10000 F!(n) end)
gc()
println("\nTEST: (for n in 1:1000 A000407(n) end)")
@time (for n in 1:1000 A000407(n) end)
gc()
println("\nTEST: (for n in 1:200, k in 1:200 RisingFactorial(n, k) end)")
@time (for n in 1:200, k in 1:200 RisingFactorial(n, k) end)
# +++ Recursive2.jl +++
gc()
println("\nTEST: L286032(10000)")
@time L286032(10000)
# +++ SelfConvolutive.jl +++
gc()
L005411(10)
println("\nTEST: L005411(500)")
@time L005411(500)
# +++ SeqBase.jl +++
# +++ StirlingLahNumbers.jl +++
gc()
println("\nTEST: StirlingSetTriangle(1000, 2)")
@time StirlingSetTriangle(1000, 2)
gc()
println("\nTEST: LahTriangle(1000, 2)")
@time LahTriangle(1000, 2)
# +++ SwingingFactorial.jl +++
gc()
println("\nTEST: (for n in 0:999 Swing(n) end)")
@time (for n in 0:999 Swing(n) end)
println("\nTEST: (for n in 0:999 A000984(n) end)")
@time (for n in 0:999 A000984(n) end)
# +++ ZumkellerNumbers.jl +++
IsZumkeller(10)
gc()
println("\nTEST: (for n in 1:2000 IsZumkeller(n) end)")
@time (for n in 1:2000 IsZumkeller(n) end)
stop = now()
tdiff = Int(stop - start) / 1000
println("\nJulia version: " * string(VERSION) )
println(start)
println("Total test time: ", tdiff, " seconds")
end # module
