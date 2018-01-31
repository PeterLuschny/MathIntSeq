# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module BernoulliNumbers
using Nemo, NumberTheory, PrimeSieve, SeqBase, Andre, Products

export BernoulliInt, BernoulliIntList, Bernoulli, BernoulliList
export A195441, A065619, A281586, A281588, A027641, L065619

# The rational Bernoulli numbers are defined here with B_(1) = 1/2.
# Why this is preferred over B_(1) = -1/2 is explained in
# http://luschny.de/math/zeta/The-Bernoulli-Manifesto.html

doc"""
Return the generalized integer Bernoulli numbers ``b_{m}(n) = n`` André``_{m}(n-1)``.

julia> BernoulliInt(3, 10)
-15130
"""
BernoulliInt(m::Int, n::Int) = n == 0 ? ZZ(0) : n * André(m, n - 1)

doc"""
Return the number of down-up permutations w on ``[n+1]`` such that
``w_2 = 1``. (D. Callan)

julia> [A065619(n) for n in 0:9]
[0, 1, 2, -3, -8, 25, 96, -427, -2176, 12465]
"""
A065619(n::Int) = BernoulliInt(2, n)

doc"""
Return the generalized integer Bernoulli numbers ``b_{3}(n) = n``
André``_{3}(n-1)``.

julia> [A281586(n) for n in 0:9]
[0, 1, 2, 3, -4, -15, -54, 133, 792, 4293]
"""
A281586(n::Int) = BernoulliInt(3, n)

doc"""
Return the generalized integer Bernoulli numbers ``b_{4}(n) = n``
André``_{4}(n-1)``.

julia> [A281588(n) for n in 0:9]
[0, 1, 2, 3, 4, -5, -24, -98, -272, 621]
"""
A281588(n::Int) = BernoulliInt(4, n)

doc"""
Return a list of length len of the integer Bernoulli numbers ``b_{m}(n)``
using Seidel's boustrophedon algorithm.
"""
function BernoulliIntList(m::Int, len::Int)
    len ≤ 0 && return fmpz[]
    R = zeros(ZZ, len)
    len == 1 && return R
    R[2] = 1
    len == 2 && return R
    A = zeros(ZZ, len)
    A[1] = 1; A[2] = 1

    for n in 1:len - 2
        if n % m ≠ 0
            for i in n:-1:1 A[i] += A[i + 1] end
            C = A[1]
        else
            C = 0
            for i in 1:(n + 2) A[i], C = C, A[i]; C = A[i] - C end
        end
        R[n + 2] = C
    end
    R
end

doc"""
Computes a list of length len of the integer Bernoulli numbers ``b_{2}(n)``
using Seidel's boustrophedon algorithm.

julia> L065619(10)
[0, 1, 2, 3, 8, 25, 96, 427, 2176, 12465]
"""
function L065619(len::Int)
    len ≤ 0  && return fmpz[]
    R = SeqArray(len)
    len == 1 && return R
    len == 2 && (R[1] = 1; return R)

    A = Dict{Int,fmpz}(-1 => 1, 0 => 0)
    k = 0; e = 1
    for i in 0:len - 1
        Am = 0; A[k + e] = 0; e = -e
        for j in 0:i
            Am += A[k]; A[k] = Am; k += e
        end
        j = e < 0 ? div(-i, 2) : div(i, 2)
        R[i] = A[j]
    end
    R
end

doc"""
Return the rational Bernoulli number ``B_n``. Cf. A027641/A027642.

julia> Bernoulli(12)
-691/2730
"""
function Bernoulli(n::Int)
    IsOdd(n) && (n == 1 ? (return fmpq(1, 2)) : (return fmpq(0, 1)))
    n == 0 && return fmpq(1)
    denom = ^(ZZ(4), n) - ^(ZZ(2), n)
    fmpq(BernoulliInt(2, n), denom)
end

doc"""
Return a list of the first len Bernoulli numbers ``B_n``. Cf. A027641/A027642.

julia> BernoulliList(10)
[1, 1//2, 1//6, 0, -1//30, 0, 1//42, 0, -1//30, 0]
"""
function BernoulliList(len::Int)
    if len ≤ 0 return fmpq[] end
    R = Array{fmpq}(len)
    R[1] = fmpq(1, 1); len == 1 && return R
    R[2] = fmpq(1, 2); len == 2 && return R

    A = Dict{Int,fmpz}(0 => 1, -2 => 0, -1 => 1, 1 => 0)
    a = fmpz(12); b = fmpz(240)
    k = e = 1

    for i in 2:len - 1
        Am = 0; A[k + e] = 0; e = -e
        for j in 0:i
            Am += A[k]; A[k] = Am; k += e
        end
        if e > 0
            R[i + 1] = fmpq(0, 1)
        else
            d = i >> 1
            R[i + 1] = IsEven(d) ? fmpq(-A[-d], a) : fmpq(A[-d], a)
            a, b = b, b << 4 + b << 2 - a << 6
        end
    end
    R
end

doc"""
Return the numerator of the Bernoulli number ``B_n``.

julia> [A027641(n) for n in 0:10]
[1, -1, 1, 0, -1, 0, 1, 0, -1, 0, 5]
"""
function A027641(n::Int)
    IsOdd(n) && (n == 1 ? (return ZZ(-1)) : return ZZ(0))
    n == 0 && return ZZ(1)
    denom = ^(ZZ(4), n) - ^(ZZ(2), n)
    Nemo.numerator(BernoulliInt(2, n) // denom)
end

# We could also define the denominator of the Bernoulli number as
#    A027642(n::Int) = denominator(Bernoulli(n))
# however A027642 is implemented more efficiently via the
# Clausen numbers (see the module "Clausen").

doc"""
Return denominator(Bernoulli ``_{n+1}(x) - `` Bernoulli ``_{n+1})``.

julia> [A195441(n) for n in 0:9]
[1, 1, 2, 1, 6, 2, 6, 3, 10, 2]
"""
function A195441(n::Int)
    n < 4 && return ZZ([1, 1, 2, 1][n + 1])
    P = Primes(2, div(n + 2, 2 + n % 2))
    ∏([p for p in P if p ≤ sum(digits(n + 1, Int(p)))])
end

end # module

module BernoulliNumbersTest
using Base.Test, SeqBase, Nemo, NumberTheory, Clausen, BernoulliNumbers
using SeqTests, OEISUtils

function test()

    @testset "BernoulliNum" begin
        @test BernoulliInt(2, 10) == 79360
        @test BernoulliInt(3, 30) == -7708110416280010548302670
        @test BernoulliInt(4, 40) == -44494882577309421077208834962882560

        @test isa(BernoulliInt(3, 30), fmpz)
        @test isa(BernoulliIntList(2, 20)[end], fmpz)
        @test isa(Bernoulli(0), fmpq)

        @test isa(A027641(10), fmpz)
        @test isa(A195441(10), fmpz)

        for m in 1:20
            @test BernoulliInt(m, 200) == BernoulliIntList(m, 200 + 1)[end]
        end

        for n in 0:20
            @test denominator(Rational(Bernoulli(2 * n))) == ClausenNumber(n)
        end

        t = fmpz(8622529719094842064796322984685715031642180319435676189471082876882585178647210)
        @test A195441(10000) == t

        # A065619 in the OEIS has an arbitrary offset of 1.
        l = SeqArray([0, 1, 2, 3, 8, 25, 96, 427, 2176])
        @test all(L065619(9) .== l)

        if oeis_isinstalled()
            A = [A195441, A281586, A281588, A027641]
            SeqTest(A, 'A')
        end
    end
end

function demo()
    for m in 1:5
        println([BernoulliInt(m, n) for n in 0:10])
    end
    for len in 1:6
        println(BernoulliIntList(3, len))
    end
    for n in 0:20
        println(fmpq(Nemo.numerator(Bernoulli(n)), Nemo.denominator(Bernoulli(n))))
    end
end

doc"""
Since the BernoulliInts are cached via the André numbers the benchmarks below
depend on the cache history.

for m in 1:10, n in 0:499 BernoulliInt(m,n) end :: 1.269732 seconds (2.02 M allocations: 46.747 MB, 8.12% gc time)
for m in 1:10 BernoulliIntList(m,500) end :: 1.217272 seconds (1.27 M allocations: 19.564 MB, 7.44% gc time)
for n in 0:1000 Bernoulli(n) end :: 1.506394 seconds (1.13 M allocations: 25.936 MB, 5.55% gc time)
BernoulliList(1000) :: 0.647778 seconds (505.27 k allocations: 7.823 MB)
for n in 0:1000 A027641(n) end :: 0.003213 seconds (5.24 k allocations: 97.422 KB)
for n in 0:10000 A195441(n) end :: 1.051287 seconds (3.69 M allocations: 344.337 MB, 13.38% gc time)
"""
function perf()
    gc()
    for m in 1:2, n in 0:1 BernoulliInt(m, n) end
    @time (for m in 1:10, n in 0:499 BernoulliInt(m, n) end)

    for m in 1:2, len in 1:2 BernoulliIntList(m, len) end
    @time (for m in 1:10 BernoulliIntList(m, 500) end)

    gc()
    for n in 0:10 Bernoulli(n) end
    @time (for n in 0:1000 Bernoulli(n) end)

    BernoulliList(10)
    @time BernoulliList(1000)

    for n in 0:10 A027641(n) end
    @time (for n in 0:1000 A027641(n) end)

    gc()
    A195441(10)
    @time (for n in 0:10000 A195441(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module

# [0] 0, 1,  2, 3,  4,  5,    6,    7,     8,     9,     10
# [1] 0, 1, -2, 3, -4,  5,   -6,    7,    -8,     9,    -10
# [2] 0, 1, 2, -3, -8,  25,  96, -427, -2176, 12465,  79360
# [3] 0, 1, 2,  3, -4, -15, -54,  133,   792,  4293, -15130
# [4] 0, 1, 2,  3,  4,  -5, -24,  -98,  -272,   621,   4960
# [5] 0, 1, 2,  3,  4,   5,  -6,  -35,  -160,  -495,  -1250
