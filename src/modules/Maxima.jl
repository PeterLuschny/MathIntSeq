# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Maxima
using NumberTheory, SeqBase

export FindMaxima, FindMaximaUpTo
export B002183, B002182, B034885, B006093, B006005, B002093
export L002183, L002182, L034885, L006093, L006005

doc"""
Return maxima relative to the predecessors and their positions of a sequence f
as a pair of lists (r, i) where r contains the maxima and i the index where
they occur.

f(0) is a maximum and f(i) is a maximum if f(j) < f(i) for all 0 < j < i.
All maxima ≤ bound are returned. In the OEIS such values are often called
'record values'.
"""
function FindMaximaUpTo(bound::Int, f::Function)

    bound == 0 && return [[], []]
    max = f(0); r = [max]; i = [0]
    bound == 1 && return [r, i]

    for n in 1:bound
        if f(n) > max
            max = f(n)
            push!(r, max)
            push!(i, n)
        end
    end
    [SeqArray(r), SeqArray(i)]
end

doc"""
Return maxima relative to the predecessors and their positions of a sequence f
as a pair of lists (r, i) where r contains the maxima and i the index where
they occur.

f(0) is a maximum and f(i) is a maximum if f(j) < f(i) for all 0 < j < i.
In the OEIS such values are often called 'record values'.
"""
function FindMaxima(len::Int, f::Function)

    len == 0 && return [[], []]
    max = f(0); r = [max]; i = [0]
    len == 1 && return [r, i]

    wm, n = 1, 1
    while wm < len
        if f(n) > max
            max = f(n)
            push!(r, max)
            push!(i, n)
            wm += 1
        end
        n += 1
    end
    [SeqArray(r), SeqArray(i)]
end

doc"""
Return the number of divisors of n-th highly composite numbers
which do not exceed the given bound.
"""
B002183(bound) = FindMaximaUpTo(bound, τ)[1]

doc"""
Return the number of divisors of n-th highly composite numbers.
"""
L002183(len) = FindMaxima(len, τ)[1]

doc"""
Return highly composite numbers which do not exceed the given bound.
"""
B002182(bound) = FindMaximaUpTo(bound, τ)[2]

doc"""
Return highly composite numbers.
"""
L002182(len) = FindMaxima(len, τ)[2]

doc"""
Return the record values of sigma(n) which do not exceed the given bound.
"""
B034885(bound) = FindMaximaUpTo(bound, σ)[1]

doc"""
Return the record values of sigma(n).
"""
L034885(len) = FindMaxima(len, σ)[1]

doc"""
Return highly abundant numbers which do not exceed the given bound.
"""
B002093(bound) = FindMaximaUpTo(bound, σ)[2]

# Version in Abundant.jl is more efficient.
#doc"""
#Return a list of length len of highly abundant numbers.
#"""
#L002093(len) = FindMaxima(len, σ)[2]

doc"""
Return the primes minus 1 which do not exceed the given bound.
"""
B006093(bound) = FindMaximaUpTo(bound, ϕ)[1]

doc"""
Return the primes minus 1.
"""
L006093(len) = FindMaxima(len, ϕ)[1]

doc"""
The odd prime numbers together with 1 which do not exceed the given bound.
"""
B006005(bound) = FindMaximaUpTo(bound, ϕ)[2]

doc"""
The odd prime numbers together with 1.
"""
L006005(len) = FindMaxima(len, ϕ)[2]

end # module

module MaximaTest
using Base.Test, Maxima, NumberTheory, SeqBase, OEISUtils

function test()

    @testset "Maxima" begin

        if oeis_isinstalled()

            B = [B002183, B002182, B034885, B006093, B006005 ]
            L = [L002183, L002182, L034885, L006093, L006005 ]

            for seq in B
                name = SeqName(seq)
                O = oeis_local(name, 12)
                # the parameter is 'search bound'.
                S = seq(300)
                # the OEIS offset is 1, therefore
                @test all(S[1:10] .== O[0:9])
            end

            for seq in L
                name = SeqName(seq)
                O = oeis_local(name, 12)
                # the parameter is 'length'
                S = seq(12)
                # the OEIS offset is 1, therefore
                @test all(S[1:10] .== O[0:9])
            end
        end
    end
end

function demo()
    println(FindMaximaUpTo(200, τ))
    println()
    println(FindMaxima(20, τ)[1])
    println(FindMaxima(20, τ)[2])
    println(FindMaximaUpTo(20, τ)[1])
    println(FindMaximaUpTo(20, τ)[2])
end

function perf()
end

function main()
    test()
    demo()
    perf()
end

main()

end # module
