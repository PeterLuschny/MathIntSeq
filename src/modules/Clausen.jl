# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Clausen
using Nemo, SeqBase, PrimeSieve, Products, NumberTheory

export ClausenNumber, ClausenNumberList
export A002445, L002445, A027642

doc"""
Return the Clausen number ``C_n`` which is the denominator of the Bernoulli
number ``B_{2n}``.
"""
function ClausenNumber(n::Int)
    n == 0 && return ZZ(1)
    m = [d + 1 for d in Divisors(2n)]
    ∏([q for q in m if IsPrime(q)])
end

doc"""
Return the list of length len of Clausen numbers which are the denominators of
the Bernoulli numbers ``B_{2n}``.
"""
function ClausenNumberList(len::Int)
    len ≤ 0 && return fmpz[]
    A = FillArray(2, len); A[0] = 1
    m = len - 1
    m == 0 && return A

    for p in Primes(3, 2m + 1)
        r = Int(div(p - 1, 2))
        for k in range(r, r, div(m, r))
            A[k] *= p
        end
    end
    A
end

doc"""
Return the Clausen number ``C(n)`` which is the denominator of the Bernoulli
number ``B_{2n}``.
"""
A002445(n::Int) = ClausenNumber(n)

doc"""
Return the list of length len of Clausen numbers which are the denominators
of the Bernoulli numbers ``B_{2n}``.
"""
L002445(len::Int) = ClausenNumberList(len)

doc"""
Return the denominator of Bernoulli number ``B_n``.
"""
function A027642(n::Int)
    IsEven(n) && return ClausenNumber(div(n, 2))
    n == 1 && return ZZ(2)
    return ZZ(1)
end

end # module

module ClausenTest
using Base.Test, SeqBase, SeqTests, OEISUtils, Nemo, Clausen

function test()

    @testset "Clausen" begin
        C = ClausenNumberList(800)
        @test C[124] == 30
        @test C[780] == 32695402455500348373810
        @test C[793] == 6

        @test isa(C[781], fmpz)
        @test isa(ClausenNumber(10), fmpz)

        if oeis_isinstalled()
            SeqTest([L002445], 'L')
            SeqTest([A002445, A027642], 'A')
        end
    end
end

function demo()
    SeqShow(ClausenNumberList(7))
    for n in 0:6
        println(n, " ↦ ", ClausenNumberList(n))
    end
    for n in 0:10
        println(n, " ↦ ",  A027642(n))
    end
end

doc"""
ClausenNumberList(10000) :: 0.015402 seconds (141.11 k allocations: 2.344 MB)
for n in 0:10000 ClausenNumber(n) end :: 1.008843 seconds (798.01 k allocations: 40.378 MB, 6.15% gc time)
"""
function perf()
    # ClausenNumberList is ~ 60 times faster than a list of ClausenNumbers!
    gc()
    @time ClausenNumberList(10000)
    gc()
    @time (for n in 0:10000 ClausenNumber(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module
