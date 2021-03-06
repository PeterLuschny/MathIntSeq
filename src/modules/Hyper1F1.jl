# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module Hyper1F1
using Nemo

export GammaHyp, A000262, A007060, A033815, A099022, A251568

# Numerical evaluation based on hypergeometric functions.
# Nemo.hyp1f1(a::acb, b::acb, x::acb)
# May fail if the required precicison is > 10000.

doc"""
Return ``Γ(a) 1F1(b, c, d).``
"""
function GammaHyp(a, b, c, d)
    prec = 500
    while prec <= 10000
        CC = AcbField(prec)
        c = gamma(CC(a)) * hyp1f1(CC(b), CC(c), CC(d))
        b, i = unique_integer(c)
        b && return i
        prec *= 2
    end
    println("$a $b $c $d gives an InexactError!")
    # throw(InexactError()) error()
end

doc"""
Return ``n!`` Hyper``1F1[1-n, 2, -1]``.
Number of partitions of ``{1,...,n}`` into any number of ordered subsets.
"""
A000262(n::Int) = n == 0 ? ZZ(1) : GammaHyp(n + 1, 1 - n, 2, -1)

doc"""
Return ``(2n)!`` Hyper``1F1[-n, -2n, -2]``.
Number of ways ``n`` couples can sit in a row without any spouses next to each other.
"""
A007060(n::Int) = n == 0 ? ZZ(1) : GammaHyp(2n + 1, -n, -2n, -2)

doc"""
Return ``(2n)!`` Hyper``1F1[-n, -2n, -1]``.
Number of acyclic orientations of the Turán graph ``T(2n,n)``. (Alois P. Heinz)
"""
A033815(n::Int) = n == 0 ? ZZ(1) : GammaHyp(2n + 1, -n, -2n, -1)

doc"""
Return ``(2n)!`` Hyper``1F1[-n, -2n, 1]``.
``\sum_{k=0..n} \binom{n}{k} (2n-k)!``.
"""
A099022(n::Int) = n == 0 ? ZZ(1) : GammaHyp(2n + 1, -n, -2n, 1)

doc"""
Return ``((2n)!/(n+1)!)`` Hyper``1F1[1-n, n+2, -1]``.
Egf. ``exp(x C(x)^2)`` where ``C(x) = 1 + xC(x)^2`` is the gf. of the Catalan numbers.
"""
function A251568(n::Int)
    n == 0 && return fmpz(1)
    prec = 500
    while prec <= 10000
        CC = AcbField(prec)
        c = gamma(CC(2 * n + 1)) * hyp1f1(CC(1 - n), CC(n + 2), CC(-1)) / gamma(CC(n + 2))
        b, i = unique_integer(c)
        b && return i
        prec *= 2
    end
    println("n = $n gives an InexactError!")
    # throw(InexactError()) error()
end

end # module

module Hyper1F1Test
using Base.Test, SeqBase, SeqTests, Nemo, OEISUtils, Hyper1F1

function test()
    @testset "Hyper1F1" begin

        for n in [10, 73, 150]
            @test A000262(n) ≠ 0
            @test A007060(n) ≠ 0
            @test A033815(n) ≠ 0
            @test A099022(n) ≠ 0
            @test A251568(n) ≠ 0
        end

        if oeis_isinstalled()

            A = [A000262, A007060, A033815, A099022, A251568]
            SeqTest(A, 'A')
        end
    end
end

function demo()
end

doc"""
(for n in 0:150 A000262(n) end) :: 0.010073 seconds (3.75 k allocations: 166.422 KB)
"""
function perf()
    gc()
    @time (for n in 0:150 A000262(n) end)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module
