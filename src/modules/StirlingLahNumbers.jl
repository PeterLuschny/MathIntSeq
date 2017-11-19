# This file is part of OLMS (Open Library of Mathematical Sequences).
# Copyright Peter Luschny. License is MIT.

module StirlingLahNumbers
using Memoize, Nemo, SeqBase

export Lah, LahTriangle, StirlingLah, StirlingLahTriangle
export StirlingSet, StirlingSetTriangle, StirlingCycle, StirlingCycleTriangle
export T038207, T048993, T132393, T271703, T269945, T269944, T268434, T269948
export T269947, T269946

#     The hierarchy of the Stirling/Lah numbers:
#     Order	Stirling set Stirling cycle Lah numbers
#      0      A007318       A007318      A038207
#      1      A048993       A132393      A271703
#      2      A269945       A269944      A268434
#      3      A269948       A269947      A269946

# Note that many OEIS triangle are (1,1)-based wheras we always assume a triangle
# to be (0,0)-based (which allows the canonical interpretation as coefficients
# of a polynomial sequence). This implies that references to the OEIS are in
# these cases hints, not identities (for example T038207).

doc"""
Compute the generalized Stirling-Lah recurrence with weight function w.
"""
function StirlingLah(n::Int, k::Int, w::Function)

    @memoize function T(n::Int, k::Int)
        if n == k return fmpz(1) end
        if k < 0 || k > n return fmpz(0) end
        T(n - 1, k - 1) + w(n, k) * T(n - 1, k)
    end
    T(n, k)
end

doc"""
Return the Stirling set numbers (a.k.a. Stirling numbers of 2nd kind)
of order ``m``. Case ``m = 1`` is the classical case.
"""
StirlingSet(n::Int, k::Int, m = 1) = StirlingLah(n, k, (n, k) -> fmpz(k^m))

doc"""
Return the Stirling cycle numbers (a.k.a. unsigned Stirling numbers of 1st kind)
of order ``m``. Case ``m = 1`` is the classical case.
"""
StirlingCycle(n::Int, k::Int, m = 1) = StirlingLah(n, k, (n, k) -> fmpz((n - 1)^m))

doc"""
Return the Lah numbers of order ``m``. Case ``m=1`` is the classical case.
"""
Lah(n::Int, k::Int, m = 1) = StirlingLah(n, k, (n, k) -> k^m + fmpz((n - 1)^m))

doc"""
Compute the generalized Stirling-Lah square matrix of dimension dim and weight
function w.
"""
function StirlingLahMatrix(dim::Int, w::Function)
    M = SeqMatrix(dim)

    for n in 0:dim - 1
        M[n, n] = fmpz(1)
        for k in 1:n - 1
            M[n, k] = w(n - 1, k - 1) * M[n - 1, k] + M[n - 1, k - 1]
        end
    end
    M
end

doc"""
Compute the generalized Stirling-Lah triangle of dimension dim and weight
function w.
"""
function StirlingLahTriangle(dim::Int, w::Function)
    T = SeqTriangle(dim)

    T[0] = h = 1
    for s in 2:dim
        l = h; h += s
        T[l] = 0; T[h - 1] = 1
        for k in h - 2:-1:l + 1
            T[k] = w(s - 1, k - l) * T[k - s + 1] + T[k - s]
        end
    end
    T
end

doc"""
Return the first dim rows of the triangle of Stirling set numbers of order ``m``.
Case ``m=1`` is the classical case.
"""
StirlingSetTriangle(dim::Int, m = 1) = StirlingLahTriangle(dim, (n, k) -> fmpz(k^m))

# doc"""
# Pascal's triangle
# """
# T007318(dim::Int) = StirlingSetTriangle(dim, 0)

doc"""
Return the first dim rows of the triangle of Stirling set numbers.
Counts partitions of an ``n``-set into ``k`` nonempty subsets.
"""
T048993(dim::Int) = StirlingSetTriangle(dim, 1)

doc"""
Return the first dim rows of the triangle of Stirling set numbers of order 2,
also known as central factorial numbers ``T(2n, 2k)``.
"""
T269945(dim::Int) = StirlingSetTriangle(dim, 2)

doc"""
Triangle of Stirling set numbers of order 3, also called 3rd central factorial
numbers.
"""
T269948(dim::Int) = StirlingSetTriangle(dim, 3)

doc"""
Return the first dim rows of the triangle of Stirling cycle numbers of order ``m``.
Case ``m=1`` is the classical case.
"""
StirlingCycleTriangle(dim::Int, m = 1) = StirlingLahTriangle(dim, (n, k) -> fmpz((n - 1)^m))

doc"""
Return the first dim rows of the triangle of Stirling cycle numbers.
Counts permutations of ``n`` objects with exactly ``k`` cycles.
"""
T132393(dim::Int) = StirlingCycleTriangle(dim, 1)

doc"""
Return the first dim rows of the triangle of Stirling cycle numbers of order 2,
also known as central factorial numbers ``|t(2n, 2k)|``.
"""
T269944(dim::Int) = StirlingCycleTriangle(dim, 2)

doc"""
Return the first dim rows of the triangle of Stirling cycle numbers of order 3.
"""
T269947(dim::Int) = StirlingCycleTriangle(dim, 3)

doc"""
Return the first dim rows of the triangle of Lah numbers of order ``m``.
Case ``m=1`` is the classical case.
"""
LahTriangle(dim::Int, m = 1) = StirlingLahTriangle(dim, (n, k) -> k^m + fmpz((n - 1)^m))

doc"""
Return the first dim rows of the triangle of the square of the Pascal triangle.
"""
T038207(dim::Int) = LahTriangle(dim, 0)

doc"""
Return the first dim rows of the triangle of unsigned Lah numbers.
Counts number of partitions of ``{1..n}`` into ``k`` ordered subset.
"""
T271703(dim::Int) = LahTriangle(dim, 1)

doc"""
Return the first dim rows of the triangle of Lah numbers of order 2.
"""
T268434(dim::Int) = LahTriangle(dim, 2)

doc"""
Return the first dim rows of the triangle of Lah numbers of order 3.
"""
T269946(dim::Int) = LahTriangle(dim, 3)

end # module

module StirlingLahNumbersTest
using Base.Test, SeqBase, SeqTests, Nemo, OEISUtils, StirlingLahNumbers

function test()
    @testset "Stirling-Lah" begin

        SST = StirlingSetTriangle(6, 2)
        a = SeqArray([0, 1, 85, 147, 30, 1])
        b = Row(SST, 5)
        @test all(a .== b)

        SCT = StirlingCycleTriangle(6, 2)
        a = SeqArray([0, 576, 820, 273, 30, 1])
        b = Row(SCT, 5)
        @test all(a .== b)

        LT = LahTriangle(6, 2)
        a = SeqArray([0, 1700, 2900, 840, 60, 1])
        b = Row(LT, 5)
        @test all(a .== b)

        if oeis_isinstalled()

            SLT = [T048993, T132393, T271703, T269945, T269944, T268434, T269948,
                T269947, T269946]
            SeqTest(SLT, 'T')
        end
    end
end

function demo()
    println(); println("generalized StirlingSet numbers")
    for m in 0:3
        SST = StirlingSetTriangle(6, m)
        Show(SST)
    end
    println(); println("generalized StirlingCycle numbers")
    for m in 0:3
        SCT = StirlingCycleTriangle(6, m)
        Show(SCT)
    end
    println(); println("generalized Lah numbers")
    for m in 0:3
        LT = LahTriangle(6, m)
        Show(LT)
    end
    println()
end

doc"""
StirlingSetTriangle(100, 2) :: 1.757261 seconds (1.50 M allocations: 26.669 MB, 7.55% gc time)
LahTriangle(1000, 2) :: 1.968734 seconds (2.00 M allocations: 34.275 MB, 8.35% gc time)
"""
function perf()
    gc()
    @time StirlingSetTriangle(1000, 2)
    gc()
    @time LahTriangle(1000, 2)
end

function main()
    test()
    demo()
    perf()
end

main()

end # module
